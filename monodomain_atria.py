from scipy import sparse
from petsc4py import PETSc
from petsc4py.PETSc import ScalarType, Mat
import pyvista as pv
import numpy as np
from tqdm import trange
import meshio
from ioniclib.ioniclib import IonicModel
import csv
import time

def Bmatrix(nodeCoords,fiber,sigma,ratio):
    e1 = (nodeCoords[1,:] - nodeCoords[0,:])/np.linalg.norm(nodeCoords[1,:] - nodeCoords[0,:])
    e2 = ((nodeCoords[2,:] - nodeCoords[0,:]) - np.dot((nodeCoords[2,:] - nodeCoords[0,:]),e1)*e1)
    e2 = e2/np.linalg.norm(e2) # normalize
    
    x21 = np.dot(nodeCoords[1,:] - nodeCoords[0,:],e1)
    x13 = np.dot(nodeCoords[0,:] - nodeCoords[2,:],e1)
    x32 = np.dot(nodeCoords[2,:] - nodeCoords[1,:],e1)
    
    y23 = np.dot(nodeCoords[1,:] - nodeCoords[2,:],e2)
    y31 = np.dot(nodeCoords[2,:] - nodeCoords[0,:],e2)
    y12 = np.dot(nodeCoords[0,:] - nodeCoords[1,:],e2)
    
    J = x13*y23 - y31*x32
    
    B = np.array([[y23, y31, y12],[x32, x13, x21]])

    f = np.array([np.dot(fiber,e1),np.dot(fiber,e2)])
    D = ratio*sigma*np.eye(2) + (1-ratio)*sigma*np.outer(f,f)
    
    return B, D, J

def localStiffnessMatrix(B,D,J):
    "Assemble the local stiffness matrix"
    return (B.T @ D @ B)/(2.*J)
    
def localMassMatrix(J):
    "Assemble the local mass matrix"
    return np.array([[2.0,1.0,1.0],
                     [1.0,2.0,1.0],
                     [1.0,1.0,2.0]])*J/24

def assembleParabolic(pts,elm,pre,t2f,sigma,ratio):
    "Assemble the global mass and stiffness matrices"

    I,J,Vm,Vk = [],[],[],[]
    for k,tri in enumerate(elm):
        j, i = np.meshgrid(tri,tri)
        I.extend(list(i.ravel()))
        J.extend(list(j.ravel()))
        fiber = t2f[k,:]
        center = np.sum(pts[tri],axis=0)/3
        if ((pre==1)&(center[0]>0.7)):
            s = 0
        else:
            s = sigma
        B, D, Jac = Bmatrix(pts[tri],fiber,s,ratio)
        Mloc = localMassMatrix(Jac)
        Kloc = localStiffnessMatrix(B,D,Jac)
        Vm.extend(list(Mloc.ravel()))
        Vk.extend(list(Kloc.ravel()))

    n = pts.shape[0]
    M = sparse.coo_matrix((Vm,(I,J)),shape=(n,n)).tocsr()
    K = sparse.coo_matrix((Vk,(I,J)),shape=(n,n)).tocsr()

    # convert to PETSc
    Kp = PETSc.Mat()
    Kp.createAIJWithArrays((n,n),(K.indptr,K.indices,K.data))
    Kp.assemble()

    Mp = PETSc.Mat()
    Mp.createAIJWithArrays((n,n),(M.indptr,M.indices,M.data))
    Mp.assemble()

    return Mp, Kp

def solver(A):

    ksp = PETSc.KSP().create()
    ksp.setOperators(A)
    ksp.setType('cg')
    ksp.setConvergenceHistory()
    ksp.getPC().setType('hypre')

    return ksp


# Read the mesh
atria = meshio.read('atria.vtk')
pts = 0.1*atria.points
elm = atria.cells_dict['triangle']
n_elms = np.shape(elm)[0]

# Matrices for the monodomain equation
t2f = np.zeros((n_elms,3))
with open('t2f_atria.csv') as f:
    reader = csv.reader(f)
    n = 0
    for row in reader:
        t2f[n,0] = row[3]
        t2f[n,1] = row[4]
        t2f[n,2] = row[5]
        n = n+1
sigma = 1.5
ratio = 0.16
M, K = assembleParabolic(pts,elm,1,t2f,sigma,ratio)

# parameters
Imax = 40.0   # [uA/cm^2]
Idur = 1.0    # [ms]
Cm   = 1.0    # [uF/cm^2]
beta = 800.0  # [1/cm]

# full matrix for backward Euler
dt = 0.02
Tend = 1000
ndt = int(np.rint(Tend/dt)+1)

A = Cm * M + dt/beta*K

# create the solution vector
vm = A.createVecRight()

# create the solver
ksp = solver(A)

# stimulation current
Istim = A.createVecRight()
Istim = Imax*((pts[:,0]>0.6)*(pts[:,0]<0.7)*(pts[:,1]>-0.5)*(pts[:,2]>2.4)).astype(float)

# ionic model
ion_model = IonicModel(vm.size,dt=dt)
ion_model.set_parameter('g_CaL',0.037)
ion_model.set_parameter('g_K1',0.18)
ion_model.set_parameter('g_to',0.026)
ion_model.set_parameter('g_Ks',0.1)
ion_model.set_parameter('g_Na',8.0)

# initialization
Isd   = np.zeros_like(vm.array)
Is    = A.createVecRight()
Iion  = A.createVecRight()
dtime = np.zeros((vm.size,),dtype=np.single)
y = np.zeros((vm.size,25),dtype=np.double)

pre = np.load('prepacing.npz')
vm_init = pre['vm']
vm.array[:] = np.tile(vm_init[0],vm.size)
y_init = pre['y']
y = np.tile(y_init[0,:],(vm.size,1))

# loop
t_unblock = 125
start = time.time()
with meshio.xdmf.TimeSeriesWriter("atria.xdmf") as writer:
    writer.write_points_cells(pts, [("triangle",elm)])
    for i in trange(ndt):
        # current time
        t = dt * i

        # remove the block line
        if t==t_unblock:
            M, K = assembleParabolic(pts,elm,0,t2f,sigma,ratio)
            A = Cm * M + dt/beta*K
            ksp = solver(A)

        # stimulus currents
        if t<=Idur:
            Is.array[:] = Istim
        else:
            Is.array[:] = 0.0

        # update gating variables and get total ion current
        Isd[:] = - Is.array
        ion_model.step(t, vm.array, y, Isd, dtime, Iion.array)

        # diffusion (implicit)
        b = M*(Cm*vm + dt*(Is - Iion))
        ksp.solve(b, vm)

        # write solution every 50 steps
        if i % 50 == 0:
            v = 1*(vm.array>=-62)
            writer.write_data(t, point_data={"v": v})
comp_time = time.time()-start
print(comp_time)
