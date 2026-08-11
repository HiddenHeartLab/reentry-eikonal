from scipy import sparse
from petsc4py import PETSc
from petsc4py.PETSc import ScalarType, Mat
import pyvista as pv
import numpy as np
from tqdm import trange
import meshio
from ioniclib.ioniclib import IonicModel
import time

def Bmatrix(nodeCoords,fl,ft,sigmal,sigmat):
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

    f_l = np.array([np.dot(fl,e1),np.dot(fl,e2)])
    f_t = np.array([np.dot(ft,e1),np.dot(ft,e2)])
    D = sigmal*np.outer(f_l,f_l)+sigmat*np.outer(f_t,f_t)

    return B, D, J

def localStiffnessMatrix(B,D,J):
    "Assemble the local stiffness matrix"
    return (B.T @ D @ B)/(2.*J)
    
def localMassMatrix(J):
    "Assemble the local mass matrix"
    return np.array([[2.0,1.0,1.0],
                     [1.0,2.0,1.0],
                     [1.0,1.0,2.0]])*J/24

def assembleParabolic(pts,elm,fl,ft,sigmal,sigmat):
    "Assemble the global mass and stiffness matrices"

    I,J,Vm,Vk = [],[],[],[]
    for k,tri in enumerate(elm):
        j, i = np.meshgrid(tri,tri)
        I.extend(list(i.ravel()))
        J.extend(list(j.ravel()))
        B,D,Jac = Bmatrix(pts[tri],fl,ft,sigmal,sigmat)
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

# Read the mesh
square = meshio.read('square.o.mesh')
pts = square.points
elm = square.cells_dict['triangle']

# Matrices for the monodomain equation
sigmal = 1.5
sigmat = 0.24
fl = np.array([1,0])
ft = np.array([0,1])
M, K = assembleParabolic(pts,elm,fl,ft,sigmal,sigmat)

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
ksp = PETSc.KSP().create()
ksp.setOperators(A)
ksp.setType('cg')
ksp.setConvergenceHistory()
ksp.getPC().setType('hypre')

# stimulation current
Istim0 = A.createVecRight()
Istim0.array[:] = Imax*(pts[:,0]<=0.1).astype(float)

Istim1 = A.createVecRight()
Istim1.array[:] = Imax*(pts[:,0]<=6)*(pts[:,1]<=12).astype(float)
t1 = 210

# ionic model
ion_model = IonicModel(vm.size,dt=dt)
ion_model.set_parameter('g_CaL',0.037)
ion_model.set_parameter('g_to',0.026)
ion_model.set_parameter('g_K1',0.18)
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
start = time.time()
with meshio.xdmf.TimeSeriesWriter("spiral.xdmf") as writer:
    writer.write_points_cells(pts, [("triangle",elm)])
    for i in trange(ndt):
        # current time
        t = dt * i

        # stimulus currents
        if t<=Idur:
            Is.array[:] = Istim0
        elif (t>=(t1+dt) and (t<=(t1+Idur))):
            Is.array[:] = Istim1
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
