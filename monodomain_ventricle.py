from scipy import sparse
from petsc4py import PETSc
from petsc4py.PETSc import ScalarType, Mat
import pyvista as pv
import numpy as np
from tqdm import trange
import meshio
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

def assembleParabolic(pts,elm,t2s,t2f,sigmah,sigmabz,ratio,pre,block):
    "Assemble the global mass and stiffness matrices"

    I,J,Vm,Vk = [],[],[],[]
    for k,tri in enumerate(elm):
        j, i = np.meshgrid(tri,tri)
        I.extend(list(i.ravel()))
        J.extend(list(j.ravel()))
        fiber = t2f[k,:]
        if t2s[k]==1:
            sigma = 0
        elif t2s[k]==0:
            sigma = sigmah
        else:
            sigma = sigmabz
        if ((pre==1)&(block[k]==1)):
            sigma = 0
        B,D,Jac = Bmatrix(pts[tri],fiber,sigma,ratio)
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
ventricle = meshio.read('meshes/ventricle.vtk')
pts = 0.1*ventricle.points
elm = ventricle.cells_dict['triangle']

# Matrices for the monodomain equation
n_elms = np.shape(elm)[0]
t2s = np.zeros(n_elms)
with open('meshes/t2s.csv') as f:
    reader = csv.reader(f)
    n = 0
    for row in reader:
        t2s[n] = row[3]
        n = n+1
t2f = np.zeros((n_elms,3))
with open('meshes/t2f_ventricle.csv') as f:
    reader = csv.reader(f)
    n = 0
    for row in reader:
        t2f[n,0] = row[3]
        t2f[n,1] = row[4]
        t2f[n,2] = row[5]
        n = n+1
sigmah = 0.46
sigmabz = 0.19
ratio = 1/9
centroids = np.zeros((n_elms,3))
for n in range(n_elms):
    centroids[n,:] = (pts[elm[n,0],:]+pts[elm[n,1],:]+pts[elm[n,2],:])/3
line = np.zeros((101,3))
for dx in np.arange(0,1.01,0.01):
    line[int(dx*100),:] = centroids[291818,:]+dx*(centroids[142635,:]-centroids[291818,:])
distance = np.zeros(n_elms)
for n in range(n_elms):
    distance[n] = np.min(np.linalg.norm(centroids[n,:]-line,axis=1))
block = (distance<=0.1)
M, K = assembleParabolic(pts,elm,t2s,t2f,sigmah,sigmabz,ratio,1,block)

# parameters
Imax = 40.0   # [uA/cm^2]
Idur = 1.0    # [ms]
Cm   = 1.0    # [uF/cm^2]
beta = 800.0  # [1/cm]

# full matrix for backward Euler
dt = 0.02
Tend = 2000
ndt = int(np.rint(Tend/dt)+1)

A = Cm * M + dt/beta*K

# create the solution vector
vm = A.createVecRight()

# create the solver
ksp = solver(A)

# stimulation current
s = pts[69948,:]
r = 0.1
Istim = A.createVecRight()
Istim = Imax*((pts[:,0]-s[0])**2+(pts[:,1]-s[1])**2+(pts[:,2]-s[2])**2<=r**2).astype(float)

# tissue types
n_pts = np.shape(pts)[0]
p2s = np.zeros(n_pts)
with open('meshes/p2s.csv') as f:
    reader = csv.reader(f)
    n = 0
    for row in reader:
        p2s[n] = row[3]
        n = n+1
inds_h = np.argwhere(p2s==0)[:,0]
inds_bz = np.argwhere(p2s!=0)[:,0]

# ionic models
h_min  = 0.3
tau_out = 1
tau_in = h_min*tau_out/4
v_gate = 0.13
v_min = -80
v_max = 0
tau_open_h = 300
tau_close_h = 150
tau_open_bz = 100
tau_close_bz = 250

# initialization
vm.array[:] = v_min*np.ones(vm.size)
h = np.ones(vm.size)
Is    = A.createVecRight()
Iion = A.createVecRight()

# loop
t_unblock = 330
start = time.time()
with meshio.xdmf.TimeSeriesWriter("output/monodomain_ventricle.xdmf") as writer:
    writer.write_points_cells(pts, [("triangle",elm)])
    for i in trange(ndt):
        # current time
        t = dt * i

        # remove the block line
        if t==t_unblock:
            M, K = assembleParabolic(pts,elm,t2s,t2f,sigmah,sigmabz,ratio,0,block)
            A = Cm * M + dt/beta*K
            ksp = solver(A)

        # stimulus currents
        if t<=Idur:
            Is.array[:] = Istim
        else:
            Is.array[:] = 0.0
        
        # update gating variables and get total ion current
        vm_h = vm.array[inds_h]
        h_h = h[inds_h]
        v_h = (vm_h-v_min)/(v_max-v_min)
        h_h = (v_h<v_gate)*(h_h+dt*(1-h_h)/tau_open_h)+(v_h>=v_gate)*(h_h-dt*h_h/tau_close_h)
        h[inds_h] = h_h
        Iion.array[inds_h] = (-h_h*v_h*v_h*(1-v_h)/tau_in+v_h/tau_out)*(v_max-v_min)

        vm_bz = vm.array[inds_bz]
        h_bz = h[inds_bz]
        v_bz = (vm_bz-v_min)/(v_max-v_min)
        h_bz = (v_bz<v_gate)*(h_bz+dt*(1-h_bz)/tau_open_bz)+(v_bz>=v_gate)*(h_bz-dt*h_bz/tau_close_bz)
        h[inds_bz] = h_bz
        Iion.array[inds_bz] = (-h_bz*v_bz*v_bz*(1-v_bz)/tau_in+v_bz/tau_out)*(v_max-v_min)

        # diffusion (implicit)
        b = M*(Cm*vm + dt*(Is - Iion))
        ksp.solve(b, vm)

        # write solution every 50 steps
        if i % 50 == 0:
            v = 1*(vm.array>=-62)-1*(p2s==1)
            writer.write_data(t, point_data={"v": v})
comp_time = time.time()-start
print(comp_time)
