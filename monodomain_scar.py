from petsc4py import PETSc
from petsc4py.PETSc import ScalarType, Mat
import pyvista as pv
import numpy as np
from tqdm import trange
import meshio
import time
from scipy import sparse

def Bmatrix(nodeCoords,fib,sigma,ratio):
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

    f = np.array([np.dot(fib,e1),np.dot(fib,e2)])
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

def assembleParabolic(pts,elm,f,sigmah,sigmabz,ratio):
    "Assemble the global mass and stiffness matrices"

    I,J,Vm,Vk = [],[],[],[]
    for k,tri in enumerate(elm):
        j, i = np.meshgrid(tri,tri)
        I.extend(list(i.ravel()))
        J.extend(list(j.ravel()))

        center = np.sum(pts[tri],axis=0)/3
        if ((center[0]>4)*(center[0]<7)*(center[1]>5)*(center[1]<10)):
            sigma = 0
        elif ((center[0]>8)*(center[0]<11)*(center[1]>5)*(center[1]<10)):
            sigma = 0
        elif ((center[0]>=7)*(center[0]<=8)*(center[1]>5)*(center[1]<10)):
            sigma = sigmabz
        else:
            sigma = sigmah
        B,D,Jac = Bmatrix(pts[tri],f,sigma,ratio)

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
square = meshio.read('meshes/square.o.mesh')
pts = square.points
elm = square.cells_dict['triangle']

# Matrices for the monodomain equation
sigmah = 0.46
sigmabz = 0.19
ratio = 1/9
f = np.array([1,0])
M, K = assembleParabolic(pts,elm,f,sigmah,sigmabz,ratio)

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
ksp = PETSc.KSP().create()
ksp.setOperators(A)
ksp.setType('cg')
ksp.setConvergenceHistory()
ksp.getPC().setType('hypre')

# stimulation current
Istim = A.createVecRight()
Istim.array[:] = Imax*((pts[:,0]>=7)*(pts[:,0]<=8)*(pts[:,1]>=10)*(pts[:,1]<=11)).astype(float)
t1 = 290

# tissue types
n_pts = np.shape(pts)[0]
marker = np.zeros(n_pts)
for n in range(n_pts):
    if ((pts[n,0]>4)*(pts[n,0]<11)*(pts[n,1]>5)*(pts[n,1]<10)):
        marker[n] = 1
    else:
        marker[n] = 0
inds_h = np.argwhere(marker==0)[:,0]
inds_bz = np.argwhere(marker==1)[:,0]

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
Iion  = A.createVecRight()

# loop
start = time.time()
with meshio.xdmf.TimeSeriesWriter("output/monodomain_scar.xdmf") as writer:
    writer.write_points_cells(pts, [("triangle",elm)])
    for i in trange(ndt):
        # current time
        t = dt * i

        # stimulus currents
        if t<=Idur:
            Is.array[:] = Istim
        elif ((t>=t1) and (t<=(t1+Idur))):
            Is.array[:] = Istim
        else:
            Is.array[:] = 0.0

        # update gating variable and get total ion current
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
            v = 1*(vm.array>=-62)-1*((pts[:,0]>=4)*(pts[:,0]<=7)*(pts[:,1]>=5)*(pts[:,1]<=10))-1*((pts[:,0]>=8)*(pts[:,0]<=11)*(pts[:,1]>=5)*(pts[:,1]<=10))
            writer.write_data(t, point_data={"v": v})
comp_time = time.time()-start
print(comp_time)