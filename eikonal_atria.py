import numpy as np
import meshio
from collections import defaultdict
from skfmm import heap
import csv
import time

# read the mesh
m = meshio.read('meshes/atria_aniso.vtk')
points = 0.1*m.points
n_points = np.shape(points)[0]
elems = m.cells_dict['triangle']
n_elems = np.shape(elems)[0]
# p2t = point to triangle
# p2e = point to opposite edges
# p2p = point to adjacent points
p2t = defaultdict(list)
p2e = defaultdict(list)
p2p = defaultdict(set)
for ti,t in enumerate(m.cells_dict['triangle']):
    for p in t:
        p2t[p].append(ti)
        e = tuple(pj for pj in t if pj != p)
        p2e[p].append(e)
        p2p[p].update(e)
max_neighs = max([len(p2p[i]) for i in range(n_points)])
max_elems = max([len(p2t[i]) for i in range(n_points)])

def solve_quadratic_v(M, U):
    "M n x d x d, U n x d"
    ones = np.ones_like(U)
    a = np.einsum('eij,ei,ej->e',M,ones,ones)
    b = np.einsum('eij,ei,ej->e',M,U,ones)
    c = np.einsum('eij,ei,ej->e',M,U,U) - 1
    delta = b**2 - a*c
    valid = (delta >= 0)
    U_z = (b[valid] + np.sqrt(delta[valid]))/a[valid]
    v = np.einsum('eij,ei->ej',M[valid,...],ones[valid,:] * U_z[:,np.newaxis] - U[valid,:])
    lmbda = v / np.sum(v,axis=1,keepdims=True)
    valid2 = np.all(lmbda > 0, axis=1)
    return U_z[valid2],lmbda[valid2,:].T


def hopf_lax_update(p,phi,cv):
    # solve faces
    neigh = np.array(list(p2p[p]))
    edg = np.array(p2e[p])
    inds = np.array([np.argwhere(neigh==edg[i,j]) for (i,j) in [(i,j) for i in range(np.shape(edg)[0]) for j in range(np.shape(edg)[1])]])
    inds = np.reshape(inds,np.shape(edg))
    UU = phi[inds]
    PP = points[edg,:] - points[p,:]
    F = cv[range(len(p2t[p]))]**2
    # DD is inverse metric
    DD = D[p,range(len(p2t[p])),:,:]/F[:,np.newaxis,np.newaxis]
    MM = np.einsum('ijp,ipq,ikq->ijk',PP,DD,PP)
    MMinv = np.linalg.inv(MM)
    Uz_faces,lmbda = solve_quadratic_v(MMinv, UU)

    # solve edges
    Uz_edges = (np.sqrt(np.einsum('eii->ei',MM)) + UU).min()

    Uz = min(Uz_edges,Uz_faces.min()) if len(Uz_faces) > 0 else Uz_edges
    return Uz

# restitution
di_min = 35
APD = np.load('rest/APD.npz')
DIres = APD['DI']
APDres = APD['APD']
CV = np.load('rest/CV.npz')
DIres = CV['DI']
CVres = CV['CV']

# fibers
t2f = np.zeros((n_elems,3))
with open('meshes/t2f_atria_aniso.csv') as f:
    reader = csv.reader(f)
    n = 0
    for row in reader:
        t2f[n,0] = row[3]
        t2f[n,1] = row[4]
        t2f[n,2] = row[5]
        n = n+1

# anisotropy tensor
D = np.zeros((n_points,max_elems,3,3))
for p in range(n_points):
    i = 0
    for t in p2t[p]:
        vl = 1
        ratio = 0.4
        f = t2f[t,:]
        M = (ratio*vl)**2*np.eye(3) + (1-ratio**2)*(vl**2)*np.outer(f,f)
        D[p,i,:,:] = np.linalg.inv(M)
        i = i+1

# block line
x_block = np.argwhere(points[:,0]>0.7)
x_block = np.array(x_block[:,0],dtype=int)
t_unblock = 125

# stimulus
x_stim = np.argwhere((points[:,0]>0.6)&(points[:,0]<0.7)&(points[:,1]>-0.5)&(points[:,2]>2.4))
x_stim = np.array(x_stim[:,0],dtype=int)

# initialization
phi = np.inf*np.ones(n_points)
phi_comp = np.inf*np.ones((n_points,max_neighs))
v = np.zeros(n_points)
di = np.inf*np.ones(n_points)
apd = np.zeros(n_points)
eps = 1e-10*np.ones(1)
cv = eps*np.ones((n_points,max_elems))
h = heap(100*n_points)
ids = -np.ones(n_points)

# loop  
dt = 1
t = 0
T = 1000
start = time.time()
with meshio.xdmf.TimeSeriesWriter("output/eikonal_atria.xdmf") as writer:
    writer.write_points_cells(points, [("triangle",elems)])
    writer.write_data(t, point_data={"v": v})
    while t<T:
        # stimulus
        if t==0:
            phi[x_stim] = t
            for x in x_stim:
                apd[x] = np.interp(di[x],DIres,APDres)
            di[x_stim] = 0
            v[x_stim] = 1
            # include the excitable neighbors in the heap
            for x in x_stim:
                for p in p2p[x]:
                    if ((phi[p]==np.inf)&(np.isin(x_block,p).any()==False)):
                        Pp = np.array(list(p2p[p])) 
                        phi_comp[p,np.argwhere(Pp==x)[0][0]] = phi[x]
                        Tx = p2t[x] 
                        Tp = p2t[p]
                        TxTp = np.intersect1d(Tx,Tp)
                        for i in range(len(TxTp)): 
                            cv[p,np.argwhere(Tp==TxTp[i])[0][0]] = np.interp(di[p],DIres,CVres)
                        if ids[p]==-1:
                            ids[p] = h.push(p,hopf_lax_update(p,phi_comp[p,:],cv[p,:]))
                        else:
                            h.update(ids[p],hopf_lax_update(p,phi_comp[p,:],cv[p,:]))
        # remove the block line
        if t==t_unblock:
            x_block = []
        # iterate within one time step
        t_prev = t
        while h.empty()==False:
            if h.peek()>t_prev+dt:
                break
            # accept one node
            xa,ta = h.pop()
            phi[xa] = ta
            ids[xa] = -1
            phi_comp[xa,:] = np.inf 
            cv[xa,:] = eps 
            apd[xa] = np.interp(di[xa],DIres,APDres)
            di[xa] = 0
            v[xa] = 1
            # include the excitable neighbors in the heap
            for p in p2p[xa]:
                if ((phi[p]==np.inf)&(np.isin(x_block,p).any()==False)):
                    Pp = np.array(list(p2p[p]))
                    phi_comp[p,np.argwhere(Pp==xa)[0][0]] = phi[xa]
                    Ta = p2t[xa]
                    Tp = p2t[p]
                    TaTp = np.intersect1d(Ta,Tp)
                    for i in range(len(TaTp)):
                        if cv[p,np.argwhere(Tp==TaTp[i])[0][0]]==eps:
                            cv[p,np.argwhere(Tp==TaTp[i])[0][0]] = np.interp(di[p],DIres,CVres)
                    if ids[p]==-1:
                        ids[p] = h.push(p,hopf_lax_update(p,phi_comp[p,:],cv[p,:]))
                    else:
                        h.update(ids[p],hopf_lax_update(p,phi_comp[p,:],cv[p,:]))
        t = t_prev+dt
        # consider the non-excitable nodes where v=0 and update the diastolic interval
        inds = (phi+apd<=t).nonzero()[0]
        di[inds] = t-(phi[inds]+apd[inds])
        v[inds] = 0
        # consider the excitable nodes and update the diastolic interval
        inds = (phi==np.inf).nonzero()[0]
        di[inds] = di[inds]+dt
        # update the non-excitable nodes with diastolic interval >= dmin to excitable
        inds = ((phi<np.inf)&(di>=di_min)).nonzero()[0]
        phi[inds] = np.inf
        writer.write_data(t, point_data={"v": v})
comp_time = time.time()-start
print(comp_time)
