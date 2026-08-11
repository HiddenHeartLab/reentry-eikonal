import numpy as np
import meshio
from collections import defaultdict
from skfmm import heap
import time
    
# read the mesh
m = meshio.read('meshes/square_anisoMS.mesh')
points = m.points
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
di_min_h = 60
APD = np.load('rest/APD_healthy.npz')
DIres_h = APD['DI']
APDres_h = APD['APD']
CV = np.load('rest/CV_healthy.npz')
DIres_h = CV['DI']
CVres_h = CV['CV']
di_min_bz = 25
APD = np.load('rest/APD_bz.npz')
DIres_bz = APD['DI']
APDres_bz = APD['APD']
CV = np.load('rest/CV_bz.npz')
DIres_bz = CV['DI']
CVres_bz = CV['CV']

# tissue types
p2s = np.zeros(n_points)
bz = np.argwhere((points[:,0]>=7)*(points[:,0]<=8)*(points[:,1]>5)*(points[:,1]<10))[:,0]
p2s[bz] = 0.5
scar1 = np.argwhere((points[:,0]>4)*(points[:,0]<7)*(points[:,1]>5)*(points[:,1]<10))[:,0]
scar2 = np.argwhere((points[:,0]>8)*(points[:,0]<11)*(points[:,1]>5)*(points[:,1]<10))[:,0]
p2s[scar1] = 1
p2s[scar2] = 1

# anisotropy tensor
f = np.array([1,0])
vl = 1
D = np.zeros((n_points,max_elems,2,2))
for p in range(n_points):
    i = 0
    for t in p2t[p]:
        center = np.sum(points[elems[t,:]],axis=0)/3
        if ((center[0]>4)*(center[0]<7)*(center[1]>5)*(center[1]<10)):
            D[p,i,:,:] = 1e10*np.eye(2)
        elif ((center[0]>8)*(center[0]<11)*(center[1]>5)*(center[1]<10)):
            D[p,i,:,:] = 1e10*np.eye(2)
        else:
            if ((center[0]>=7)*(center[0]<=8)*(center[1]>5)*(center[1]<10)):
                ratio = 0.63
            else:
                ratio = 0.61
            M = (ratio*vl)**2*np.eye(2) + (1-ratio**2)*(vl**2)*np.outer(f,f)
            D[p,i,:,:] = np.linalg.inv(M)
        i = i+1

# stimuli
x_stim = []
x_stim.append(np.argwhere((points[:,0]>=7)*(points[:,0]<=8)*(points[:,1]>=10)*(points[:,1]<=11))[:,0])
x_stim.append(np.argwhere((points[:,0]>=7)*(points[:,0]<=8)*(points[:,1]>=10)*(points[:,1]<=11))[:,0])
t_stim = np.array([0,290,int(1e5)])

#initialization
phi = np.inf*np.ones(n_points)
phi_comp = np.inf*np.ones((n_points,max_neighs))
v = np.zeros(n_points)
scar = np.argwhere(p2s==1)[:,0]
v[scar] = -1
di = np.inf*np.ones(n_points)
apd = np.zeros(n_points)
eps = 1e-10*np.ones(1)
cv = eps*np.ones((n_points,max_elems))
h = heap(100*n_points) 
ids = -np.ones(n_points)
    
# loop
dt = 1
t = t_stim[0]
T = 2000
stim = 0
scale = 0.9
start = time.time()
with meshio.xdmf.TimeSeriesWriter("output/eikonal_scar.xdmf") as writer:
    writer.write_points_cells(points, [("triangle",elems)])
    writer.write_data(t, point_data={"v": v})
    while t<T:
        # stimuli
        if t>=t_stim[stim]:
            xs = x_stim[stim]
            xs_inds = xs[np.argwhere(phi[xs]==np.inf)][:,0]
            phi[xs_inds] = t
            for x in xs_inds:
                if p2s[x]==0:
                    apd[x] = np.interp(di[x],DIres_h,APDres_h)
                else:
                    apd[x] = np.interp(di[x],DIres_bz,APDres_bz)
            di[xs_inds] = 0
            v[xs_inds] = 1
            # include the excitable neighbors in the heap
            for x in xs_inds:
                for p in p2p[x]:
                    if phi[p]==np.inf:
                        Pp = np.array(list(p2p[p])) 
                        phi_comp[p,np.argwhere(Pp==x)[0][0]] = phi[x]
                        Tx = p2t[x] 
                        Tp = p2t[p]
                        TxTp = np.intersect1d(Tx,Tp)
                        for i in range(len(TxTp)): 
                            if cv[p,np.argwhere(Tp==TxTp[i])[0][0]]==eps:
                                if p2s[p]==0:
                                    cv[p,np.argwhere(Tp==TxTp[i])[0][0]] = np.interp(di[p],DIres_h,scale*CVres_h)
                                else:
                                    cv[p,np.argwhere(Tp==TxTp[i])[0][0]] = np.interp(di[p],DIres_bz,scale*CVres_bz)
                        if ids[p]==-1:
                            ids[p] = h.push(p,hopf_lax_update(p,phi_comp[p,:],cv[p,:]))
                        else:
                            h.update(ids[p],hopf_lax_update(p,phi_comp[p,:],cv[p,:]))
            stim = stim+1
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
            if p2s[xa]==0:
                apd[xa] = np.interp(di[xa],DIres_h,APDres_h)
            else:
                apd[xa] = np.interp(di[xa],DIres_bz,APDres_bz)
            di[xa] = 0
            v[xa] = 1
            # include the excitable neighbors in the heap
            for p in p2p[xa]:
                if phi[p]==np.inf:
                    Pp = np.array(list(p2p[p]))
                    phi_comp[p,np.argwhere(Pp==xa)[0][0]] = phi[xa]
                    Ta = p2t[xa]
                    Tp = p2t[p]
                    TaTp = np.intersect1d(Ta,Tp)
                    for i in range(len(TaTp)):
                        if cv[p,np.argwhere(Tp==TaTp[i])[0][0]]==eps:
                            if p2s[p]==0:
                                cv[p,np.argwhere(Tp==TaTp[i])[0][0]] = np.interp(di[p],DIres_h,scale*CVres_h)
                            else:
                                cv[p,np.argwhere(Tp==TaTp[i])[0][0]] = np.interp(di[p],DIres_bz,scale*CVres_bz)
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
        inds = ((p2s==0)*(phi<np.inf)*(di>=di_min_h)+(p2s!=0)*(phi<np.inf)*(di>=di_min_bz)).nonzero()[0]
        phi[inds] = np.inf
        writer.write_data(t, point_data={"v": v})
comp_time = time.time()-start
print(comp_time)
