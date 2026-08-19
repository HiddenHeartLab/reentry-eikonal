import matplotlib.pyplot as plt
import numpy as np
from ioniclib import IonicModel
from tqdm.auto import trange

nnodes = 1

vm  = np.zeros((nnodes,),dtype=np.double)
Isd = np.zeros_like(vm)
Imi = np.zeros_like(vm)
dtime = np.zeros((nnodes,),dtype=np.single)

dt = 0.01
ion_model = IonicModel(nnodes,dt=dt,model='cm98')
ion_model.set_parameter('g_CaL',0.037)
ion_model.set_parameter('g_K1',0.18)
ion_model.set_parameter('g_to',0.026)
ion_model.set_parameter('g_Ks',0.1)
ion_model.set_parameter('g_Na',8.0)

ion_model.set_initial_conditions(vm)

bcl = 400.0

sol = []
ndt = 180_000
for k in trange(ndt):
    t = k * dt
    Isd[:] = 30.0 * (1.0 <= np.mod(t, bcl) <= 2.0)
    ion_model.step(t, vm, ion_model._yyy, Isd, dtime, Imi)
    vm += dt*(Isd - Imi)
    sol.append(vm[0])

sol = np.array(sol)
tt = dt*np.arange(ndt)
plt.plot(tt, sol)
plt.show()
print(sol)