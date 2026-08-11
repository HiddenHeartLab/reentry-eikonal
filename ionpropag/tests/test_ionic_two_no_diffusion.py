from ioniclib import IonicModel
import numpy as np

nnodes1 = 2
nnodes2 = 2

vm1  = np.zeros((nnodes1,),dtype=np.double)
Isd1 = np.zeros_like(vm1)
Imi1 = np.zeros_like(vm1)
dtime1 = np.zeros((nnodes1,),dtype=np.single)

vm2  = np.zeros((nnodes2,),dtype=np.double)
Isd2 = np.zeros_like(vm2)
Imi2 = np.zeros_like(vm2)
dtime2 = np.zeros((nnodes2,),dtype=np.single)

dt = 0.01
ion_model1 = IonicModel(nnodes1,dt=dt,model='cm98')
ion_model1.set_initial_conditions(vm1)

ion_model2 = IonicModel(nnodes2,dt=dt,model='tnnp')
ion_model2.set_initial_conditions(vm2)

print(np.r_[vm1,vm2])

for _ in range(10):
    ion_model1.step(0.0, vm1, ion_model1._yyy, Isd1, dtime1, Imi1)
    ion_model2.step(0.0, vm2, ion_model1._yyy, Isd2, dtime2, Imi2)
    vm1 += dt*Imi1
    vm2 += dt*Imi2
    print(np.r_[vm1,vm2])
