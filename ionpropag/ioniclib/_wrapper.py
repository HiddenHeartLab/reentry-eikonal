from ._findlib import findlib
from ctypes import Structure, POINTER, CDLL, byref
from ctypes import c_int, c_char, c_void_p, c_ubyte, c_float, c_double
from numpy.ctypeslib import ndpointer
import json
import os
import numpy as np

class Membrane_info(Structure):
    _fields_ = [("Nsvar",c_int),
                ("Nvar",c_int),
                ("Nstypes",c_int),
                ("Tname",POINTER(POINTER(c_char))),
                ("Nparam",c_int),
                ("param",POINTER(c_int)),
                ("info",c_void_p),
                ("init",c_void_p),
                ("infinite",c_void_p),
                ("step",c_void_p),
                ("get_status",c_void_p),
                ("ccode",c_void_p),
                ("can_use_cuda",c_int),
                ("init_cuda",c_void_p)]

class Membrane_cell_info(Structure):
    _fields_ = [("mcode",c_ubyte),
                ("ccode",c_ubyte),
                ("param",POINTER(c_float))]


class IonicModel(object):
    
    def __init__(self, nnodes=10, dt=0.01, model='cm98'):

        # load the library
        pylib = findlib("propag_wrapper")
        lib = CDLL(pylib)

        # generate parameters
        param_json = os.path.join(os.path.dirname(__file__), f'{model}.json')
        with open(param_json,"r") as fi:
            prm = json.load(fi)

        ctype_map = { 'int': c_int, 'float': c_float, 'double': c_double, 'boolean': c_int }
        fields = [ (m[0], ctype_map[m[1]] if m[2]==0 else ctype_map[m[1]] * m[2])
                   for m in prm['members'] ]
        param_cls = type(prm['name'],(Structure,),{'_fields_': fields})
        params = param_cls()

        # set default values of parameters
        for m in prm['members']:
            val = m[3]
            if isinstance(val,list):
                val = (ctype_map[m[1]] * m[2])(*(val))
            elif not isinstance(val,list) and m[2]>0:
                val = (ctype_map[m[1]] * m[2])(*([float(val)]*m[2]))
            #if m[2] > 0:
            #    val = (ctype_map[m[1]] * m[2])(*([float(val)]*m[2]))
            setattr(params,m[0],val)

        self.params = params

        # types
        yyy_t   = ndpointer(dtype=np.double,ndim=2,flags='C')
        vm_t    = ndpointer(dtype=np.double,ndim=1,flags='C')
        dtime_t = ndpointer(dtype=np.single,ndim=1,flags='C')

        # library
        getattr(lib,f'{model}_info').argtypes = [ POINTER(Membrane_cell_info), 
                                   POINTER(Membrane_info),
                                   POINTER(param_cls) ]
        getattr(lib,f'{model}_info').restype  = None

        getattr(lib,f'{model}_init').argtypes = [ c_float ]
        getattr(lib,f'{model}_init').restype  = None

        getattr(lib,f'{model}_get_status').argtypes = [ c_float, # vm
                                                        yyy_t,   # cell_status
                                                        POINTER(Membrane_cell_info),
                                                        POINTER(c_float), # Stats
                                                        POINTER(c_int),   # Nstats
                                                        POINTER(c_char) ] # Names
        getattr(lib,f'{model}_get_status').restype  = None

        lib.ion_infinite.argtypes = [ POINTER(Membrane_info), 
                                      POINTER(Membrane_cell_info),
                                      vm_t,
                                      yyy_t,
                                      c_int ]
        lib.ion_infinite.restype  = None

        lib.ion_step.argtypes = [ POINTER(Membrane_info), 
                                  POINTER(Membrane_cell_info),
                                  vm_t,     # Vm
                                  yyy_t,    # yyy
                                  vm_t,     # Isd
                                  dtime_t,  # dtime
                                  vm_t,     # Imi (output)
                                  c_double, # dt
                                  c_double, # simtime
                                  c_int ]   # nnodes
        lib.ion_step.restype  = None

        # initialize the model
        ct  = Membrane_cell_info()
        ifo = Membrane_info()
        getattr(lib,f'{model}_info')(byref(ct),byref(ifo),byref(params))
        getattr(lib,f'{model}_init')(dt)

        # initialize internal variables
        self._yyy = np.zeros((nnodes,ifo.Nsvar),dtype=np.double)
        self._dt  = dt
        self._ct  = ct
        self._ifo = ifo
        self._lib = lib


    def set_initial_conditions(self, vm):

        self._lib.ion_infinite(byref(self._ifo),
                               byref(self._ct),
                               vm,
                               self._yyy,
                               self._yyy.shape[0])
           
    
    #def step(self, t, vm, Isd, dtime, Imi):

    #    self._lib.ion_step(byref(self._ifo),
    #                       byref(self._ct),
    #                       vm,
    #                       self._yyy,
    #                       Isd,
    #                       dtime,
    #                       Imi,
    #                       self._dt,
    #                       t,
    #                       self._yyy.shape[0])
        
    #    return self._yyy

    def step(self, t, vm, y, Isd, dtime, Imi):

        self._lib.ion_step(byref(self._ifo),
                           byref(self._ct),
                           vm,
                           y,
                           Isd,
                           dtime,
                           Imi,
                           self._dt,
                           t,
                           y.shape[0])


    def get_parameter(self, pname):
        val = getattr(self.params,pname)
        if hasattr(val, '__len__'):
            return list(val)
        else:
            return val
            
    def set_parameter(self, pname, pvalue):
        val = getattr(self.params, pname)
        if hasattr(val, '__len__'):
            if hasattr(pvalue, '__len__'):
                setattr(self.params, pname, type(val)(*pvalue))
            else:
                setattr(self.params, pname, type(val)(*[pvalue]*len(val)))
        else:
            setattr(self.params, pname, pvalue)

