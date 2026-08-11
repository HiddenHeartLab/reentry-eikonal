import os
from distutils import sysconfig

def findlib(libname):
    "Find dynamic library in common paths."

    suffix = sysconfig.get_config_var('EXT_SUFFIX')
    libpaths = [os.path.dirname(os.path.realpath(__file__)),
                os.path.join(os.path.dirname(__file__), '..'),
                os.path.join(os.path.dirname(__file__), '..', '..')]

    pylib = None
    for p in libpaths:
        tmp = os.path.join(p,libname + suffix)
        if not os.path.isfile(tmp): continue
        pylib = tmp

    return pylib
