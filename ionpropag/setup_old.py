import setuptools
from setuptools import Extension
import os
from ioniclib._prm import parse_prs
import json
import re

PROPAG_DIR  = os.path.expanduser("src")

# models
models = ['tp06','cm98','tnnp','ms']
models_src = [ os.path.join(PROPAG_DIR, f"{m}.c")
               for m in models ]

print(models_src)

# generate PRM files
for model in models:
    print(f"GENERATING {model} parameters...")
    prs_ifile = os.path.join(PROPAG_DIR, f"{model}.prs")
    prs_ofile = f"ioniclib/{model}.json"

    with open(prs_ifile,"r") as fi:
        prs = fi.read()
    prs = parse_prs(prs)

    params = {'name': prs.name, 'members': []}
    for m in prs.members:
        val = m.default
        typ = m.type
        mtc = re.match(r'^\$(float|double|int|boolean)(?:|\[([0-9]+)\])$',typ)
        if mtc:
            typ = mtc.group(1)
            if typ == 'int':
                val = [int(v) for v in val] if isinstance(val,list) else int(val)
            elif typ == 'float':
                val = [float(v) for v in val] if isinstance(val,list) else float(val)
            elif typ == 'double':
                val = [float(v) for v in val] if isinstance(val,list) else float(val)
            elif typ == 'boolean':
                val = 1 if val == "TRUE" else 0
            arr = int(mtc.group(2)) if mtc.group(2) is not None else 0
        else:
            raise RuntimeError(f"Cannot parse {m.type} in PRM file")

        params['members'].append( (m.name, typ, arr, val) )

    # NOTE: prm reverses the order of the parameters in the .h file
    params['members'] = params['members'][::-1]
    
    with open(prs_ofile,"w") as fo:
        json.dump(params, fo)

# ioniclib = Extension("propag_wrapper",
#                      models_src + ["src/ionlib.c"],
#                      include_dirs=[PROPAG_DIR,PRM_INCLUDE],
#                      define_macros=[('finite','isfinite')])

ioniclib = Extension(
    "propag_wrapper",
    sources=models_src + ["src/ionlib.c"],
    include_dirs=[PROPAG_DIR, PRM_INCLUDE, "/opt/homebrew/opt/libomp/include"],
    define_macros=[("finite", "isfinite")],
    extra_compile_args=["-Xpreprocessor", "-fopenmp"],
    extra_link_args=["-lomp", "-L/opt/homebrew/opt/libomp/lib"],
)

setuptools.setup(
    name="IonPropag",
    version="0.5",
    ext_modules=[ioniclib])
