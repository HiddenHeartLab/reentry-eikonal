from pathlib import Path
import sys
import platform
from setuptools import Extension, setup

ROOT = Path(__file__).resolve().parent
# Local build helpers. These are part of the source tree, not
# dependencies of the installed package.
sys.path.insert(0, str(ROOT))

from build_tools.prm import generate_parameters

SRC = ROOT / "src"
PACKAGE = ROOT / "ioniclib"

MODELS = ["tp06", "cm98", "tnnp", "ms"]

for model in MODELS:
    generate_parameters(SRC / f"{model}.prs", PACKAGE / f"{model}.json")

def extension():
    sources = [
        str(SRC / f"{model}.c")
        for model in MODELS
    ]
    sources.append(str(SRC / "ionlib.c"))

    include_dirs = [
        str(SRC),
    ]

    extra_compile_args = []
    extra_link_args = []

    # Apple Clang does not ship OpenMP itself.
    if platform.system() == "Darwin":
        libomp = Path("/opt/homebrew/opt/libomp")

        if libomp.exists():
            include_dirs.append(str(libomp / "include"))
            extra_compile_args += [
                "-Xpreprocessor",
                "-fopenmp",
            ]
            extra_link_args += [
                "-lomp",
                f"-L{libomp / 'lib'}",
            ]
    else:
        extra_compile_args.append("-fopenmp")
        extra_link_args.append("-fopenmp")

    return Extension(
        "propag_wrapper",
        sources = [*(f"src/{model}.c" for model in MODELS),
                   "src/ionlib.c"],
        include_dirs=include_dirs,
        define_macros=[
            ("finite", "isfinite"),
        ],
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
    )

setup(
    ext_modules=[extension()],
)