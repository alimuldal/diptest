import numpy as np
from distutils.core import setup
from Cython.Distutils import Extension, build_ext


diptest = Extension(
    name="diptest._diptest",
    sources=["diptest/_dip.c", "diptest/_diptest.pyx"],
    extra_compile_args=['-O3', '-std=c99'],
)

setup(
    name='diptest',
    cmdclass={'build_ext': build_ext},
    ext_modules=[diptest],
    packages=['diptest'],
    package_data={'diptest': ['dip_crit.txt']},
    include_dirs = [np.get_include()]
)