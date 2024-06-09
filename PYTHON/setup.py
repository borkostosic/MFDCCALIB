import torch
import setuptools
from setuptools import setup
from torch.utils.cpp_extension import BuildExtension, CUDAExtension

setup(
    name='mfdcca',
    ext_modules=[
        CUDAExtension(
            name='mfdcca',
            sources=['mfdcca.cpp'],
            extra_compile_args={'cxx':['-O3'],})
    ],
    cmdclass={
        'build_ext': BuildExtension
})
