from setuptools import setup
from Cython.Build import cythonize

setup(
    name='Genome Rearrangements',
    ext_modules=cythonize(["structures/*.py", "algorithms/*.py"]),
    zip_safe=False,
)