from distutils.core import setup, Extension

setup(
    ext_modules=[Extension("VectorMath", ["VectorMath.c", "backend\\vector.c"])],
    include_dirs=None,
    )
