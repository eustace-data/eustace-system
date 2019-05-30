"""Setup script required to build C extension for fast gridding.
   To compile: python2.7 setup.py build_ext --inplace"""

from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

include_dirs = [numpy.get_include()]

extensions_aggregate = [
    Extension(
        'eustace.satgrid.aggregate', 
        sources=['eustace/satgrid/aggregate.c'], 
        include_dirs=include_dirs
    )
]

setup(
    name='eustace.satgrid.aggregate', 
    ext_modules=extensions_aggregate
)

extensions_gridbin = [ 
    Extension(
        'eustace.satgrid_iris.gridbin.operation.sumcountmaxmin', 
        sources=['eustace/satgrid_iris/gridbin/operation/sumcountmaxmin.pyx'], 
        include_dirs=include_dirs
    ),
    Extension(
        'eustace.satgrid_iris.gridbin.operation.count', 
        sources=['eustace/satgrid_iris/gridbin/operation/count.pyx'], 
        include_dirs=include_dirs
    ),
    Extension(
        'eustace.satgrid_iris.gridbin.operation.sumsq', 
        sources=['eustace/satgrid_iris/gridbin/operation/sumsq.pyx'], 
        include_dirs=include_dirs
    ),
    Extension(
        'eustace.satgrid_iris.gridbin.operation.sumsqdev', 
        sources=['eustace/satgrid_iris/gridbin/operation/sumsqdev.pyx'], 
        include_dirs=include_dirs
    )
]

setup(
    name='eustace.satgrid_iris.gridbin', 
    ext_modules=cythonize(extensions_gridbin, compiler_directives={'embedsignature': True})
)

extensions_analysis = [
    Extension(
        'eustace.analysis.diagnosticgrid_aggregate', 
        sources=['eustace/analysis/diagnosticgrid_aggregate.pyx'], 
        include_dirs=include_dirs
    ),
    Extension(
        'eustace.analysis.mesh.locator_search_planes', 
        sources=['eustace/analysis/mesh/locator_search_planes.pyx'], 
        include_dirs=include_dirs
    )
]

setup(
    name='eustace.analysis', 
    ext_modules=cythonize(extensions_analysis, compiler_directives={'embedsignature': True})
)

extensions_analysis_mesh = [
    Extension(
        'eustace.analysis.mesh.locator_search_planes', 
        sources=['eustace/analysis/mesh/locator_search_planes.pyx'], 
        include_dirs=include_dirs
    )
]

setup(
    name='eustace.analysis.mesh', 
    ext_modules=cythonize(extensions_analysis_mesh, compiler_directives={'embedsignature': True})
)

extensions_analysis_linalg = [
    Extension(
        'eustace.analysis.advanced_standard.linalg.qinv', 
        sources=['eustace/analysis/advanced_standard/linalg/qinv.pyx'], 
        include_dirs=include_dirs
    )
]

setup(
    name='eustace.analysis.advanced_standard.linalg', 
    ext_modules=cythonize(extensions_analysis_linalg, compiler_directives={'embedsignature': True})
)
