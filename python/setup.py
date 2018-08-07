#!/usr/bin/env python
import os
import sys
import json
import glob
from setuptools import setup
if sys.version_info.major >= 3:
    from setuptools import Extension
    try:
        # for pip < 10.0
        from pip import locations
    except ImportError:
        # for pip >= 10.0
        from pip._internal import locations

version = '1.2.0'

basedir = os.path.abspath(os.path.dirname(__file__))

args = {
    'name': 'ngt',
    'version': version,
    'author': 'ngt-pj',
    'author_email': 'https://www.yahoo-help.jp/',
    'url': 'https://github.com/yahoojapan/NGT',
    'license': 'Apache License Version 2.0',
    'packages': ['ngt']
}

if sys.version_info.major >= 3:
    module1 = Extension('ngtpy', 
                        include_dirs=['/usr/local/include', 
                                      os.path.dirname(locations.distutils_scheme('pybind11')['headers']),
                                      os.path.dirname(locations.distutils_scheme('pybind11', True)['headers'])],
                        library_dirs=['/usr/local/lib', '/usr/local/lib64'],
                        libraries=['ngt'],
                        extra_compile_args=['-std=c++11', '-mavx', '-Ofast', '-march=native', '-lrt', '-DNDEBUG'],
                        sources=['src/ngtpy.cpp'])
    args['ext_modules'] = [module1]

setup_arguments = args

if os.path.isdir('scripts'):
    setup_arguments['scripts'] = [
        os.path.join('scripts', f) for f in os.listdir('scripts')
    ]

if __name__ == '__main__':
    setup(**setup_arguments)
