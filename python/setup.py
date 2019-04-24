#!/usr/bin/env python
import os
import sys
import json
import glob
import setuptools

static_library_option = '--static-library'
static_library = False
if static_library_option in sys.argv:
    print('use the NGT static library')
    sys.argv.remove(static_library_option)
    static_library = True

if sys.version_info.major >= 3:
    from setuptools import Extension
    try:
        # for pip < 10.0
        from pip import locations
    except ImportError:
        # for pip >= 10.0
        from pip._internal import locations

version = '1.3.1'

if static_library:
    with open('../VERSION', 'r') as fh:
        version = fh.read()

basedir = os.path.abspath(os.path.dirname(__file__))

with open('README.md', 'r', encoding='utf-8') as fh:
    long_description = fh.read()

args = {
    'name': 'ngt',
    'version': version,
    'author': 'Yahoo! JAPAN research',
    'author_email': 'miwasaki@yahoo-corp.jp',
    'url': 'https://github.com/yahoojapan/NGT',
    'description': 'python NGT',
    'long_description': long_description,
    'long_description_content_type': 'text/markdown',
    'license': 'Apache License Version 2.0',
    'packages': ['ngt'],
    'install_requires': ['numpy']
}

if sys.version_info.major >= 3:
    params = {
        'include_dirs': ['/usr/local/include', 
                        os.path.dirname(locations.distutils_scheme('pybind11')['headers']),
                        os.path.dirname(locations.distutils_scheme('pybind11', True)['headers'])],
        'extra_compile_args': ['-std=c++11', '-Ofast', '-march=native', '-lrt', '-DNDEBUG'],
        'sources': ['src/ngtpy.cpp']
    }
    dynamic_lib_params = {
        'library_dirs': ['/usr/local/lib', '/usr/local/lib64'],
        'libraries': ['ngt']
    }
    static_lib_params = {
        'extra_objects': ['../build/lib/NGT/libngt.a'],
        'libraries': ['gomp']
    }
    if static_library:
        params.update(static_lib_params)
    else:
        params.update(dynamic_lib_params)

    module1 = Extension('ngtpy', **params)
    args['ext_modules'] = [module1]

setup_arguments = args

if os.path.isdir('scripts'):
    setup_arguments['scripts'] = [
        os.path.join('scripts', f) for f in os.listdir('scripts')
    ]

if __name__ == '__main__':
    setuptools.setup(**setup_arguments)
