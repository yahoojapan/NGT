#!/usr/bin/env python
import os
import sys
import json
import shutil
import glob
import setuptools
import pybind11

static_library_option = '--static-library'
included_library_option = '--included-library'
version_file = 'VERSION'

static_library = False
if static_library_option in sys.argv:
    print('use the NGT static library')
    sys.argv.remove(static_library_option)
    static_library = True

included_library = False
if included_library_option in sys.argv:
    print('use the NGT included library')
    sys.argv.remove(included_library_option)
    included_library = True

if sys.version_info.major >= 3:
    from setuptools import Extension

if os.path.isfile('../' + version_file):
    shutil.copyfile('../' + version_file, version_file)

with open(version_file, 'r') as fh:
    version = fh.read().rstrip('\n')

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
    'install_requires': ['numpy', 'pybind11']
}

if sys.version_info.major >= 3:
    if static_library or included_library:
        params = {
            'include_dirs': ['/usr/local/include',
                             pybind11.get_include(True),
                             pybind11.get_include(False)],
            'extra_compile_args': ['-std=c++11', '-Ofast', '-fopenmp', '-lrt', '-DNDEBUG'],
            'sources': ['src/ngtpy.cpp']
        }
    else:
        params = {
            'include_dirs': ['/usr/local/include',
                             pybind11.get_include(True),
                             pybind11.get_include(False)],
            'extra_compile_args': ['-std=c++11', '-Ofast', '-fopenmp', '-march=native', '-lrt', '-DNDEBUG'],
            'sources': ['src/ngtpy.cpp']
        }

    dynamic_lib_params = {
        'library_dirs': ['/usr/local/lib', '/usr/local/lib64'],
        'libraries': ['ngt', 'gomp']
    }
    included_lib_params = {
        'library_dirs': ['/usr/local/lib', '/usr/local/lib64'],
        'libraries': ['ngt', 'gomp'],
        'extra_link_args': ['-static-libstdc++']
    }
    static_lib_params = {
        'extra_objects': ['../build/lib/NGT/libngt.a'],
        'libraries': ['gomp']
    }
    if static_library:
        params.update(static_lib_params)
    elif included_library:
        params.update(included_lib_params)
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
