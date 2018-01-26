#!/usr/bin/env python
import os
import json
import glob
from setuptools import setup

basedir = os.path.abspath(os.path.dirname(__file__))

version = '1.0.0'

# Create a dictionary of our arguments, this way this script can be imported
#  without running setup() to allow external scripts to see the setup settings.
args = {
    'name': 'ngt',
    'version': version,
    'author': 'ngt-pj',
    'author_email': 'https://www.yahoo-help.jp/',
    'url': 'https://github.com/yahoojapan/NGT',
    'license': 'Apache License Version 2.0',
    'packages': ['ngt'],
    #'namespace_packages': ['ngt'],
}
setup_arguments = args

# Add any scripts we want to package
if os.path.isdir('scripts'):
    setup_arguments['scripts'] = [
        os.path.join('scripts', f) for f in os.listdir('scripts')
    ]


if __name__ == '__main__':
    # We're being run from the command line so call setup with our arguments
    setup(**setup_arguments)
