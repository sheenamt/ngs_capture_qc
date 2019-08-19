import glob
import os
import subprocess

from distutils.core import Command
from setuptools import setup, find_packages

subprocess.call(
    ('mkdir -p ngs_capture_qc/data && '
     'git describe --tags --dirty > ngs_capture_qc/data/ver.tmp '
     '&& mv ngs_capture_qc/data/ver.tmp ngs_capture_qc/data/ver '
     '|| rm -f ngs_capture_qc/data/ver.tmp'),
    shell=True, stderr=open(os.devnull, "w"))

from ngs_capture_qc import __version__


class CheckVersion(Command):
    description = 'Confirm that the stored package version is correct'
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        with open('ngs_capture_qc/data/ver') as f:
            stored_version = f.read().strip()

        git_version = subprocess.check_output(
            ['git', 'describe', '--tags', '--dirty']).strip()

        assert stored_version == git_version
        print('the current version is', stored_version)


package_data = ['data/*']

params = {'author': 'Sheena Todhunter',
          'author_email': 'sheena.todhunter@gmail.com',
          'description': 'QC scripts for NGS probe files',
          'name': 'ngs_capture_qc',
          'packages': find_packages(),
          'package_dir': {'ngs_capture_qc': 'ngs_capture_qc'},
          'entry_points': {
              'console_scripts': ['capqc = ngs_capture_qc.scripts.main:main']
          },
          'version': __version__,
          'package_data': {'ngs_capture_qc': package_data},
          'test_suite': 'tests',
          'cmdclass': {'check_version': CheckVersion}
          }

setup(**params)
