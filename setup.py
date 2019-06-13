#!/usr/bin/env python

from distutils.core import setup

try:
    # for pip >= 10
    from pip._internal.req import parse_requirements
except ImportError:
    # for pip <= 9.0.3
    from pip.req import parse_requirements


def load_requirements(fname):
    reqs = parse_requirements(fname, session="test")
    return [str(ir.req) for ir in reqs]


setup(name='dna-pipeline',
      version='1.0dev',
      description='Workflow for single-cell CNV data analysis',
      author='Mustafa Anil Tuncel',
      author_email='tuncel.manil@gmail.com',
      # url='https://www.python.org/sigs/distutils-sig/',
      packages=['pipeline', 'pipeline.secondary_analysis'],
      install_requires=load_requirements("requirements.txt")

     )