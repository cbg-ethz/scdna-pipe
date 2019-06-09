#!/usr/bin/env python

from distutils.core import setup

setup(name='dna-pipeline',
      version='1.0dev',
      description='Workflow for single-cell CNV data analysis',
      author='Mustafa Anil Tuncel',
      author_email='tuncel.manil@gmail.com',
      # url='https://www.python.org/sigs/distutils-sig/',
      packages=['pipeline', 'pipeline.secondary_analysis'],
     )