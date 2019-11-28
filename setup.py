from setuptools import setup, find_packages

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name="scgenpy",
    version='1.0dev',
    description='Package for single-cell CNV data analysis',
    author=['Mustafa Anil Tuncel', 'Pedro Falé Ferreira'],
    author_email=['tuncel.manil@gmail.com', 'pedro.ferreira@bsse.ethz.ch'],
    packages=find_packages(),
    install_requires=requirements
)
