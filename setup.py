from setuptools import setup, find_packages

with open("requirements.txt") as f:
    requirements = f.read().splitlines()

setup(
    name="dnapipe",
    version="1.0dev",
    description="Workflow for single-cell CNV data analysis",
    author="Mustafa Anil Tuncel",
    author_email="tuncel.manil@gmail.com",
    packages=find_packages(),
    install_requires=requirements,
)
