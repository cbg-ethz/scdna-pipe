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
    dependency_links=['git+git://github.com/anilbey/PhenoGraph.git@7ef72746688217b2feed2f3f8226da2f304c532c#egg=Phenograph']
)
