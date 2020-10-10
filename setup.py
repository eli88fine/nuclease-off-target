# -*- coding: utf-8 -*-
"""Setup configuration."""
from setuptools import find_packages
from setuptools import setup


setup(
    name="nuclease-off-target",
    version="0.2.1",
    description="Analyzing Nuclease Off-Target Sites and Activity",
    url="https://github.com/eli88fine/nuclease-off-target",
    author="Eli Fine",
    author_email="ejfine@gmail.com",
    license="MIT",
    packages=find_packages("src"),
    package_dir={"": "src"},
    install_requires=[
        "requests>=2.24.0",
        "biopython>=1.78",
        "beautifulsoup4>=4.9.3",
        "parasail>=1.2",
        "stdlib_utils>=0.3.1",
    ],
    zip_safe=False,
    include_package_data=True,
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering",
    ],
)
