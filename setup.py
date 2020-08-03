#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open("README.rst") as readme_file:
    readme = readme_file.read()

with open("HISTORY.rst") as history_file:
    history = history_file.read()

requirements = [
    "Click>=7.0",
    "pysam>=0.15.4",
    "dataclasses",
]

setup(
    version="0.5.7",
    author="Warren W. Kretzschmar",
    author_email="winni@warrenwk.com",
    python_requires=">=3.6",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    description="Create multi-sample text-pileups of streaming SAM/BAM files.",
    entry_points={"console_scripts": ["spileup=streaming_pileupy:main",],},
    install_requires=requirements,
    license="MIT license",
    long_description=readme + "\n\n" + history,
    include_package_data=True,
    keywords="streaming_pileupy",
    name="streaming_pileupy",
    py_modules=["streaming_pileupy"],
    package_dir={"": "src"},
    url="https://github.com/winni2k/streaming_pileupy",
    zip_safe=False,
)
