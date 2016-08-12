#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'numpy',# TODO: put package requirements here
]

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='ifu_met_parser',
    version='0.1.0',
    description="A compilation of parsers for meteorological observation data",
    long_description=readme + '\n\n' + history,
    author="Christian Chwala",
    author_email='christian.chwala@kit.edu',
    url='https://github.com/cchwala/ifu_met_parser',
    packages=[
        'ifu_met_parser',
    ],
    package_dir={'ifu_met_parser':
                 'ifu_met_parser'},
    include_package_data=True,
    install_requires=requirements,
    license="MIT license",
    zip_safe=False,
    keywords='ifu_met_parser',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    test_suite='tests',
    tests_require=test_requirements
)
