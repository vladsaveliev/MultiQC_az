#!/usr/bin/env python
"""
Plugin for MultiQC, providing additional tools which are
specific to AstraZeneca Oncology.
For more information about MultiQC, see http://multiqc.info
"""
import sys
from setuptools import setup, find_packages

version = open('VERSION.txt').read().strip().split('\n')[0]


setup(
    name = 'multiqc_az',
    version = version,
    author = 'Vlad Saveliev',
    author_email = 'vladislav.sav@gmail.com',
    description = "MultiQC plugin for the AstraZeneca",
    long_description = __doc__,
    keywords = 'bioinformatics',
    url = 'https://github.com/AstraZeneca-NGS/NGS_Reporting',
    license = 'MIT',
    packages = find_packages(),
    include_package_data = True,
    install_requires = [
        'simplejson',
        'pyyaml',
    ],
    entry_points = {
        'multiqc.modules.v1': [
            'targqc      = multiqc_az.modules.targqc:MultiqcModule',
            'rnaseq_az   = multiqc_az.modules.rnaseq_az:MultiqcModule',
        ],
        'multiqc.templates.v1': [
            'az = multiqc_az.templates.az',
        ],
        'multiqc.hooks.v1': [
            'config_loaded            = multiqc_az.multiqc_az:config_loaded',
            'execution_start          = multiqc_az.multiqc_az:execution_start',
            'after_modules            = multiqc_az.multiqc_az:before_set_general_stats_html',
            'before_report_generation = multiqc_az.multiqc_az:after_set_general_stats_html',
        ]
    },
    classifiers = [
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Environment :: Web Environment',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: JavaScript',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
    ],
)
