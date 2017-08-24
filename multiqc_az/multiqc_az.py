#!/usr/bin/env python
""" MultiQC hook functions - we tie into the MultiQC
core here to add in extra functionality. """

from collections import OrderedDict
import logging
import json
import os
import re
import requests
import shutil
import socket
import subprocess
import sys
import yaml
from os.path import join, dirname

from pkg_resources import get_distribution
__version__ = get_distribution("multiqc_az").version

from multiqc.utils import report, util_functions, config

from .utils import format_bed_info

log = logging.getLogger('multiqc.multiqc_az')


class config_loaded:
    def __init__(self):
        log.debug("Running config_loaded hook v{}. Loading specific settings and metadata".format(__version__))

        with open(join(dirname(__file__), 'multiqc_config.yaml')) as f:
            cfg = yaml.load(f)
        config.update_dict(config.__dict__, cfg)


class load_metadata:
    def __init__(self):
        if config.kwargs.get('az_metadata'):
            with open(config.kwargs['az_metadata']) as f:
                metadata = yaml.load(f)
            report.az = metadata
            config.title = report.az['project_name']

            if 'coverage_thresholds' in report.az and 'coverage_thresholds_hidden' in report.az:
                if 'qualimap_config' not in config.__dict__:
                    config.qualimap_config = dict()
                config.qualimap_config['general_stats_coverage'] = report.az['coverage_thresholds']
                config.qualimap_config['general_stats_coverage_hidden'] = report.az['coverage_thresholds_hidden']
            
            if 'preseq' in report.az:
                config.preseq = config.update_dict(config.__dict__.get('preseq', {}), report.az['preseq'])
                
            if report.az.get('is_rnaseq', False) is True:
                config.table_columns_visible['bcbio']['Mapped_reads'] = True
                config.table_columns_visible['QualiMap']['median_coverage'] = True
                config.table_columns_visible['FastQC']['percent_gc'] = True
                config.table_columns_visible['Samtools Stats']['error_rate'] = True

        else:
            log.warn('Warning: --az-metadata option is not provided to MultiQC')


# TODO: move descriptions to proper places according to changes in MultiQC v1.2
general_stats_descriptions = {
    '5_3_bias':
        '5\'/3\' bias. Significant deviation of 5\'/3\' ratios from the ' +
        'value of 1 could indicate biases which may have been introduced ' +
        'during the library preparation, RNA degradation, etc.'
}


class before_set_general_stats_html:
    def __init__(self):
        for header in report.general_stats_headers:
            for k, v in general_stats_descriptions.items():
                if k in header:
                    header[k]['description'] = v

        if 'az' in report.__dict__:
            if 'gender_by_sample' in report.az:
                log.info('Adding Gender metrics from NGS_Reports')
                report.general_stats_data.append({
                    sname: {'gender' : data} for sname, data
                    in report.az['gender_by_sample'].items()
                })
                report.general_stats_headers.append({
                    'gender': {
                        'title': 'Gender',
                        'description': 'Gender based on coverage of specific chrY regions (if they overlap with target)',
                    }
                })


class after_set_general_stats_html:
    def __init__(self):
        if 'az' in report.__dict__ and 'ngs_report_by_sample' in report.az:
            log.info('Adding NGS repots links')
            for sname in report.az['ngs_report_by_sample']:
                if report.az['ngs_report_by_sample'][sname] is not None:
                    report.general_stats_html = report.general_stats_html \
                        .replace(
                            '>' + sname + '<',
                            '><a href="' + report.az['ngs_report_by_sample'][sname] + '">' + sname + '</a><')
            report.ngs_reports_added = True

            if len(report.az['ngs_report_by_sample'].items()) >= config.max_table_rows:
                new_general_stats_html = ''
                for l in report.general_stats_html.split('\n'):
                    if 'http://multiqc.info/docs/#tables--beeswarm-plots' in l:
                        l = ('<p class="text-muted"><span class="glyphicon glyphicon-exclamation-sign" '
                             'data-toggle="tooltip"></span> A <a href="http://multiqc.info/docs/#tables--beeswarm-plots"> '
                             'beeswarm</a> plot has been generated instead because of the large number of samples. '
                             'Showing {} samples:<br>').format(len(report.az['ngs_report_by_sample']))
                        sample_lines = []
                        for sname in report.az['ngs_report_by_sample']:
                            if report.az['ngs_report_by_sample'][sname] is not None:
                                sample_lines.append('<a href="' + report.az['ngs_report_by_sample'][sname] + '">' + sname + '</a>')
                            else:
                                sample_lines.append(' <span>' + sname + '</span>')
                        l += ',&nbsp'.join(sample_lines) + '</p>'
                    new_general_stats_html += l
                report.general_stats_html = new_general_stats_html

                # TODO beeswarm plot:
                # 1. hide hidden columns
                # 2. point on sample -> highlight dot in chart
                # 3. point on dot in chart -> show sample name in a tooltip

            # fixed_stats_data = []
            # for section in report.general_stats_data:
            #     fixed_section = dict()
            #     for sname, data in section.items():  # section = {sample: {metric: name, metric: name, ..}, sample: {}, ...}
            #         # fixing sample names to link to Oncology NGS reports
            #         if sname in report.az['ngs_report_by_sample']:
            #             sname = '<a href="' + report.az['ngs_report_by_sample'][sname] + '">' + sname + '</a>'
            #         fixed_section[sname] = data
            #     fixed_stats_data.append(fixed_section)
            # report.general_stats_data = fixed_stats_data
