#!/usr/bin/env python
""" MultiQC hook functions - we tie into the MultiQC
core here to add in extra functionality. """

import logging
import yaml
from os.path import join, dirname
from multiqc.utils import report, config

from pkg_resources import get_distribution
__version__ = get_distribution("multiqc_az").version

log = logging.getLogger('multiqc.multiqc_az')


class config_loaded:
    def __init__(self):
        log.debug("Running config_loaded hook v{}. Loading specific settings and metadata".format(__version__))

        with open(join(dirname(__file__), 'multiqc_config.yaml')) as f:
            cfg = yaml.load(f)
        config.update_dict(config.__dict__, cfg)


class execution_start:
    def __init__(self):
        az_conf = config.__dict__.get('az')
        if az_conf:
            report.az = az_conf
            if az_conf.get('is_rnaseq', False) is True:
                config.table_columns_visible['bcbio']['Mapped_reads'] = True
                config.table_columns_visible['QualiMap']['median_coverage'] = True
                config.table_columns_visible['FastQC']['percent_gc'] = True
                config.table_columns_visible['Samtools Stats']['error_rate'] = True


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

        az_conf = config.__dict__.get('az')

        if az_conf and 'gender_by_sample' in az_conf:
            log.info('Adding Gender metrics from NGS_Reports')
            report.general_stats_data.append({
                sname: {'gender': data} for sname, data
                in az_conf['gender_by_sample'].items()
            })
            report.general_stats_headers.append({
                'gender': {
                    'title': 'Gender',
                    'description': 'Gender based on coverage of specific chrY regions (if they overlap with target)',
                }
            })

        if az_conf and 'prealign' in az_conf:
            config.report_header_info = [
                d for d in config.report_header_info
                if list(d.keys())[0] != 'Target for var. calling:']


class after_set_general_stats_html:
    def __init__(self):
        az_conf = config.__dict__.get('az')
        if az_conf and 'ngs_report_by_sample' in az_conf:
            log.info('Adding NGS repots links')
            for sname in az_conf['ngs_report_by_sample']:
                if az_conf['ngs_report_by_sample'][sname] is not None:
                    report.general_stats_html = report.general_stats_html \
                        .replace(
                            '>' + sname + '<',
                            '><a href="' + az_conf['ngs_report_by_sample'][sname] + '">' + sname + '</a><')
            report.ngs_reports_added = True

            if len(az_conf['ngs_report_by_sample'].items()) >= config.max_table_rows:
                new_general_stats_html = ''
                for l in report.general_stats_html.split('\n'):
                    if 'http://multiqc.info/docs/#tables--beeswarm-plots' in l:
                        l = ('<p class="text-muted"><span class="glyphicon glyphicon-exclamation-sign" '
                             'data-toggle="tooltip"></span> A <a href="http://multiqc.info/docs/#tables--beeswarm-plots"> '
                             'beeswarm</a> plot has been generated instead because of the large number of samples. '
                             'Showing {} samples:<br>').format(len(az_conf['ngs_report_by_sample']))
                        sample_lines = []
                        for sname in az_conf['ngs_report_by_sample']:
                            if az_conf['ngs_report_by_sample'][sname] is not None:
                                sample_lines.append('<a href="' + az_conf['ngs_report_by_sample'][sname] + '">' + sname + '</a>')
                            else:
                                sample_lines.append(' <span>' + sname + '</span>')
                        l += ',&nbsp'.join(sample_lines) + '</p>'
                    new_general_stats_html += l
                report.general_stats_html = new_general_stats_html
                report.beeswarm_renderred = True

                # TODO beeswarm plot:
                # 1. hide hidden columns
                # 2. point on sample -> highlight dot in chart
                # 3. point on dot in chart -> show sample name in a tooltip
