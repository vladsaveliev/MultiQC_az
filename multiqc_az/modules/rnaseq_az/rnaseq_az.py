#!/usr/bin/env python

""" MultiQC module to add link to Bcl2fastq reports """

import logging
from os.path import join

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc import config
from multiqc.plots import scatter

from ngs_utils.file_utils import remove_quotes


# Initialise the logger
log = logging.getLogger(__name__.replace('multiqc_az', 'multiqc'))


standard_colors = [
    '#0000FF',
    '#008000',
    '#FFA500',
    '#FF00FF',
    '#CCCC00',
    '#800000',
    '#00CCCC',
    '#808080',
    '#800080',
    '#808000',
    '#000080',
    '#008080',
    '#00FF00',
]

def parse_pca_data(pca_fpath):
    pca_data = dict()
    color_by_sample = dict()
    color_by_cond = dict()
    conditions = []
    variances = []
    with open(pca_fpath) as f:
        for i, l in enumerate(f):
            if i == 0:
                continue
            if l.startswith('#'):
                variances = l.split(':')[1].split(',')
                continue
            name, pc1, pc2, group, condition, name = [remove_quotes(field) for field in l.split()]
            pc1 = float(pc1)
            pc2 = float(pc2)
            pca_data[name] = [{
                'x': pc1,
                'y': pc2,
                'name': name
            }]
            if condition not in conditions:
                conditions.append(condition)
            color_by_cond[condition] = color_by_sample[name] = standard_colors[conditions.index(condition) % len(standard_colors)]

        if len(pca_data) == 0:
            log.debug("Couldn't parse contents of PCA data file {}".format(f['fn']))
            return None
    return pca_data, color_by_sample, color_by_cond, variances


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Principal Components Analysis', anchor='rnaseq_az')

        rnaseq_pca_files = self.find_log_files('rnaseq_az/pca_data', filecontents=False)
        rnaseq_pca_files = [f for f in rnaseq_pca_files if f]

        if not rnaseq_pca_files:
            log.debug("Could not find the PCA data file in {}".format(config.analysis_dir))
            raise UserWarning
        if len(rnaseq_pca_files) > 1:
            log.warning("More than 1 PCA data file found in {}".format(config.analysis_dir))
            raise UserWarning
        rnaseq_pca_file = rnaseq_pca_files[0]
        pca_dirpath, pca_fname = rnaseq_pca_file['root'], rnaseq_pca_file['fn']
        pca_fpath = join(pca_dirpath, pca_fname)
        pca_data, color_by_sample, color_by_cond, variances = parse_pca_data(pca_fpath)
        pca_data = self.ignore_samples(pca_data)

        description = ("<p>PCA is a popular method that is based on the principles of dimensional reduction. "
            "Below is a PCA plot of the samples within the space of the first two principal components that explain the most variation in the data. "
            "These were calculated using the read counts of the top 1000 most variable genes within the dataset.</p>")

        legend = ''
        if color_by_cond:
            label_style = 'font-family: \'Lucida Grande\', \'Lucida Sans Unicode\', Arial, Helvetica, sans-serif; ' \
                          'font-size: 12px; ' \
                          'font-weight: bold; '
            legend += '<center><div>'
            legend += '<span style="' + label_style + ' margin-right: 10px;">Conditions: </span>'
            for cond, color in color_by_cond.items():
                legend += '<span style="white-space: nowrap;">'
                legend += '<span style="display: inline-block; width: 16px; height: 12px; ' + \
                          '            margin-bottom: -1px; margin-right: 1px; background-color: ' + color + '"></span>'
                legend += '<span style="' + label_style + ' margin-right: 20px; white-space: normal;"> ' + cond + '</span>'
                legend += '</span>'
            legend += '</div></center>'

        self.add_section(
            name='Principal Components Analysis',
            anchor='rnaseq_az-pca',
            content=description + legend + scatter.plot(pca_data, {
                'title': 'Principal Components Analysis',
                'xlab': 'PC1: ' + variances[0] + '% variance',
                'ylab': 'PC2: ' + variances[1] + '% variance',
                'colors': color_by_sample,
                'tt_label': 'PC1: {point.x}<br/>PC2: {point.y}',
            })
        )
