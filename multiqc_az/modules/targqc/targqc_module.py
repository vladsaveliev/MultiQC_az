#!/usr/bin/env python

""" MultiQC module to parse output from TargQC """

from __future__ import division
from __future__ import absolute_import
from __future__ import print_function
from collections import defaultdict, OrderedDict
import logging
import os
from os.path import basename, dirname, join, splitext

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import linegraph, table

from ngs_utils.Sample import BaseSample

from .rendering import calc_cell_contents
from .reporting import FullReport, Metric

log = logging.getLogger(__name__.replace('multiqc_az', 'multiqc'))


targqc_repr              = 'TargQC'
targqc_name              = 'targqc'

qualimap_name                   = 'qualimap'
qualimap_report_fname           = 'qualimapReport.html'
qualimap_genome_results_fname   = 'genome_results.txt'
qualimap_raw_data_dirname       = 'raw_data_qualimapReport'

qualimap_ishist_fname           = 'insert_size_histogram.txt'
qualimap_covhist_fname          = 'coverage_histogram.txt'
qualimap_gchist_fname           = 'mapped_reads_gc-content_distribution.txt'

picard_name              = 'picard'
picard_ishist_fname      = 'picard_ishist.txt'

fastqc_name              = 'fastqc'
dedup_bam                = 'dedup'  # work/post_processing/dedup and -ready.dedup.bam

fastqc_repr              = 'FastQC'
fastqc_report_fname      = 'fastqc_report.html'


class Sample(BaseSample):
    def __init__(self, name, dirpath, *args, **kwargs):
        BaseSample.__init__(self, name, dirpath, targqc_dirpath=dirpath, *args, **kwargs)
        self.targqc_norm_depth_vcf_txt   = None
        self.targqc_norm_depth_vcf_tsv   = None

        self.targqc_txt_fpath            = join(self.targqc_dirpath, 'summary.txt')
        self.targqc_html_fpath           = join(self.targqc_dirpath, 'summary.html')
        self.targqc_json_fpath           = join(self.targqc_dirpath, 'summary.json')
        self.targqc_region_txt           = join(self.targqc_dirpath, 'regions.txt')
        self.targqc_region_tsv           = join(self.targqc_dirpath, 'regions.tsv')

        self.qualimap_dirpath = join(self.targqc_dirpath, 'qualimap')
        self.qualimap_html_fpath            = join(self.qualimap_dirpath, qualimap_report_fname)
        self.qualimap_genome_results_fpath  = join(self.qualimap_dirpath, qualimap_report_fname)
        self.qualimap_raw_dirpath           = join(self.qualimap_dirpath, qualimap_raw_data_dirname)
        self.qualimap_ins_size_hist_fpath   = join(self.qualimap_raw_dirpath, qualimap_ishist_fname)
        self.qualimap_cov_hist_fpath        = join(self.qualimap_raw_dirpath, qualimap_covhist_fname)
        self.qualimap_gc_hist_fpath         = join(self.qualimap_raw_dirpath, qualimap_gchist_fname)

        self.picard_dirpath                 = join(self.targqc_dirpath, picard_name)
        self.picard_ins_size_hist_txt_fpath = join(self.picard_dirpath, picard_ishist_fname)
        self.picard_ins_size_hist_pdf_fpath = join(self.picard_dirpath, splitext(picard_ishist_fname)[0] + '.pdf')


class OldStyleSample(BaseSample):
    def __init__(self, name, dirpath, *args, **kwargs):
        BaseSample.__init__(self, name, dirpath, targqc_dirpath=dirpath, *args, **kwargs)
        self.targqc_norm_depth_vcf_txt   = None
        self.targqc_norm_depth_vcf_tsv   = None

        self.targqc_txt_fpath         = join(self.targqc_dirpath, name + '.targetSeq.txt')
        self.targqc_html_fpath        = join(self.targqc_dirpath, name + '.targetSeq.html')
        self.targqc_json_fpath        = join(self.targqc_dirpath, name + '.targetSeq.json')
        self.targqc_detailed_txt      = join(self.targqc_dirpath, name + '.targetSeq.details.gene.txt')
        self.targqc_detailed_tsv      = join(self.targqc_dirpath, name + '.targetSeq.details.gene.tsv')

        self.qualimap_dirpath = join(self.targqc_dirpath, 'qualimap')
        self.qualimap_html_fpath            = join(self.qualimap_dirpath, qualimap_report_fname)
        self.qualimap_genome_results_fpath  = join(self.qualimap_dirpath, qualimap_report_fname)
        self.qualimap_raw_dirpath           = join(self.qualimap_dirpath, qualimap_raw_data_dirname)
        self.qualimap_ins_size_hist_fpath   = join(self.qualimap_raw_dirpath, qualimap_ishist_fname)
        self.qualimap_cov_hist_fpath        = join(self.qualimap_raw_dirpath, qualimap_covhist_fname)
        self.qualimap_gc_hist_fpath         = join(self.qualimap_raw_dirpath, qualimap_gchist_fname)

        self.picard_dirpath                 = join(self.targqc_dirpath, picard_name)
        self.picard_ins_size_hist_txt_fpath = join(self.picard_dirpath, picard_ishist_fname)
        self.picard_ins_size_hist_pdf_fpath = join(self.picard_dirpath, splitext(picard_ishist_fname)[0] + '.pdf')


class MultiqcModule(BaseMultiqcModule):
    # noinspection PyMethodMayBeStatic
    def get_s_name(self, f):
        return f['s_name']

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='TargQC', anchor='targqc',
            href="https://github.com/vladsaveliev/TargQC",
            info="- Alignment target coverage analysis tool")

        # self.general_stats_headers = OrderedDict()
        # self.general_stats_data = defaultdict(lambda: dict())
        self.coverage_hist = OrderedDict()
        self.insert_size_hist = OrderedDict()
        self.gc_content_dist = OrderedDict()

        old_style = False
        samples = []
        for f in self.find_log_files('targqc/summary'):
            if f:
                f['s_name'] = basename(f['root'])
                if f['s_name'] not in [s.name for s in samples]:
                    s = Sample(f['s_name'], f['root'], bam=None)
                    samples.append(s)
                    self.add_data_source(f, section='targqc_summary')

        # Exit if we didn't find anything
        if not samples:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        for module, value_dict in {
            'FastQC': {
                'percent_gc': False,
                'percent_duplicates': False,
                'total_sequences': False,
            }}.items():
            if module not in config.table_columns_visible:
                config.table_columns_visible[module] = dict()
            config.table_columns_visible[module].update(value_dict)
            config.skip_generalstats = True

        # Set up the ge2neral stats table
        targqc_full_report = _make_tarqc_html_report(samples, None, tag_by_sample=None)

        # Common metrics
        self.add_common_info(targqc_full_report, old_style)

        # # Add to the General Stats table (has to be called once per MultiQC module)
        self.add_section(**self.make_table(targqc_full_report))
        # general_stats_data = defaultdict(lambda: dict())
        # for s_rep in targqc_full_report.sample_reports:
        #     for rec in s_rep.records:
        #         if rec.value is not None:
        #             general_stats_data[s_rep.sample.name][rec.metric.name] = rec.value
        # general_stats_headers = make_general_stats_headers(targqc_full_report)
        # general_stats_addcols(general_stats_data, general_stats_headers)

        # Make the plots for the report
        for s in samples:
            self.coverage_hist[s.name] = parse_coverage_hist(s.qualimap_cov_hist_fpath)
            self.add_data_source(s_name=s.name, source=s.qualimap_cov_hist_fpath, section='coverage_histogram')

            self.gc_content_dist[s.name], self.human_gc_dist = parse_gc_content_hist(s.qualimap_gc_hist_fpath)
            self.add_data_source(s_name=s.name, source=s.qualimap_gc_hist_fpath, section='mapped_gc_distribution')

            self.insert_size_hist[s.name] = parse_insert_size_hist(s.qualimap_ins_size_hist_fpath)
            self.add_data_source(s_name=s.name, source=s.qualimap_ins_size_hist_fpath, section='insert_size_histogram')

        self.make_report_sections()

    @staticmethod
    def make_table(targqc_full_report):
        data = defaultdict(lambda: dict())
        for s_rep in targqc_full_report.sample_reports:
            for rec in s_rep.records:
                if rec.value is not None:
                    data[s_rep.sample.name][rec.metric.name] = rec.value

        headers = make_general_stats_headers(targqc_full_report)

        return {'name': 'Target coverage statistics',
                'anchor': 'targqc',
                'plot': table.plot(data, headers)}

    def add_common_info(self, targqc_full_report, old_style=False):
            # Metric('Target',                  short_name='Target',                            common=True),
            # Metric('Reference size',          short_name='Reference bp', unit='bp',           common=True),
            # Metric('Regions in target',       short_name='Regions in target',                 common=True),
            # Metric('Bases in target',         short_name='Target bp', unit='bp',              common=True),
            # Metric('Percentage of reference', short_name='Percentage of reference', unit='%', common=True),
            # Metric('Genes in target',         short_name='Genes in target',                   common=True),
        recs = targqc_full_report.get_common_records()

        def _find_rec(_metric_name):
            res = next((r for r in recs if r.metric.name == _metric_name), None)
            assert res is not None, 'Record of metric "' + _metric_name + '" is not found in ' + str([r.metric.name for r in recs])
            return res.format()

        if _find_rec('Scope' if not old_style else 'Target') == 'WGS':
            self.intro += (
                '<b>Target: </b>' + _find_rec('Target') + '<br>' +
                '<b>Genome size: </b>' + _find_rec('Reference size') + '<br>')
        else:
            self.intro += (
                '<b>Target: </b>' + _find_rec('Target') + '<br>' +
                '<b>Target size: </b>' + _find_rec('Bases in target') + ', ' + _find_rec('Percentage of reference') + ' <b>of the reference genome</b>.<br>' +
                '<b>Regions in target: </b>' + _find_rec('Regions in target') + ', <b>genes in target: </b>' + _find_rec('Genes in target') + '.')

    def make_report_sections(self):
        if self.coverage_hist:
            # Chew back on histogram to prevent long flat tail
            # (find a sensible max x - lose 1% of longest tail)
            max_x = 0
            total_bases_by_sample = dict()
            for s_name, d in self.coverage_hist.items():
                total_bases_by_sample[s_name] = sum(d.values())
                cumulative = 0.0
                for count in sorted(d.keys(), reverse=True):
                    cumulative += d[count]
                    if cumulative / total_bases_by_sample[s_name] > 0.01:
                        max_x = max(max_x, count)
                        break

            rates_within_threshs = OrderedDict()
            for s_name, hist in self.coverage_hist.items():
                total = total_bases_by_sample[s_name]
                rates_within_threshs[s_name] = _calc_bases_within_threshs(hist, total, range(max_x + 1))

            self.sections.append({
                'name': 'Cumulative coverage genome fraction',
                'anchor': 'qualimap-genome-fraction-coverage',
                'plot': linegraph.plot(rates_within_threshs, {
                    'title': 'Genome fraction covered by at least X reads (clipped at 1% of the longest tail)',
                    'ylab': 'Fraction of reference (%)',
                    'xlab': 'Coverage (X)',
                    'ymax': 100,
                    'ymin': 0,
                    'xmin': 0,
                    'xmax': max_x,
                    'tt_label': '<b>{point.x}X</b>: {point.y:.2f}%',
                })
            })
            self.sections.append({
                'name': 'Coverage histogram',
                'anchor': 'qualimap-coverage-histogram',
                'plot': linegraph.plot(self.coverage_hist, {
                    'title': 'Coverage histogram (clipped at 1% of cumulative fraction of the longest tail)',
                    'ylab': 'Number of genomic locations',
                    'xlab': 'Coverage (X)',
                    'ymin': 0,
                    'xmin': 0,
                    'xmax': max_x,
                    'xDecimals': False,
                    'tt_label': '<b>{point.x}X</b>: {point.y}',
                })
            })

        # Insert size histogram
        if self.insert_size_hist:
            self.sections.append({
                'name': 'Insert size histogram',
                'anchor': 'qualimap-insert-size-histogram',
                'plot': linegraph.plot(self.insert_size_hist, {
                    'title': 'Insert size histogram',
                    'ylab': 'Fraction of reads',
                    'xlab': 'Insert Size (bp)',
                    'ymin': 0,
                    'xmin': 0,
                    'tt_label': '<b>{point.x} bp</b>: {point.y}',
                })
            })

        # GC-content distribution
        if self.gc_content_dist:
            if self.human_gc_dist:
                self.gc_content_dist['Human reference'] = self.human_gc_dist

            self.sections.append({
                'name': 'GC content distribution',
                'anchor': 'qualimap-gc-distribution',
                'plot': linegraph.plot(self.gc_content_dist, {
                    'title': 'GC content distribution',
                    'ylab': 'Fraction of reads',
                    'xlab': 'GC content (%)',
                    'ymin': 0,
                    'xmin': 0,
                    'xmax': 100,
                    'tt_label': '<b>{point.x}%</b>: {point.y:.3f}',
                })
            })

        # # Duplication Rate Histogram
        # if len(self.qualimap_bamqc_dup_rate_hist) > 0:
        #     self.sections.append({
        #         'name': 'Duplication Rate Histogram',
        #         'anchor': 'dup-rate-hist',
        #         'plot': linegraph.plot(self.qualimap_bamqc_dup_rate_hist, {
        #             'title': 'Duplication Rate Histogram',
        #             'ylab': 'Number of genomic locations',
        #             'xlab': 'Duplication rate (%)',
        #             'ymin': 0,
        #             'xmin': 0,
        #             'xmax': 100,
        #             'tt_label': '<b>{point.x}%</b>: {point.y}',
        #         })
        #     })



def parse_insert_size_hist(qualimap_insert_size_fpath):
    d = dict()
    zero_insertsize = 0
    with open(qualimap_insert_size_fpath) as f:
        for l in f:
            if not l.strip() or l.startswith('#'):
                continue
            insertsize, count = l.split(None, 1)
            insertsize = int(round(float(insertsize)))
            count = float(count) / 1000000
            if insertsize == 0:
                zero_insertsize = count
            else:
                d[insertsize] = count
    return d

def parse_gc_content_hist(qualimap_gc_content_fpath):
    d = dict()
    human_d = dict()
    with open(qualimap_gc_content_fpath) as f:
        for l in f:
            if l.startswith('#'):
                continue
            sections = l.split(None, 3)
            gc = int(round(float(sections[0])))
            content = float(sections[1])
            human_content = float(sections[2])
            d[gc] = content
            human_d[gc] = human_content
    return d, human_d


def parse_coverage_hist(qualimap_coverage_hist_fpath):
    d = dict()
    with open(qualimap_coverage_hist_fpath) as f:
        for l in f:
            if l.startswith('#'):
                continue
            coverage, count = l.split(None, 1)
            coverage = int(round(float(coverage)))
            count = float(count)
            d[coverage] = count

        if len(d) == 0:
            log.debug("Couldn't parse contents of coverage histogram file {}".format(f['fn']))
            return None
    return d


def make_general_stats_headers(report):
    for section in report.metric_storage.sections:
        rows = report.get_rows_of_records(sections=[section])
        calc_cell_contents(report, rows)

    data_by_metric = OrderedDict()

    ms = sorted(report.metric_storage.get_metrics(skip_general_section=True), key=lambda m: m.multiqc.get('order', 99999))
    for metric in ms:
        d = {
            'title': metric.get_display_name(),
            'description': metric.description or metric.name,
            'suffix': metric.unit,
            # 'scale': 'Blues',
            # 'shared_key': 'read_count',
        }
        d.update(metric.multiqc)

        if any(isinstance(val, float) for val in metric.values):
            precision = 2
            for i in range(10, 1, -1):
                if abs(metric.min) < 1./(10**i):
                    precision = i + 1
                    break
            d['format'] = '{:.' + str(precision) + 'f}'

            if '%' in metric.unit:
                d['modify'] = lambda x: x * 100

        elif any(isinstance(val, int) for val in metric.values):
            if metric.max > 100000:
                d['modify'] = lambda x: float(x) / 1000000
                d['format'] = '{:.1f}<span class="hs"></span>M'
            else:
                d['format'] = '{:.0f}'

        if metric.numbers and metric.unit:
            if '%' in metric.unit:
                d['format'] += '<span class="hs"></span>'
            else:
                d['format'] += '<span class="rhs">&nbsp;</span>'
            d['format'] += metric.unit

        if metric.bottom is not None: d['min'] = metric.bottom
        elif '%' in metric.unit:      d['min'] = 0
        if metric.top is not None:    d['top'] = metric.top
        elif '%' in metric.unit:      d['max'] = 100

        if metric.multiqc.get('kind') == 'reads':
            d['scale'] = 'Blues'
        elif metric.multiqc.get('kind') == 'cov':
            d['scale'] = 'Greens'
        elif metric.multiqc.get('kind') == 'trg':
            d['scale'] = 'Oranges'

        data_by_metric[metric.name] = d
    return data_by_metric

def _make_tarqc_html_report(samples, output_dir, tag_by_sample=None):
    jsons_by_sample = {s.name: s.targqc_json_fpath for s in samples if os.path.isfile(s.targqc_json_fpath)}
    htmls_by_sample = {s.name: s.targqc_html_fpath for s in samples if os.path.isfile(s.targqc_html_fpath)}

    if not jsons_by_sample or not htmls_by_sample:
        return None

    targqc_full_report = FullReport.construct_from_sample_report_jsons(samples, output_dir, jsons_by_sample, htmls_by_sample)

    for sample_report in targqc_full_report.sample_reports:
        if tag_by_sample:
            sample_report.set_project_tag(tag_by_sample[sample_report.sample.name])

    return targqc_full_report


def _calc_bases_within_threshs(bases_by_depth, total_size, depth_thresholds):
    bases_within_threshs = OrderedDict((depth, 0) for depth in depth_thresholds)
    rates_within_threshs = OrderedDict((depth, None) for depth in depth_thresholds)

    for depth, bases in bases_by_depth.items():
        for t in depth_thresholds:
            if depth >= t:
                bases_within_threshs[t] += bases
    for t in depth_thresholds:
        bs = bases_within_threshs[t]
        if total_size > 0:
            rate = 100.0 * bases_within_threshs[t] / total_size
            assert rate <= 100, 'Error: rate is > 1: rate = ' + str(rate) + ', bases = ' + str(bs) + ', size = ' + str(total_size)
            rates_within_threshs[t] = rate

    return rates_within_threshs
