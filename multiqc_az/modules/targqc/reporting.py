# coding: utf-8

from __future__ import absolute_import

import os
from collections import OrderedDict
from os.path import join, relpath, dirname, abspath, basename
from json import load, dump


class Record:
    def __init__(
            self,
            metric=None,
            value=None,
            meta=None,
            html_fpath=None,
            url=None,
            parse=True,
            sort_as=None,
            id=None,

            num=None,
            cell_contents=None,
            frac_width=None,
            right_shift=None,
            color=None,
            text_color=None,
            show_content=True,
            rowspan=None):  # TODO: get rid of those

        self.metric = metric
        if self.metric.parse:
            self.set_value(value)
        else:
            self.value = value
        self.meta = meta or dict()
        self.html_fpath = html_fpath
        self.url = url
        self.id = id

        self.num = num
        # if sort_as is not None:
        #     self.num = sort_as
        #     self.

        self.cell_contents = None
        self.frac_width = None
        self.right_shift = None
        self.color = color
        self.text_color = text_color
        self.show_content = show_content
        self.rowspan = rowspan
        # self.color = lambda: self._color
        # self.text_color = lambda: self._text_color

    def get_value(self):
        return self.value

    def set_value(self, value):
        if value is None:
            pass
        else:
            if isinstance(value, str):
                try:
                    value = int(value)
                except ValueError:
                    try:
                        value = float(value)
                    except ValueError:
                        pass
        self.value = value

    def del_value(self):
        del self.value

    value = property(get_value, set_value, del_value, 'value')

    def format(self, human_readable=True):
        return self.metric.format(self.value, human_readable=human_readable)

    def format_html(self):
        if not self.show_content:
            return ''
        if self.html_fpath:
            if isinstance(self.html_fpath, dict):
                val = ''
                if self.value:
                    val += self.value + ': '
                val += ', '.join(
                    '<a href="' + v + '">' + self.metric.format(k, human_readable=True, is_html=True) + '</a>'
                        for k, v in self.html_fpath.items())
            else:
                val = '<a href="' + self.html_fpath + '">' + self.value + '</a>'
        else:
            val = self.metric.format(self.value, human_readable=True, is_html=True)
        return val

    def __repr__(self):
        return self.metric.name + ' ' + self.metric.format(self.value, human_readable=True)

    @staticmethod
    def load(data):
        data['metric'] = Metric.load(data['metric'])
        return Record(**data)


class Metric:
    skip_when_loading = [
        'numbers',
        'values',
        'min',
        'max',
        'med',
        'low_outer_fence',
        'low_inner_fence',
        'top_inner_fence',
        'top_outer_fence',
        'top_outer_fence',
        'all_values_equal'
    ]

    def __init__(
            self,
            name=None,
            short_name=None,
            description=None,
            quality='More is better',  # "More is better", "Less is better", "Equal"
            unit='',
            common=False,
            ok_threshold=None,
            bottom=None,
            top=None,
            is_hidden=False,
            is_mandatory=False,
            with_heatmap=True,
            style='',
            td_style='',
            class_='',
            td_class='',
            align=None,
            width=None,
            max_width=None,
            min_width=None,
            sort_direction=None,
            parse=True,

            sort_by=None,  # legacy
            header_length=1,

            numbers=None,
            values=None,
            min=None,
            max=None,
            med=None,
            low_outer_fence=None,
            low_inner_fence=None,
            top_inner_fence=None,
            top_outer_fence=None,
            all_values_equal=False,

            section_name=None,
            multiqc=None
        ):

        self.name = name
        self.short_name = short_name
        self.description = description
        self.quality = quality
        self.common = common
        self.unit = unit
        self.ok_threshold = ok_threshold
        self.bottom = bottom
        self.top = top
        self.is_hidden = is_hidden
        self.with_heatmap = with_heatmap
        self.style = style
        self.td_style = td_style
        self.class_ = class_
        self.td_class = td_class
        self.align = align
        self.max_width = max_width
        self.min_width = min_width
        self.sort_direction = sort_direction
        self.parse = parse
        self.header_length = header_length

        self.numbers = []
        self.values = []
        self.min = None
        self.max = None
        self.med = med
        self.low_outer_fence = low_inner_fence
        self.low_inner_fence = low_inner_fence
        self.top_inner_fence = top_inner_fence
        self.top_outer_fence = top_outer_fence
        self.all_values_equal = False

        self.section_name = section_name
        self.multiqc = multiqc or dict()

    def get_display_name(self):
        if self.short_name is not None:
            n = self.short_name
        else:
            n = self.name
        return n

    def format(self, value, human_readable=True, is_html=True):
        return Metric.format_value(value, unit=self.unit, human_readable=human_readable, is_html=is_html)

    @staticmethod
    def format_value(value, unit='', human_readable=True, is_html=True):
        if value is None:
            return ''

        unit_str = unit
        if unit and is_html:
            unit_str = '<span class=\'rhs\'>&nbsp;</span>' + unit

        if isinstance(value, str):
            if human_readable:
                return '{value}{unit_str}'.format(**locals())
            else:
                return value

        elif isinstance(value, int):
            if value == 0:
                return '0'
            if human_readable:
                if value <= 9999:
                    return str(value)
                else:
                    v = '{value:,}{unit_str}'.format(**locals())
                    if is_html:
                        v = v.replace(',', '<span style="margin-left: .2em"></span>')
                    return v
            else:
                return str(value)

        elif isinstance(value, float):
            if value == 0.0:
                return '0'
            if human_readable:
                if unit == '%':
                    value *= 100
                precision = 2
                for i in range(10, 1, -1):
                    if abs(value) < 1./(10**i):
                        precision = i + 1
                        break
                return '{value:.{precision}f}{unit_str}'.format(**locals())
            else:
                return str(value)

        if isinstance(value, list):
            return ', '.join(Metric.format_value(v, unit, human_readable, is_html) for v in value)

        return '.'

    def __repr__(self):
        return self.name

    @staticmethod
    def load(data):
        data = {k: v for k, v in data.items() if k not in Metric.skip_when_loading}
        return Metric(**data)


class BaseReport:
    def __init__(self, sample=None, html_fpath=None, url=None, json_fpath=None,
                 records=None, plots=None, metric_storage=None, display_name='',
                 report_name='', caller_tag=None, project_tag=None, expandable=False,
                 unique=False, keep_order=False, heatmap_by_rows=False, vertical_sample_names=False, **kwargs):
        self.sample = sample
        self.html_fpath = html_fpath
        self.plots = plots or []  # TODO: make real JS plots, not just included PNG
        self.json_fpath = json_fpath
        self.url = url
        self.records = [r for r in records if r.metric] if records else []
        self.metric_storage = metric_storage

        self.report_name = report_name
        self.display_name = display_name
        if not display_name:
            if sample:
                if isinstance(sample, str):
                    self.display_name = sample
                else:
                    self.display_name = sample.name
        self.caller_tag = None
        self.set_caller_tag(caller_tag)
        self.project_tag = None
        self.set_project_tag(project_tag)

        self.expandable = expandable

        self.unique = unique

        self.keep_order = keep_order
        self.heatmap_by_rows = heatmap_by_rows
        self.vertical_sample_names = vertical_sample_names

    def set_project_tag(self, tag):
        if not self.project_tag and tag:
            self.display_name = '[' + tag + ']' + ' ' + self.display_name
            self.project_tag = tag

    def set_caller_tag(self, tag):
        if not self.caller_tag and tag:
            self.display_name = self.display_name + ' ' + tag
            self.project_tag = tag

    def flatten(self, sections=None, human_readable=True):
        raise NotImplementedError()

    def get_rows_of_records(self, sections=list()):  # TODO: move logic from flatten here, use this method both in flatten and save_html
        raise NotImplementedError()

    @staticmethod
    def find_record(records, metric_name):
        try:
            rec = next(r for r in records if r.metric.name == metric_name)
        except StopIteration:
            return None  # if no record for the metric
        except AttributeError:
            return None  # record with no metric
        else:
            return rec

    def get_common_records(self):
        raise NotImplementedError()


class Row:
    def __init__(self, parent_report, records=None, highlighted=False, color=None, class_='', hidden=False):
        self.__parent_report = parent_report
        self.records = records or []
        self.highlighted = highlighted
        self.color = color
        self.class_ = class_
        self.hidden = hidden

        self.numbers = []
        self.min = None
        self.max = None
        self.med = None
        self.low_outer_fence = None
        self.low_inner_fence = None
        self.top_inner_fence = None
        self.top_outer_fence = None
        self.all_values_equal = False

    def add_record(self, metric_name, value, **kwargs):
        metric = self.__parent_report.metric_storage.find_metric(metric_name)
        assert metric, metric_name
        rec = Record(metric=metric, value=value, **kwargs)
        self.records.append(rec)
        return rec


class SampleReport(BaseReport):
    def __init__(self, *args, **kwargs):
        BaseReport.__init__(self, *args, **kwargs)

    def get_common_records(self):
        common_records = []
        for record in self.records:
            if record.metric.common:
                common_records.append(record)
        return common_records

    def get_rows_of_records(self, sections=None):  # TODO: move logic from flatten here, use this method both in flatten and save_html
        if sections:
            recs = []
            for metric in self.metric_storage.get_metrics(sections=sections, skip_general_section=True):
                if not metric.name == 'Sample':
                    rec = BaseReport.find_record(self.records, metric.name)
                    if rec:
                        recs.append(rec)
                    else:
                        recs.append(Record(metric=metric, value=None))
        else:
            recs = self.records
        return [Row(parent_report=self, records=recs)]

    def add_record(self, metric_name, value, silent=False, **kwargs):
        metric = self.metric_storage.find_metric(metric_name.strip())
        if not metric:
            # err('Could not find metric ' + metric_name)
            return None
        rec = Record(metric, value, **kwargs)
        self.records.append(rec)
        # if not silent:
        #     info(metric_name + ': ' + rec.format(human_readable=True))
        return rec

    def flatten(self, sections=None, human_readable=True):
        header_rows = []
        for m in self.metric_storage.general_section.metrics:
            r = BaseReport.find_record(self.records, m.name)
            if r:
                if human_readable:
                    header_rows.append(['## ' + m.name + ': ' + r.format(human_readable=True)])
                else:
                    header_rows.append(['##' + m.name + '=' + r.format(human_readable=False)])

        flat_rows = [['Sample', self.display_name]]
        for metric in self.metric_storage.get_metrics(sections, skip_general_section=True):
            vals = [metric.name]
            rec = BaseReport.find_record(self.records, metric.name)
            vals.append(rec.format(human_readable=human_readable) if rec is not None else '.')
            flat_rows.append(vals)
        return header_rows, flat_rows

    def __repr__(self):
        return self.display_name + (', ' + self.report_name if self.report_name else '')

    @staticmethod
    def load(data, sample=None, bcbio_structure=None):
        data['sample'] = sample or Sample.load(data['sample'], bcbio_structure)
        data['records'] = [Record.load(d) for d in data['records']]
        data['metric_storage'] = MetricStorage.load(data['metric_storage'])

        rep = SampleReport(**data)
        for rec in rep.records:
            rec.metric = rep.metric_storage.find_metric(rec.metric.name)
        return rep


class PerRegionSampleReport(SampleReport):
    def get_common_records(self):
        return []

    def __init__(self, **kwargs):
        SampleReport.__init__(self, **kwargs)
        self.records = []
        self.rows = []

    def add_row(self):
        row = Row(parent_report=self)
        self.rows.append(row)
        return row

    def flatten(self, sections=None, human_readable=True):
        header_rows = []
        for m in self.metric_storage.general_section.metrics:
            rec = BaseReport.find_record(self.records, m.name)
            if rec:
                if human_readable:
                    header_rows.append(['## ' + m.name + ': ' + rec.format(human_readable=True)])
                else:
                    header_rows.append(['##' + m.name + '=' + rec.format(human_readable=False)])

        flat_rows = []
        header_row = []
        metrics = self.metric_storage.get_metrics(sections, skip_general_section=True)
        for i, m in enumerate(metrics):
            header_row.append(('#' if i == 0 else '') + m.get_display_name())
        flat_rows.append(header_row)

        for row in self.rows:
            flat_row = []
            for m in self.metric_storage.get_metrics(sections, skip_general_section=True):
                rec = BaseReport.find_record(row.records, m.name)
                if rec:
                    flat_row.append(rec.format(human_readable=human_readable))
                else:
                    pass
            flat_rows.append(flat_row)

        # next_list_value = next((r.value for r in self.records if isinstance(r.value, list)), None)
        # if next_list_value:
        #     for i in range(len(next_list_value)):
        #         header_row = []
        #         for m in self.metric_storage.get_metrics(sections):
        #             try:
        #                 r = next((r for r in self.records if r.metric.name == m.name), None)
        #                 if r:
        #                     if m.name in self.metric_storage.general_section.metrics_by_name:
        #                         val = r.value
        #                     elif not r.value:
        #                         val = None
        #                     else:
        #                         val = r.value[i]
        #                     header_row.append(r.metric.format(val))
        #             except StopIteration:
        #                 header_row.append('-')  # if no record for the metric
        #         rows.append(header_row)

        return header_rows, flat_rows

    def get_rows_of_records(self, sections=None):  # TODO: move logic from flatten here, use this method both in flatten and save_html
        if not sections:
            sections = []

        elif isinstance(sections, ReportSection):
            sections = [sections]

        rows = []
        for i, row in enumerate(self.rows):
            recs = []
            for metric in self.metric_storage.get_metrics(sections=sections, skip_general_section=True):
                rec = BaseReport.find_record(row.records, metric.name)
                if rec:
                    recs.append(rec)
                else:
                    recs.append(Record(metric=metric, value=None))
            row.records = recs
            rows.append(row)
        return rows

    def save_html(self, output_fpath, caption='', #type_=None,
                  extra_js_fpaths=None, extra_css_fpaths=None,
                  tmpl_fpath=None, data_dict=None, debug=False):
        return None
        # sample_reports = []
        # fr = FullReport(self.report_name, sample_reports, self.metric_storage)

        # for i in range(len(next(r for r in self.records if isinstance(r.value, list)).value)):row = []
        #     row = []
        #     for m in self.metric_storage.get_metrics(self.metric_storage.sections):
        #         try:
        #             r = next((r for r in self.records if r.metric.name == m.name), None)
        #             if r:
        #                 if m.name in self.metric_storage.general_section.metrics_by_name:
        #                     val = r.value
        #                 elif not r.value:
        #                     val = None
        #                 else:
        #                     val = r.value[i]
        #                 row.append(r.metric.format(val))
        #         except StopIteration:
        #             row.append('-')  # if no record for the metric

            # records = []
            # sr = SampleReport(self.sample, self.html_fpath, records=None, metric_storage=None,
            #      report_name='', plots=None, json_fpath=None,
            #      )
            # return Report.save_html(self, output_dirpath, base_fname, caption=caption, type_='SquareSampleReport')

    # def add_record(self, metric_name, value, meta=None):
    #     raise NotImplementedError
    #
    # def add_row(self, row):
    #     # self.metric_storage.get_metric(metric_name.strip())
    #     assert metric, metric_name
    #     rec = Record(metric, value, meta)
    #     self.records.append(rec)
    #     info(metric_name + ': ' + rec.format())
    #     return rec


class FullReport(BaseReport):
    def __init__(self, name='',
                 sample_reports=None,
                 metric_storage=None,
                 general_records=None,
                 base_fname=None, **kwargs):
        BaseReport.__init__(self, **kwargs)
        self.name = name
        self.sample_reports = sample_reports or []
        self.metric_storage = metric_storage
        self._general_records = general_records
        self.base_fname = base_fname

        if metric_storage:
            for sample_report in sample_reports:
                sample_report.metric_storage = metric_storage

        elif sample_reports and sample_reports[0].metric_storage:
            self.metric_storage = sample_reports[0].metric_storage
            for sample_report in sample_reports:
                sample_report.metric_storage = metric_storage

        self.sample_metric = Metric(name='Sample', with_heatmap=False, align='left', parse=False)
        if self.metric_storage:
            for section in self.metric_storage.sections:
                # section.get_metrics = lambda: [self.sample_metric] + section.metrics
                section.add_metric(self.sample_metric, prepend=True)

    def get_common_records(self):
        common_records = []

        if self._general_records:
            common_records.extend(self._general_records)

        if self.sample_reports:
            sample_report = self.sample_reports[0]
            for record in sample_report.records:
                if record.metric and record.metric.common:  #TODO: why can record.metric be None?
                    common_records.append(record)
        return common_records

    def flatten(self, sections=None, human_readable=True):
        if len(self.sample_reports) == 0:
            # err('No sample reports found: summary will not be produced.')
            return []

        header_rows = []
        some_rep = self.sample_reports[0]
        for m in self.metric_storage.general_section.metrics:
            rec = BaseReport.find_record(some_rep.records, m.name)
            if rec:
                if human_readable:
                    header_rows.append(['## ' + m.name + ': ' + rec.format(human_readable=True)])
                else:
                    header_rows.append(['##' + m.name + '=' + rec.format(human_readable=False)])

        rows = [['Sample'] + [rep.display_name for rep in self.sample_reports]]
        # rows_of_records = self.get_rows_of_records(sections) # TODO: use logic from get_rows_of_records
        for metric in self.metric_storage.get_metrics(sections, skip_general_section=True):
            if not metric.name == 'Sample':
                row = [metric.name]
                for sr in self.sample_reports:
                    rec = BaseReport.find_record(sr.records, metric.name)
                    row.append(rec.format(human_readable=human_readable) if rec is not None else '.')
                rows.append(row)
        return header_rows, rows

    def get_rows_of_records(self, sections=None):  # TODO: move logic from flatten here, use this method both in flatten and save_html
        if not sections:
            sections = []

        elif isinstance(sections, ReportSection):
            sections = [sections]

        rows = []
        for i, sr in enumerate(self.sample_reports):
            recs = []
            recs.append(Record(metric=self.sample_metric, value=sr.display_name,
                url=sr.url, html_fpath=sr.html_fpath, num=len(self.sample_reports) - i))

            for metric in self.metric_storage.get_metrics(sections=sections, skip_general_section=True):
                if not metric.name == 'Sample':
                    rec = BaseReport.find_record(sr.records, metric.name)
                    if rec:
                        recs.append(rec)
                    else:
                        recs.append(Record(metric=metric, value=None))
            rows.append(Row(parent_report=self, records=recs))
        return rows

    def get_values_by_metric(self, metric):
        recs = []
        for sr in self.sample_reports:
            for record in sr.records:
                if record.metric.name == metric.name:
                    recs.append(record)
        return recs

    @staticmethod
    def _sync_sections(dst_section, src_section):
        for metric in src_section.metrics:
            if not dst_section.find_metric(metric.name):
                dst_section.add_metric(metric)

    @staticmethod
    def construct_from_sample_report_jsons(samples, output_dirpath,
            jsons_by_sample, htmls_by_sample, bcbio_structure=None):
        ms = None
        sample_reports = []
        for sample in samples:
            if sample.name in jsons_by_sample:
                with open(jsons_by_sample[sample.name]) as f:
                    data = load(f, object_pairs_hook=OrderedDict)
                    sr = SampleReport.load(data, sample, bcbio_structure)
                    sr.html_fpath = htmls_by_sample.get(sample.name)
                    sample_reports.append(sr)
                    if ms is None:
                        ms = sr.metric_storage
                    else:  # Maximize metric storage
                        FullReport._sync_sections(ms.general_section, sr.metric_storage.general_section)
                        for section in sr.metric_storage.sections:
                            if section.name not in ms.sections_by_name:
                                ms.add_section(section)
                            else:
                                FullReport._sync_sections(ms.sections_by_name[section.name], section)

        for sr in sample_reports:
            sr.metric_storage = ms
            for rec in sr.records:
                rec.metric = ms.find_metric(rec.metric.name)
            sr.records = [r for r in sr.records if r.metric is not None]

        fr = FullReport(sample_reports=sample_reports, metric_storage=ms)

        return fr

    def __repr__(self):
        return self.name


def parse_tsv(tsv_fpath):
    values = []

    with open(tsv_fpath, 'r') as f:
        for line in f:
            values.append(line.split('\t'))

    if not values:
        return None
        # critical('Data not found in ' + tsv_fpath)

    return values


class ReportSection:
    def __init__(self, name='', title='', metrics=None, **kwargs):
        self.name = name
        self.title = title
        self.metrics = metrics or []
        self.metrics_by_name = dict((m.name, m) for m in metrics)
        for m in metrics:
            m.section_name = name

    def add_metric(self, metric, prepend=False):
        if not prepend:
            self.metrics.append(metric)
        else:
            self.metrics = [metric] + self.metrics
        self.metrics_by_name[metric.name] = metric
        metric.section_name = self.name

    def get_metrics(self):
        return self.metrics

    def find_metric(self, metric_name):
        return self.metrics_by_name.get(metric_name, None)
        # return next((m for m in self.get_metrics() if m.name == metric_name), None)

    def remove_metric(self, metric_name):
        metric = self.metrics_by_name[metric_name]
        self.metrics.remove(metric)
        del self.metrics_by_name[metric_name]

    def __repr__(self):
        return self.name + ', ' + self.title

    @staticmethod
    def load(data):
        data['metrics'] = [Metric.load(m_data) for m_data in data['metrics']]
        return ReportSection(**data)


# class FullReportSection(ReportSection):
#     def __init__(self, samples, *args, **kwargs):
#         ReportSection.__init__(self, *args, **kwargs)


class MetricStorage:
    def __init__(self,
                 metrics_list=None,
                 general_section=None,
                 sections=None,
                 sections_by_name=None,
                 **kwargs):
        self.sections_by_name = OrderedDict()
        self.sections = []
        self.general_section = general_section or ReportSection('general_section', '', [])
        self.general_section.name = self.general_section.name or 'general_section'

        if sections:
            self.sections = sections
            for section in sections:
                self.sections_by_name[section.name] = section

        elif sections_by_name:
            for name, section in sections_by_name.items():
                self.sections_by_name[name] = section
                self.sections.append(section)

        elif metrics_list:
            section = ReportSection('metrics_list', '', metrics_list)
            self.sections_by_name[''] = section
            self.sections.append(section)

    def add_section(self, section):
        self.sections.append(section)
        self.sections_by_name[section.name] = section

    def add_metric(self, metric, section_name=''):
        self.sections_by_name[section_name].add_metric(metric)

    def find_metric(self, metric_name):
        return next((m for m in self.get_metrics() if m.name == metric_name), None)

    def remove_metric(self, metric_name):
        for section in self.sections:
            section.remove_metric(metric_name)

    def get_metrics(self, sections=None, skip_general_section=False):
        metrics = []
        for section in ((sections or self.sections) + ([] if skip_general_section else [self.general_section])):
            if sections:
                pass
            for m in section.metrics:
                metrics.append(m)
        return metrics

    def __repr__(self, *args, **kwargs):
        return self.sections_by_name.__repr__(*args, **kwargs)

    @staticmethod
    def load(data):
        return MetricStorage(
            sections=[ReportSection.load(d) for d in data['sections']],
            general_section=ReportSection.load(data.get('general_section') or data.get('common_for_all_samples_section')))


def read_sample_names(sample_fpath):
    sample_names = []

    with open(sample_fpath) as f:
        for line in f:
            sample_name = line.strip()
            sample_names.append(sample_name)

    return sample_names


def load_records(json_fpath):
    with open(json_fpath) as f:
        data = load(f, object_pairs_hook=OrderedDict)
        return [Record.load(d) for d in data['records']]

