from __future__ import absolute_import

from math import floor
from .reporting import BaseReport


# Color heat map
BLUE_HUE = 240
BLUE_OUTER_BRT = 55
BLUE_INNER_BRT = 65

GREEN_HUE = 120
GREEN_OUTER_BRT = 50
GREEN_INNER_BRT = 60

RED_HUE = 0
RED_OUTER_BRT = 50
RED_INNER_BRT = 60

MIN_NORMAL_BRT = 80
MEDIAN_BRT = 100  # just white.

def hsl2rgb(h, s, l):
    r, g, b = None, None, None

    if s == 0:
        r = g = b = l  # achromatic
    else:
        q = l * (1 + s) if l < 0.5 else l + s - l * s
        p = 2 * l - q
        r = hue2rgb(p, q, h + 1./3)
        g = hue2rgb(p, q, h)
        b = hue2rgb(p, q, h - 1./3)

    return map(int, [round(r * 255), round(g * 255), round(b * 255)])

def hue2rgb(p, q, t):
    if t < 0: t += 1
    if t > 1: t -= 1
    if t < 1./6: return p + (q - p) * 6 * t
    if t < 1./2: return q
    if t < 2./3: return p + (q - p) * (2./3 - t) * 6
    return p

def get_color(hue, lightness):
    lightness = lightness or 92
    # lightness = Math.round (Math.pow hue - 75, 2) // 350 + 35
    rgb = hsl2rgb(float(hue) / 360, 0.8, float(lightness) / 100)
    hex_rgb = [hex(c)[2:] for c in rgb]
    return '#' + ''.join(hex_rgb)

def mean(values):
    return float(sum(values)) / len(values) if len(values) > 0 else float('nan')

def calc_heatmap_stats(metric):
    numbers = sorted([v for v in metric.numbers],
                     key=lambda a: a if a is not None else -1)  # None is always less than anything
    l = len(numbers)

    if 'Percentage of genome within' in metric.name:
        pass

    metric.min = numbers[0]
    metric.max = numbers[l - 1]
    metric.all_values_equal = metric.min == metric.max
    if metric.med is None:
        metric.med = numbers[(l - 1) // 2] if l % 2 != 0 else mean([numbers[l // 2], numbers[(l // 2) - 1]])
    q1 = numbers[int(floor((l - 1) // 4))]
    q3 = numbers[int(floor((l - 1) * 3 // 4))]

    d = q3 - q1
    metric.low_outer_fence = metric.low_outer_fence if metric.low_outer_fence is not None else q1 - 3   * d
    metric.low_inner_fence = metric.low_inner_fence if metric.low_inner_fence is not None else q1 - 1.5 * d
    metric.top_inner_fence = metric.top_inner_fence if metric.top_inner_fence is not None else q3 + 1.5 * d
    metric.top_outer_fence = metric.top_outer_fence if metric.top_outer_fence is not None else q3 + 3   * d
    return metric

def calc_cell_contents(report, rows):
    if not rows:
        return report

    metrics = report.metric_storage.get_metrics()

    max_frac_widths_by_metric = dict()

    # First round: calculatings max/min integral/fractional widths (for decimal alingment) and max/min values (for heatmaps)
    for row in rows:
        for rec in row.records:
            # if rec.metric.name in section.metrics_by_name:
            _calc_record_cell_contents(rec)

        for rec in row.records:
            if 'of mean depth' in rec.metric.name:
                pass
            # if rec.metric.name in section.metrics_by_name:
            if rec.metric.name not in max_frac_widths_by_metric or \
                            rec.frac_width > max_frac_widths_by_metric[rec.metric.name]:
                max_frac_widths_by_metric[rec.metric.name] = rec.frac_width
            if rec.num is not None:
                row.numbers.append(rec.num)
                rec.metric.numbers.append(rec.num)
            # elif rec.sort_as:
            #     rec.metric.numbers.append(rec.sort_as)
            if rec.value is not None or rec.html_fpath is not None:
                rec.metric.values.append(rec.value)
    # else:
    #     for rec in report.records:
    #         if rec.metric.name in section.metrics_by_name:
    #             if rec.num:
    #                 if not rec.metric.values:
    #                     rec.metric.values = []
    #                 rec.metric.values.append(rec.num)
    if report.heatmap_by_rows:
        for row in rows:
            if row.numbers:
                row = calc_heatmap_stats(row)
    else:
        for metric in metrics:
            if metric.numbers and metric.with_heatmap:
                # For metrics where we know the "normal value" - we want to color everything above normal white,
                #   and everything below - red, starting from normal, finishing with bottom
                if metric.ok_threshold is not None:
                    if isinstance(metric.ok_threshold, int) or isinstance(metric.ok_threshold, float):
                        metric.med = metric.ok_threshold
                        if metric.bottom is not None:
                            metric.low_outer_fence = metric.bottom
                metric = calc_heatmap_stats(metric)

    # Second round: setting shift and color properties based on max/min widths and vals
    for row in rows:
        for rec in row.records:
            # Padding based on frac width
            if rec.frac_width:
                rec.right_shift = max_frac_widths_by_metric[rec.metric.name] - rec.frac_width

            metric = rec.metric
            heatmap_stats = rec.metric if not report.heatmap_by_rows else row

            # Color heatmap
            if rec.num is not None and metric.with_heatmap:
                [top_hue, inner_top_brt, outer_top_brt] = [BLUE_HUE, BLUE_INNER_BRT, BLUE_OUTER_BRT]
                [low_hue, inner_low_brt, outer_low_brt] = [RED_HUE, RED_INNER_BRT, RED_OUTER_BRT]

                if metric.quality == 'Less is better':  # then swap colors
                    [top_hue, low_hue] = [low_hue, top_hue]
                    [inner_top_brt, inner_low_brt] = [inner_low_brt, inner_top_brt]
                    [outer_top_brt, outer_low_brt] = [outer_low_brt, outer_top_brt]

                if metric.ok_threshold is not None:
                    if isinstance(metric.ok_threshold, int) or isinstance(metric.ok_threshold, float):
                        if rec.num >= metric.ok_threshold:
                            continue  # white on black

                        # rec_to_align_with = sample_report.find_record(sample_report.records, metric.threshold)
                        # if rec_to_align_with:
                        #     rec.text_color = lambda: rec_to_align_with.text_color()
                        #     continue

                if not heatmap_stats.all_values_equal:
                    rec.text_color = 'black'

                    # Low outliers
                    if rec.num < heatmap_stats.low_outer_fence and rec.num < heatmap_stats.med:
                        rec.color = get_color(low_hue, outer_low_brt)
                        rec.text_color = 'white'

                    elif rec.num < heatmap_stats.low_inner_fence and rec.num < heatmap_stats.med:
                        rec.color = get_color(low_hue, inner_low_brt)

                    # Normal values
                    elif rec.num < heatmap_stats.med:
                        try:
                            k = float(MEDIAN_BRT - MIN_NORMAL_BRT) / (heatmap_stats.med - heatmap_stats.low_inner_fence)
                        except:
                            pass
                        else:
                            brt = round(MEDIAN_BRT - (heatmap_stats.med - rec.num) * k)
                            rec.color = get_color(low_hue, brt)

                    # High outliers
                    elif rec.num > heatmap_stats.top_inner_fence and rec.num > heatmap_stats.med:
                        rec.color = get_color(top_hue, inner_top_brt)

                    elif rec.num > heatmap_stats.top_outer_fence and rec.num > heatmap_stats.med:
                        rec.color = get_color(top_hue, outer_top_brt)
                        rec.text_color = 'white'

                    elif rec.num > heatmap_stats.med:
                        k = float(MEDIAN_BRT - MIN_NORMAL_BRT) / (heatmap_stats.top_inner_fence - heatmap_stats.med)
                        brt = round(MEDIAN_BRT - (rec.num - heatmap_stats.med) * k)
                        rec.color = get_color(top_hue, brt)

        for rec in row.records:
            metric = rec.metric

            if metric.ok_threshold is not None:
                if isinstance(metric.ok_threshold, str):
                    rec_to_align_with = BaseReport.find_record(row.records, metric.ok_threshold)
                    if rec_to_align_with:
                        rec.text_color = rec_to_align_with.text_color
                        rec.color = rec_to_align_with.color
    return report

def _calc_record_cell_contents(rec):
    rec.cell_contents = rec.format_html()
    # if rec.metric.sort_by and rec.value:
    #     try:
    #         rec.num = rec.metric.sort_by(rec.value)
    #     except:
    #         pass
    #     rec.metric.with_heatmap = False

    if rec.value is not None and (isinstance(rec.value, int) or isinstance(rec.value, float)):
        rec.num = rec.value

    #TODO: intPartTextWidth
    return rec
