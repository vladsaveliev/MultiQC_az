#!/usr/bin/env python

""" MultiQC module to add link to Bcl2fastq reports """

import plotly as py
import plotly.graph_objs as go
from os.path import join, dirname, abspath
import pandas as pd
import numpy as np


from multiqc.modules.base_module import BaseMultiqcModule

from multiqc.plots import table, scatter, heatmap
from collections import OrderedDict


class MultiqcModule(BaseMultiqcModule):

    def addVolcano(self, de):

        data = []
        buttons = []
        i=0
        first_title=''
        for c in de:
            vis=False
            if i==0:
                first_title = c
                vis = True

            low_p = de[c].loc[de[c]['p']<=0.1]
            hi_p = de[c].loc[de[c]['p'] > 0.1]

            data.append(go.Scatter(x=low_p['lfc'].tolist(),
                                   y=low_p['p'].tolist(),
                                   mode = 'markers',
                                   text=de[c]['gene'].tolist(),
                                   name='p value < 0.1',
                                   marker={'color':'#33CFA5'},
                                   visible=vis))

            data.append(go.Scatter(x=hi_p['lfc'].tolist(),
                                   y=hi_p['p'].tolist(),
                                   mode = 'markers',
                                   text=de[c]['gene'].tolist(),
                                   name='p value > 0.1',
                                   marker={'color': '#F06A6A'},
                                   visible=vis))

            vis_list = [False] * 2 * len(de)
            vis_list[i * 2]=True
            vis_list[i * 2 + 1] = True

            i+=1

            buttons.append({'label':c, 'method':'update', 'args':[{'visible':vis_list},{'title':c}]})

        updatemenus = [{'type':"buttons",
                        'active':-1,
                        'buttons':buttons}]



        layout = go.Layout(xaxis=dict(title='Log2 fold change'),
                           yaxis=dict(title='p value'),
                           height=700,
                           updatemenus=updatemenus,
                           title=first_title)
        fig = go.Figure(data=data, layout=layout)

        tab_content = py.offline.plot(fig, auto_open=False, output_type='div', include_plotlyjs=False)

        self.add_section(
            name='Volcano',
            anchor='Volcano',
            content=tab_content
        )

    def addNumDE_perContrast(self, de):
        de_num = {}
        p_val = 0.1
        logFold = 2
        for cont in de:
            curr_de = de[cont]
            num = len(curr_de.loc[(curr_de['p'] > 1) & (abs(curr_de['lfc']) > 0.5)])
            de_num[cont] = {'numbef_of_de_genes': num}


        print(de_num)

        headers = OrderedDict()
        headers['Sample Name'] = {'title': 'Contrast'}
        headers['numbef_of_de_genes'] = {'title': 'number of DE genes'}

        self.add_section(
            name='DEs per contrast',
            anchor='DEs percontrast',

            content=table.plot(de_num, headers, {'col1_header': 'Contrast'})
        )

        return de_num

    def addDE_overlap(self, de):

        de_overlap = {}

        for i in de:
            de_overlap[i] = {}
            for j in de:
                de_i = de[i].loc[(de[i]['p'] > 1) & (abs(de[i]['lfc']) > 0.5)]
                de_j = de[j].loc[(de[j]['p'] > 1) & (abs(de[j]['lfc']) > 0.5)]
                i_de_genes = set(de_i.index.tolist())
                j_de_genes = set(de_j.index.tolist())

                de_overlap[i][j] = len(i_de_genes & j_de_genes)

        de_overlap = pd.DataFrame(de_overlap)

        hmdata = de_overlap.values.tolist()
        names = de_overlap.index.tolist()

        hm_html = heatmap.plot(hmdata, names)

        self.add_section(
            name='DE overlap',
            anchor='DE overlap',
            content=hm_html,
            description='Table of numbers of overlapping genes across contrasts'
        )

    def addTopGenes(self, de):

        tab_header = '<ul>'
        tab_content = '<div>'

        for c in de:
            print(de[c].head())
            print(list(de[c]))
            de[c] = de[c].drop('gene_id', axis=1)

            de[c].sort_values(by='p', inplace=True, ascending=False)
            top = de[c][0:20]
            # top = top.set_index('gene')

            tab_header += '<li>' + c + '</li>'


            tab_content += '<div>' + top.to_html(columns=['gene', 'p', 'lfc', 'baseMean'],classes = 'de_top', index=False) + '</div>'

        tab_header += '</ul>'
        tab_content += '</div>'

        html_table = '<div class="tabs">' + tab_header + tab_content + '</div>'

        style_file = open(join(dirname(abspath(__file__)), "style.txt"), 'r')
        style = style_file.read()

        script_file = open(join(dirname(abspath(__file__)), "js.txt"), 'r')
        script = script_file.read()

        self.add_section(
            name='Top DE genes',
            anchor='Top DE genes',
            content=style + '\n' + script + '\n' + html_table + '\n' + '<script> $(document).ready(function(){ $(".tabs").lightTabs(); }); </script>',
        )

    def addBaseMeanPlot(self, de):

        data = []
        buttons = []
        i=0
        first_title=''
        for c in de:
            vis=False
            if i==0:
                first_title = c
                vis = True

            data.append(go.Scatter(x=de[c]['baseMean'].tolist(),
                                   y=de[c]['lfc'].tolist(),
                                   mode = 'markers',
                                   text=de[c]['gene'].tolist(),
                                   name='shrunken',
                                   marker={'color':'#33CFA5'},
                                   visible=vis))

            data.append(go.Scatter(x=de[c]['baseMean'].tolist(),
                                   y=de[c]['lfc_un'].tolist(),
                                   mode = 'markers',
                                   text=de[c]['gene'].tolist(),
                                   name='unshrunken',
                                   marker={'color': '#F06A6A'},
                                   visible=vis))

            vis_list = [False] * 2 * len(de)
            vis_list[i * 2]=True
            vis_list[i * 2 + 1] = True

            i+=1

            buttons.append({'label':c, 'method':'update', 'args':[{'visible':vis_list},{'title':c}]})

        updatemenus = [{'type':"buttons",
                        'active':-1,
                        'buttons':buttons}]



        layout = go.Layout(xaxis=dict(type='log', autorange=True, title='Mean counts'),
                           yaxis=dict(title='Log2 Fold change'),
                           height=700,
                           updatemenus=updatemenus,
                           title=first_title)
        fig = go.Figure(data=data, layout=layout)

        tab_content = py.offline.plot(fig, auto_open=False, output_type='div', include_plotlyjs=False)

        self.add_section(
            name='MA plot',
            anchor='MA plot',
            content=tab_content
        )

    def __init__(self):
        mod_name = 'RNAseqDE'
        super(MultiqcModule, self).__init__(name='RNA Differential Expression', anchor=mod_name)
        # make dict of de-tables per contrast
        de = {}
        for f in self.find_log_files(mod_name, filecontents=False):
            print(f)
            dirpath, fname = f['root'], f['fn']
            if f['s_name'] == 'de_gene_key':
                de_path = join(dirpath, fname)
                contrast = f['root'].split('/')
                data = pd.read_csv(de_path)
                de[contrast[-1]] = data

        if len(de)==0:
            print('no DE files')
        else:
            self.addVolcano(de)
            self.addNumDE_perContrast(de)
            self.addDE_overlap(de)
            self.addTopGenes(de)
            self.addBaseMeanPlot(de)



        #
        #
        #
        # def addHeatMap(HM_data):
        #
        #     gene_names = HM_data['aux_sorted$dataHUGO'].tolist()
        #     HM_data.drop(columns=['Unnamed: 0', 'aux_sorted$dataHUGO'], inplace=True)
        #
        #     val = HM_data.values
        #
        #     lis = val.tolist()
        #     names = list(HM_data.columns.values)
        #
        #     pconfig = {
        #         'xTitle': 'Sample Name',
        #         'yTitle': 'Gene id'
        #     }
        #
        #     hm_html = heatmap.plot(lis, names, gene_names, pconfig)
        #     self.add_section(
        #         name='Heatmap',
        #         anchor='Heatmap',
        #         content=hm_html,
        #         description='This plot shows counts of differentially expressed genes on a per-sample basis. We have scaled the data by row'
        #     )
        #
        # # addVolcanoPlot(data)

        # addBaseMeanPlot(data)
        # addTopGenes(data)
        #
        # # HM_data = pd.read_csv(hm_path)
        # # addHeatMap(HM_data)