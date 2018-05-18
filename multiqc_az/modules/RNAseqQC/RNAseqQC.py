


""" MultiQC module to parse output from bcbioRNASeq Quality control """

from multiqc.modules.base_module import BaseMultiqcModule
import plotly as py
import pandas as pd
import plotly.graph_objs as go

import scipy.stats as st
from multiqc.plots import scatter, heatmap
import logging
from os.path import join


# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        mod_name = 'RNAseqQC'

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='bcbiornaseqqc', anchor=mod_name)
        file_names, roots = [], []
        for f in self.find_log_files(mod_name, filecontents=False):
            print(f)
            dirpath, fname = f['root'], f['fn']
            if f['s_name'] == 'rawCounts':
                raw_counts = pd.read_csv(join(dirpath, fname))
                # print(raw_counts.head())
            if f['s_name'] == 'normalizedCounts':
                norm_counts = pd.read_csv(join(dirpath, fname))
                # print(norm_counts.head())
            if f['s_name'] == 'corMatrix':
                raw_data = pd.read_csv(join(dirpath, fname))
                # print(raw_data.head())
            if f['s_name'] == 'pca':
                pca_data = pd.read_csv(join(dirpath, fname), index_col=[0])
            if f['s_name'] == 'tpm':
                tpm = pd.read_csv(join(dirpath, fname),index_col=[0])
            if f['s_name'] == 'gene2biotype':
                biotype = pd.read_csv(join(dirpath, fname),index_col=[0])



        col_names = list(raw_counts)[1:]
        group_num = len(col_names)

        raw_counts['sum'] = raw_counts.sum(axis=1)

        self.plot_correlation_heatmap(raw_counts, norm_counts, col_names, group_num)
        #self.plot_mean_sd(raw_counts, norm_counts, col_names, group_num, vst, rlog, combined_counts)
        #self.plot_disp_ests(combined_counts, genes_est,genes_final,genes_fitted)
        self.plot_covariates(raw_data)
        self.plot_pca(pca_data)
        self.tpm_perbiotype(tpm,biotype)

    def plot_pca(self, pca_data):

        print(pca_data)
        data =[]
        for name, row in pca_data.iterrows():
            data.append(go.Scatter(x=[row['pc1']],
                                   y=[row['pc2']],
                                   mode = 'markers',
                                   marker={'size': 10},
                                   text=[name],
                                   name=name))

        layout = go.Layout(height=700)



        fig = go.Figure(data=data, layout=layout)

        tab_content = py.offline.plot(fig, auto_open=False, output_type='div', include_plotlyjs=True)

        self.add_section (
            name = 'PCA plot',
            anchor = 'pca',
            description = 'PCA is a popular method that is based on the principles of dimensional reduction. Below is a PCA plot of the samples within the space of the first two principal components that explain the most variation in the data. These were calculated using the read counts of the top 1000 most variable genes within the dataset.',

            content = tab_content
        )

    def plot_covariates(self, raw_data):

        # find number of PCs
        pc_num = raw_data['covar'].value_counts().tolist()
        pc_num = pc_num[0]

        comp_names = []
        for iter, row in raw_data.iterrows():
            comp_names.append(row['compare'])
            if iter % pc_num == 2:
                break

        cor_names = []
        for _, row in raw_data.iterrows():
            if not (row['covar'] in cor_names):
                cor_names.append(row['covar'])

        hmdata = [[None for _ in range(len(cor_names))] for _ in range(3)]
        for iter,row in raw_data.iterrows():
            j = int(iter / pc_num)
            i = iter % pc_num
            if row['fdr'] < 0.1:
                hmdata[i][j] = row['r']

        pconfig = {
            'title': "bcbioRNASeq Quality Control: PCA Covariates",                 # Plot title - should be in format "Module Name: Plot Title"
            'square': False,                # Force the plot to stay square? (Maintain aspect ratio)
            'reverseColors': False,        # Reverse the order of the colour axis
            'decimalPlaces': 2,            # Number of decimal places for tooltip
            'legend': True,                # Colour axis key enabled or not
            'borderWidth': 0,              # Border width between cells
            'datalabels': True,            # Show values in each cell. Defaults True when less than 20 samples.
            'datalabel_colour': '<auto>',  # Colour of text for values. Defaults to auto contrast.
        }

        hm_html = heatmap.plot(hmdata, cor_names, comp_names, pconfig)

        self.add_section (
            name = 'PCA Covariates',
            anchor = 'covar_section',
            description = 'When multiple factors may influence the results of a given experiment, it is useful to assess which of them is responsible for the most variance as determined by PCA. '
                          'We adapted the method described by Daily et al. where they integrated a method to correlate covariates with principal components values to determine the importance of each factor. '
                          'Here we are showing the correlational analysis of the rlog transformed count data’s principal components with the metadata covariates of interest. Significant correlations (FDR < 0.1) are shaded from blue (anti-correlated) to red (correlated), with non-significant correlations set to zero and shaded in gray.',
            #helptext = 'Inter-correlation analysis (ICA) is another way to look at how well samples cluster by plotting the correlation between the expression profiles of the samples. Pearson`s correlation coefficient is a measure of how well your data would be fitted by a linear regression.',
            plot = hm_html
        )

    def plot_correlation_heatmap(self, raw_counts, norm_counts, col_names, group_num):
        norm_counts = norm_counts[raw_counts['sum'] > 0]
        hmdata = [[st.pearsonr(norm_counts[col_names[i]], norm_counts[col_names[j]])[0]
                    for i in range(group_num)] for j in range(group_num)]

        pconfig = {
            'title': "bcbioRNASeq Quality Control: Correlation Heatmap",
            'xTitle': None,                # X-axis title
            'yTitle': None,                # Y-axis title
            'min': None,                   # Minimum value (default: auto)
            'max': None,                   # Maximum value (default: auto)
            'square': False,               # Force the plot to stay square? (Maintain aspect ratio)
            'reverseColors': False,        # Reverse the order of the colour axis
            'decimalPlaces': 2,            # Number of decimal places for tooltip
            'legend': True,                # Colour axis key enabled or not
            'borderWidth': 0,              # Border width between cells
            'datalabels': True,            # Show values in each cell. Defaults True when less than 20 samples.
            'datalabel_colour': '<auto>',  # Colour of text for values. Defaults to auto contrast.
        }

        hm_html = heatmap.plot(hmdata, col_names, col_names, pconfig)

        self.add_section (
            name = 'Correlation Heatmap',
            anchor = 'heatmap_section',
            description = 'This heatmap shows Pearson`s correlation values between groups.',
            helptext = 'Inter-correlation analysis (ICA) is another way to look at how well samples cluster by plotting the correlation between the expression profiles of the samples. Pearson`s correlation coefficient is a measure of how well your data would be fitted by a linear regression.',
            plot = hm_html
        )

    def tpm_perbiotype(self, tpm, biotype):
        samples = list(tpm)

        biotype_dict = biotype.gene_biotype.to_dict()

        biotype_list = []
        for gene in tpm.index.tolist():
            if gene in biotype_dict:
                biotype_list.append(biotype_dict[gene])
            else:
                # print('Warning! ' + gene + ' not in biotype_dict')
                biotype_list.append('NA')

        se = pd.Series(biotype_list)



        tpm.insert(0, 'biotype', se.values)
        tpm = tpm.dropna()

        types = set(tpm.biotype.tolist())
        groupped = tpm.groupby(by='biotype')
        tpm_per_biotype_per_sample = {}

        for type in types:
            g = groupped.get_group(type)
            tpm_per_sample = {}
            for sample in samples:
                tpm_per_sample[sample] = g[sample].tolist()

            tpm_per_biotype_per_sample[type] = tpm_per_sample

        layout = go.Layout(yaxis=dict(type='log', autorange=True, ), showlegend=False, height = 500)
        link = []
        for sample in samples:
            print(sample)
            data = []
            for type in types:
                data.append(go.Box(y=tpm_per_biotype_per_sample[type][sample], name=type))

            fig = go.Figure(data=data, layout=layout)
            fig['layout'].update(title=sample)


            link.append(py.offline.plot(fig, auto_open=False, output_type='div', include_plotlyjs=False))

        html_string=''
        for l in link:
            html_string += '''<div style="width:50%; display:inline-block" frameborder="0" seamless="seamless" scrolling="no">''' + l + '''</div>'''



        self.add_section (
            name = 'TPM per biotype',
            anchor = 'TPM_per_biotype',
            content = html_string
        )