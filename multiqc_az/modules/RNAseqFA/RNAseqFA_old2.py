import xml.etree.ElementTree as ET
import pandas as pd
from os.path import join, dirname, isfile
from collections import OrderedDict
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import table, heatmap


class MultiqcModule(BaseMultiqcModule):

    def pathway_enrichment_heatmap(self, pw):
        data = pd.DataFrame({'A': []})

        contrast_names = []
        for cont in pw:
            data = pd.concat([data, pw[cont]['enrichmentScore']], axis=1)
            contrast_names.append(cont[6:])

        data.drop('A', axis=1, inplace=True)
        data.fillna(0, inplace=True)
        print(data)

        hmdata = data.values.tolist()
        pathway_names = data.index.tolist()

        print(contrast_names)
        print(pathway_names)
        hm_html = heatmap.plot(hmdata, contrast_names, pathway_names)

        self.add_section(
            name='Heatmap of enriched pathways',
            anchor='Heatmap of enriched pathways',
            content=hm_html,
            description='Shows enrichment score for different contrasts'
        )

    def pathway_graphs(self, pw):

        for cont in pw:
            p = pw[cont]

            print(p)

    def __init__(self):

        mod_name = 'RNAseqFA'

        super(MultiqcModule, self).__init__(name='DE genes pathways', anchor=mod_name)
        # make dict of pathway tables per contrast
        pw = {}

        for f in self.find_log_files(mod_name, filecontents=False):
            print(f)
            dirpath, fname = f['root'], f['fn']
            if f['s_name'] == 'pathway_table':
                pw_path = join(dirpath, fname)
                contrast = f['root'].split('/')
                data = pd.read_csv(pw_path, usecols=[0,1,3,6], index_col=[0])
                pw[contrast[-1]] = data
                print(data)

        if len(pw) > 0:
            self.pathway_enrichment_heatmap(pw)
            self.pathway_graphs(pw)
        else:
            print('pw is empty')

        # def addPathWayTable(path_table):
        #
        #
        #     data = path_table.to_dict(orient='index')
        #
        #     headers = OrderedDict()
        #     headers['Description'] = {'title': 'Description'}
        #     headers['setSize'] = {'title': 'Set Size', 'format': '{:,.0f}' }
        #     headers['enrichmentScore'] = {'title': 'Enrichment Score', 'format': '{:,.2f}' }
        #     headers['p.adjust'] = {'title': 'p-value adj.', 'format': '{:,.2f}' }
        #
        #     table_config = {'col1_header': 'Pathway KEGG ID'}
        #
        #     html_content = table.plot(data, headers, table_config)
        #
        #     self.add_section(
        #         name='Pathways of DE-genes',
        #         anchor='PW',
        #         content=html_content,
        #         description='Table shows metabolic pathways, which contains sets of differentially expressed genes'
        #     )
        #
        # addPathWayTable(path_table)
        #
        # def generate_path_graph(pw_index, elem, cy_style, pw_dir):
        #
        #     path_view_gene = pw_dir + pw_index + '_pathway.csv'
        #     if isfile(path_view_gene):
        #         path_genes = pd.read_csv(path_view_gene)
        #     else:
        #         return ''
        #
        #     for g in path_genes.iterrows():
        #         node_id = str(g[1].T['Unnamed: 0'])
        #         x = str(g[1].T['x'])
        #         y = str(g[1].T['y'])
        #         lab = g[1].T['labels']
        #         color = g[1].T['mol.col']
        #         wid = str(g[1].T['width'])
        #         hei = str(g[1].T['height'])
        #
        #         shape = 'rectangle'
        #         tvalign = 'center'
        #         if g[1].T['type'] == 'compound':
        #             shape = 'ellipse'
        #             tvalign = 'top'
        #
        #         elem += '{ ' + \
        #                 'data: {id: "' + pw_index + '_' + node_id + '"},' + \
        #                 'classes: "graph_' + pw_index + '",' + \
        #                 'position:{ x:' + x + ',y:' + y + '}' + \
        #                 '},'
        #
        #         cy_style += '{' \
        #                     'selector: "#' + pw_index + '_' + node_id + '", ' + \
        #                     'style: {' + \
        #                     '"content": "' + lab + '",' + \
        #                     '"width": "' + wid + '",' + \
        #                     '"height": "' + hei + '",' + \
        #                     '"border-width": 0.5,' + \
        #                     '"font-size": 10,' + \
        #                     '"border-style": "solid",' + \
        #                     '"border-color": "#333",' + \
        #                     '"background-color": "' + color + '",' + \
        #                     '"shape": "'+shape+'",' + \
        #                     '"text-valign": "' +tvalign+ '"' + \
        #                     '}},'
        #
        #     xml_path = pw_dir + pw_index + '.xml'
        #     tree = ET.parse(xml_path)
        #     root = tree.getroot()
        #
        #
        #     for el in root:
        #         if el.tag == 'relation':
        #             try:
        #                 el_cont = iter(el).__next__()
        #
        #                 b = pw_index + '_' + el.attrib['entry1']
        #                 e = pw_index + '_' + el.attrib['entry2']
        #
        #                 if el_cont.attrib['name']=='inhibition':
        #                     ln_style ='solid'
        #                     color = '#777'
        #                     arr_shape = 'tee'
        #
        #                 elif el_cont.attrib['name']=='indirect effect':
        #                     ln_style ='dashed'
        #                     color = '#bbb'
        #                     arr_shape = 'triangle'
        #
        #                 else:
        #                     ln_style ='solid'
        #                     color = '#777'
        #                     arr_shape = 'triangle'
        #
        #                 elem += '{ ' + \
        #                         'classes: "relation graph_' + pw_index + '", ' + \
        #                         'data: {id: "' + b + '_' + e + '", source: "' + b + '",' + 'target: "' + e + '"}' + \
        #                         '},'
        #
        #                 cy_style += '{' + \
        #                             'selector: "#' + b + '_' + e + '",' + \
        #                             'style: {' + \
        #                             '"width": 1,' + \
        #                             '"curve-style": "bezier",' + \
        #                             '"line-color": "'+color+'",' + \
        #                             '"target-arrow-color": "'+color+'",' + \
        #                             '"line-style": "'+ln_style+'",' + \
        #                             '"target-arrow-shape": "'+arr_shape+'"' + \
        #                             '}},'
        #             except StopIteration:
        #                 continue
        #
        #         if el.attrib['type'] == 'map':
        #             el_cont = iter(el).__next__()
        #
        #             node_id = el.attrib['id']
        #             x = el_cont.attrib['x']
        #             y = el_cont.attrib['y']
        #             wid = el_cont.attrib['width']
        #             hei = el_cont.attrib['height']
        #             if 'name' in el_cont.attrib:
        #                 lab = el_cont.attrib['name']
        #
        #             if lab[0:5]=='TITLE':
        #                 lab=lab[6:].upper()
        #
        #
        #             elem += '{ ' + \
        #                     'data: {id: "' + pw_index + '_' + node_id + '"},' + \
        #                     'classes: "graph_' + pw_index + '",' + \
        #                     'position:{ x:' + x + ',y:' + y + '}' + \
        #                     '},'
        #
        #             cy_style += '{' \
        #                         'selector: "#' + pw_index + '_' + node_id + '", ' + \
        #                         'style: {' + \
        #                         '"content": "' + lab + '",' + \
        #                         '"width": "' + wid + '",' + \
        #                         '"height": "' + hei + '",' + \
        #                         '"border-width": 1,' + \
        #                         '"font-size": 10,' + \
        #                         '"border-style": "solid",' + \
        #                         '"border-color": "#333",' + \
        #                         '"shape": "roundrectangle",' + \
        #                         '"text-wrap": "wrap",' + \
        #                         '"background-color": "#adf",' + \
        #                         '"text-max-width": '+ str(0.9 * int(wid)) + ',' + \
        #                         '"text-valign": "center"' + \
        #                         '}},'
        #
        #         if el.attrib['type'] == 'group':
        #             el_cont = iter(el).__next__()
        #
        #             node_id = el.attrib['id']
        #             x = el_cont.attrib['x']
        #             y = el_cont.attrib['y']
        #             wid = el_cont.attrib['width']
        #             hei = el_cont.attrib['height']
        #
        #             elem += '{ ' + \
        #                     'data: {id: "' + pw_index + '_' + node_id + '"},' + \
        #                     'classes: "graph_' + pw_index + '",' + \
        #                     'position:{ x:' + x + ',y:' + y + '}' + \
        #                     '},'
        #
        #             cy_style += '{' \
        #                         'selector: "#' + pw_index + '_' + node_id + '", ' + \
        #                         'style: {' + \
        #                         '"width": "' + wid + '",' + \
        #                         '"background-opacity": 0,' + \
        #                         '"border-opacity": 0,' + \
        #                         '"height": "' + hei + '",' + \
        #                         '}},'
        #
        #
        #     return elem, cy_style
        #
        # cy_style = 'style: ['
        # elem = 'elements: ['
        #
        # pw_list = path_table.index.tolist()
        # btnHtml = ''
        # btnFun = ''
        # hidGr = ''
        #
        # print(pw_list)
        #
        # pw_dir = f['root'] + '/'
        #
        # for pw_it in pw_list:
        #     pw = str(pw_it)
        #     elem, cy_style = generate_path_graph(pw, elem, cy_style, pw_dir)
        #     btnHtml += '<button id="' + pw + '" onclick="show_' + pw + '()" class="inact">' + pw + '</button>'
        #     hidGr += 'var hid_' + pw + ' = cy.remove(".graph_' + pw + '");'
        #     btnFun += 'function show_' + pw + '(){ cy.remove(":inside"); hid_' + pw + '.restore(); \
        #     var x = document.getElementsByClassName("act"); var i; for (i = 0; i < x.length; i++)  \
        #       { x[i].className = "inact";} document.getElementById("'+pw+'").className = "act";    };'
        #
        # elem = elem[:-1]
        # cy_style = cy_style[:-1]
        # elem += ']'
        # cy_style += ']'
        #
        # cyto_scape_path = '<script src="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.2.7/cytoscape.js"></script>'
        # div_style = "<style> #cy { width: 100%; height: 900px; position: relative; top: 0px; left: 0px; } </style>"
        # btn_style = '<style>.inact {  font-size: 14px; margin: 2px;   background-color: #e7e7e7;    border: none;    color: black;    padding: 6px 32px;    text-align: center;    text-decoration: none;    display: inline-block;} \
        #              .act {  font-size: 14px;  margin: 2px;  background-color: #4CAF50;    border: none;    color: white;    padding: 6px 32px;    text-align: center;    text-decoration: none;    display: inline-block; }</style>'
        # layout = 'layout: {    name: "preset"  }'
        #
        # open_first_graph = 'hid_' + pw_list[0] + '.restore(); document.getElementById("'+pw_list[0]+'").className = "act";'
        #
        # graph_html = div_style +btn_style+ btnHtml + '<div id="cy"></div>' + cyto_scape_path + \
        #              '<script> var cy = cytoscape({container: document.getElementById("cy"), ' + \
        #              'minZoom: 0.5, maxZoom: 2,' + \
        #              elem + ',' + layout + ',' + cy_style + ' });' + hidGr + btnFun + open_first_graph+'</script>'
        #
        #
        #
        # self.add_section(
        #     name='Pathway detalization',
        #     anchor='Pathview',
        #     content=graph_html,
        # )
