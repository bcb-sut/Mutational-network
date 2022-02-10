import matplotlib.pyplot as plt
import pandas as pd
from plotly.subplots import make_subplots

import community_analysis
import config
from old4 import data_loader_old
import plotly_template
from utils import Dictionary


class Figures:
    # def __init__(self):
    # self.coms = data_loader.load_coms(0.15)
    # self.patient_gene = data_loader.load_patient_gene()
    # self.patients = list(self.patient_gene.keys())
    def count_covering_coms(self, patient):
        cnt = 0
        for com in self.coms:
            if community_analysis.com_cover_patient(com, self.patient_gene[patient]):
                cnt += 1
        return cnt

    def count_covering_patients(self, com):
        cnt = 0
        for patient in self.patients:
            if community_analysis.com_cover_patient(com, self.patient_gene[patient]):
                cnt += 1
        return cnt

    def sample_coverages(self):
        covered_patient_cnts = []
        for com in self.coms:
            covered_patient_cnts.append(self.count_covering_patients(com))

        com_sizes = [len(com) for com in self.coms]
        covered_com_cnt = [self.count_covering_coms(patient) for patient in self.patients]

        plt.title('number of patients covered by each community')
        plt.xlabel('community size')
        plt.ylabel('number of covered patients')
        plt.plot(com_sizes, covered_patient_cnts)
        plt.show()

        plt.plot(covered_com_cnt)
        plt.show()

    def assortativity(self):
        df = data_loader_old.load_assortativity_table()
        self.save_fig(plotly_template.table(df.index, df.values), 'assortativity')

    def degree_dist(self):
        thrs = [0.1, 0.15, 0.2]
        figs = []
        for thr in thrs:
            net = data_loader_old.load_network(thr)
            degrees = [net.degree(node) for node in net.nodes()]
            # degress = [deg / max(deg) for deg in degrees]
            hist = plotly_template.histogram(degrees, 100)
            figs.append(hist)

        fig = make_subplots(rows=1, cols=3,
                            vertical_spacing=0.05, horizontal_spacing=0.1,
                            subplot_titles=[f'threshold : {thr}' for thr in thrs])
        # subplot_titles=['a', 'b', 'c', 'd', 'e'])
        col_idx = 1
        for dist in figs:
            fig.add_trace(dist, row=1, col=col_idx)
            fig.update_xaxes(title_text="Node degree", row=1, col=col_idx)
            # fig.update_yaxes(type='log', row=1, col=col_idx)
            col_idx += 1
        fig.update_yaxes(title_text='Number of nodes (log)', row=1, col=1)

        yrange = dict(range=[0, 3])
        fig.update_yaxes(yrange, type='log')
        fig.update_layout(showlegend=False)
        self.save_fig(fig, 'degree-dist')

    def go_process(self):
        pass

    def go_function(self):
        pass

    def go_component(self):
        pass

    def cancer_pathways_clustermap(self):
        pass

    def get_key_val_list(self, dic):
        return list(dic.keys()), list(dic.values())

    def gene_subtypes_pie(self):
        net = data_loader_old.load_name_network(0.15)
        nodes = net.nodes
        types = [Dictionary().get_gene_type_fast(node) for node in nodes]
        subtypes = [Dictionary().get_subtype_fast(node) for node in nodes]

        df = pd.DataFrame({'gene': nodes, 'type': types, 'subtype': subtypes})
        coding_df = df[df['type'] == 'coding']
        non_coding = df[df['type'] != 'coding']

        coding_subtypes = set(coding_df['subtype'])
        non_coding_subtypes = set(non_coding['subtype'])

        coding_cnt = [list(coding_df['subtype']).count(subtype) for subtype in coding_subtypes]
        non_coding_cnt = [list(non_coding['subtype']).count(subtype) for subtype in non_coding_subtypes]

        self.save_fig(plotly_template.pie(coding_cnt, list(coding_subtypes), colors=['1f77b4', '17becf']),
                      'coding-subtypes-pie')
        # self.save_fig(plotly_template.pie(non_coding_cnt, list(non_coding_subtypes),
        #                                   colors = px.colors.sequential.RdBu), 'non-coding-subtypes-pie')

    def gene_type_piechart(self):
        count_df = {'non-coding': 199, 'coding': 762}
        keys, vals = self.get_key_val_list(count_df)
        self.save_fig(plotly_template.pie(vals, keys), 'type-pie')

    def count_gene_vals(self, genes, func):
        gene_vals = [func(gene) for gene in genes]
        types_set = set(gene_vals)
        count_dic = {}
        for type_ in types_set:
            count_dic[type_] = gene_vals.count(type_)
        return self.get_key_val_list(count_dic)

    def gene_type_table(self):
        count_df = {'non-coding': 199, 'coding': 762, 'total': 961}
        vals = list(count_df.values())
        self.save_fig(plotly_template.table(list(count_df.keys()), vals), 'genes-types-table')

    def link_type_table(self):
        net = data_loader_old.load_name_network(0.15)
        count_dic = {'cc': 0, 'nn': 0, 'cn': 0}
        for u, v in net.edges:
            u_type = Dictionary().get_gene_type_fast(u)
            v_type = Dictionary().get_gene_type_fast(v)
            if u_type != v_type:
                count_dic['cn'] += 1
            elif u_type == 'coding':
                count_dic['cc'] += 1
            else:
                count_dic['nn'] += 1
        count_dic['total'] = len(net.edges)
        key = list(count_dic.keys())
        vals = list(count_dic.values())
        self.save_fig(plotly_template.table(key, vals), 'edge-types')

    def ppi_vs_com_links(self):
        pass

    def coms_ppi_enrichment(self):
        pass

    def cis_trans(self):
        pass

    def save_fig(self, fig, name):
        fig.write_image(config.figure_path + f'python/svg/{name}.svg')
        fig.write_image(config.figure_path + f'python/png/{name}.png')
        print(f'figure {name} saved.')

    # table = fig.assortativity()
    # table.write_image('../figures/table.svg')


if __name__ == '__main__':
    print(config.data_path)
    config.data_path = '../../data/'
    fig = Figures()

    fig.degree_dist()
