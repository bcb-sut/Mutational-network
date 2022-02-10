import pickle as pkl

import networkx as nx
import pandas as pd
from gene_dic import GeneDic
from copy import deepcopy
import ast

class DataLoader:
    def __init__(self, data_folder, res_folder):
        self.data_folder = data_folder
        self.res_folder = res_folder
        self.dic = GeneDic(data_folder)

    def load_network(self, thr):
        net = pkl.load(open(self.data_folder + f'/networks/network{thr}.pkl', 'rb'))
        return nx.relabel_nodes(net, self.dic.get_gene_name)
        # return net

    def load_com(self, thr):
        return pkl.load(open(self.res_folder + f'coms/com{thr}.pkl', 'rb'))

    def save_com_table(self, thr, coms):
        df = pd.DataFrame([])

        df['size'] = [len(com) for com in coms]
        df['genes'] = coms
        df.to_excel(self.res_folder + f'/coms/com{thr}.xlsx')

    def save_com_dump(self, thr, coms):
        pkl.dump(coms, open(self.res_folder + f'/coms/com{thr}.pkl', 'wb'))

    def save_com(self, thr, coms):
        self.save_com_table(thr, coms)
        self.save_com_dump(thr, coms)

    def load_gene_patient(self):
        print('loading gene-patient')
        dic = pkl.load(open(self.data_folder + '/dump/gene_patient.pkl', 'rb'))
        print('done')
        return dic

    def load_patient_gene(self):
        print('loading patient-gene')
        dic = pkl.load(open(self.data_folder + '/dump/patient_gene.pkl', 'rb'))
        print('done')
        return dic

    def load_cancer_pathways_genes(self):
        df = pd.read_excel(self.data_folder + 'kegg/cancer_pathways_genes.xlsx')
        df['genes'] = [ast.literal_eval(genes) for genes in df['genes']]
        return df

if __name__ == '__main__':
    d = DataLoader('../../data/', '')
    d.load_gene_patient()
