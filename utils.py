import pandas as pd

import config
from old4 import data_loader_old
import pickle as pkl


class Singleton(type):
    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]


class Dictionary(metaclass=Singleton):
    def __init__(self):
        self.gene_name = data_loader_old.load_gene_name_dic()
        self.gene_prot, self.prot_gene = None, None
        self.gene_chromosome = None
        self.gene_type = None
        self.coding_subtypes = ['protein_coding', 'processed_transcript']
        self.fast_gene_types = {}
        self.fast_subtypes = {}
        print('---gene name dictionary loaded---')

    def get_gene_name(self, gene_id):
        if gene_id in self.gene_name:
            return self.gene_name[gene_id]['name']
        return gene_id

    def get_coms_names(self, coms):
        coms_names = []
        for com in coms:
            coms_names.append(map(self.get_gene_name, com))
        return coms_names

    def get_prot_gene(self, prot):

        if not self.prot_gene:
            print('loading prot-gene dic')
            self.gene_prot, self.prot_gene = data_loader_old.load_gene_prot_dics()

        return get_if_has(self.prot_gene, prot)

    def map_prot_name(self, prot):

        gene_id = self.get_prot_gene(prot)

        gene_name = self.get_gene_name(gene_id)
        return gene_name

    def get_gene_chromosome(self, gene_name):

        self.gene_chromosome, chrom = load_find(self.gene_chromosome, data_loader_old.load_gene_list,
                                                'gene_symbol', gene_name, 'chr')
        return chrom

    def get_gene_subtype(self, gene):
        self.gene_type, subtype = load_find(self.gene_type, data_loader_old.load_gene_list
                                            , 'gene_symbol', gene, 'Gene_class')
        return subtype

    def get_gene_type(self, gene):
        subtype = self.get_gene_subtype(gene)

        if subtype in self.coding_subtypes:
            return 'coding'
        else:
            return 'non-coding'

    def get_subtype_fast(self, gene):
        if len(self.fast_subtypes) == 0:
            self.fast_subtypes = data_loader_old.load_obj(config.fast_subtype_name)
        return self.fast_subtypes[gene]

    def get_gene_type_fast(self, gene):
        if len(self.fast_gene_types) == 0:
            self.fast_gene_types = data_loader_old.load_obj(config.fast_gene_type_name)
        return self.fast_gene_types[gene]

    def generate_network_dic(self, net, node_function):
        nodes = net.nodes
        res_dic = {}
        for node in nodes:
            res_dic[node] = node_function(node)
        return res_dic

    def generate_subtype_dic(self, net):
        return self.generate_network_dic(net, self.get_gene_subtype)

    def generate_node_type_dic(self, net):
        return self.generate_network_dic(net, self.get_gene_type)


def add_dic_set(dic, k, v, score):
    if k not in dic:
        dic[k] = {}
    dic[k][v] = score


def get_if_has(dic, k):
    try:
        return dic[k]
    except:
        return 'not-found-' + k


def find_com(coms, gene):
    for com in coms:
        if gene in com:
            return com
    return 'not found'


def load_find(ref_df, load_func, ref_col, target, target_col):
    if type(ref_df) != type(pd.DataFrame([])):
        print(f'loading file...')
        ref_df = load_func()
        print('---done---')

    try:
        res = ref_df[ref_df[ref_col] == f'"{target}"'][target_col].iloc[0]
    except:
        print(f'not found')
    return ref_df, res


def save_obj(obj, name):
    pkl.dump(obj, open(config.dump_path + f'{name}.pkl', 'wb'))


if __name__ == '__main__':
    dic = Dictionary()
    net = data_loader_old.load_name_network(0.15)
    node_types = dic.generate_subtype_dic(net)
    save_obj(node_types, 'net0.15-node-subtypes')
