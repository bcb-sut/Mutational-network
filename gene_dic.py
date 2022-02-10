import pandas as pd
import pickle as pkl


class GeneDic:
    def __init__(self, data_folder):
        self.data_folder = data_folder
        self.dic = pkl.load(open(data_folder + 'dump/gene_dic.pkl', 'rb'))

    def generate_dic(self):
        print('generating gene-name dic...')
        df = pd.read_excel(self.data_folder + 'genes_list_hg19.xlsx')
        dic = {}
        name_set = set()
        for _, row in df.iterrows():
            id_ = row['gene_name']
            name = row['gene_symbol'][1:-1]
            gene_subtype = row['Gene_class']
            gene_type = self.get_type(gene_subtype)
            if name in name_set:
                print(f'{name}, {id_}, {gene_type}, {dic[name]}')
                # dic[id_] = {'name_id': name, 'type': gene_type}
                continue
            name_set.add(name)
            dic[id_] = {'name_id': name, 'type': gene_type}
            dic[name] = {'name_id': id_, 'type': gene_type}
        pkl.dump(dic, open(self.data_folder + 'dump/gene_dic.pkl', 'wb'))
        print('done')

    @staticmethod
    def get_type(subtype: str):
        coding_subtypes = ['protein_coding', 'processed_transcript']
        if subtype.lower() in coding_subtypes:
            return 'coding'
        return 'non-coding'

    def get_gene_name(self, gene):
        try:
            return self.dic[gene]['name_id']
        except:
            print(f'{gene} not found')
            return gene

    def get_gene_type(self, gene):
        if gene in self.dic:
            return self.dic[gene]['type']
        return 'not-found'


if __name__ == '__main__':
    dic = GeneDic('../../data/')
    dic.generate_dic()
