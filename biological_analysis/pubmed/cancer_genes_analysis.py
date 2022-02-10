import pandas as pd
from gene_dic import GeneDic


class CancerGene:
    def __init__(self, data_folder):
        self.cancer_genes = pd.read_excel(data_folder + 'pubmed/cancer_genes.xlsx')
        self.breast_cancer_genes = pd.read_excel(data_folder + 'pubmed/breast_cancer_genes.xlsx')
        self.gene_dic = GeneDic('../../data/')

    def com_cancer_breast_cancer_genes(self, com):
        com_cancer = set(self.cancer_genes['gene']).intersection(set(com))
        com_breast_cancer = set(self.breast_cancer_genes['gene']).intersection(set(com))
        novel_genes = set(com) - com_cancer - com_breast_cancer
        return com_cancer, com_breast_cancer, novel_genes

    def separate_coding_non_coding(self, genes):
        coding, non_coding = [], []
        for gene in genes:
            gene_type = self.gene_dic.get_gene_type(gene)
            if gene_type == 'coding':
                coding.append(gene)
            else:
                non_coding.append(gene)
        return coding, non_coding

    def generate_row(self, genes, prefix):
        cancer, breast, novel = self.com_cancer_breast_cancer_genes(genes)
        row = {f'{prefix}_cancer': cancer, f'{prefix}_cancer_cnt': len(cancer),
               f'{prefix}_breast_cancer': breast, f'{prefix}_breast_cnt': len(breast),
               f'{prefix}_novel': novel, f'{prefix}_novel_cnt': len(novel)}
        return row

    def cancer_genes_analysis(self, coms):
        res = []
        cancer_coms_cnt = 0
        for com in coms:
            coding, non_coding = self.separate_coding_non_coding(com)
            coding_row = self.generate_row(coding, 'coding')
            row = self.generate_row(non_coding, 'non_coding')
            row.update(coding_row)

            cancer_cnt = row['coding_cancer_cnt'] + row['non_coding_cancer_cnt']
            if cancer_cnt > 0:
                cancer_coms_cnt += 1
            row['cancer_ratio'] = cancer_cnt / len(com)
            res.append(row)

        print(f'{cancer_coms_cnt / len(coms) * 100} percent of coms have cancer-related genes')
        return pd.DataFrame(res)
