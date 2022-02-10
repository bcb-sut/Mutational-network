import pandas as pd


class CancerPathways:
    def __init__(self, data_loader):
        self.cancer_pathway_genes = data_loader.load_cancer_pathways_genes()

    def find_com_cancer_pathways(self, com):
        com_pathways = []
        cancer_ptw_genes = set([])
        for _, row in self.cancer_pathway_genes.iterrows():
            common_genes = set(com).intersection(set(row['genes']))
            if len(common_genes) > 0:
                com_pathways.append({'pathway': row['pathway'],
                                     'common_genes': common_genes,
                                     'common_genes_cnt': len(common_genes)})
            cancer_ptw_genes = cancer_ptw_genes.union(set(common_genes))
        return pd.DataFrame(com_pathways), cancer_ptw_genes

    def concat_pathway_genes(self, pathway, genes_arr):
        genes = ', '.join(genes_arr)
        return f'{pathway}({genes})'

    def analyze_coms(self, coms):
        print('cancer pathway analysis...')
        res = []
        ptw_dic = {}
        for i in range(len(coms)):
            com = coms[i]
            pathways, cancer_ptw_genes = self.find_com_cancer_pathways(com)

            if len(pathways) == 0:
                pathways_and_genes = ''
            else:
                pathways_and_genes = [self.concat_pathway_genes(pathway, genes)
                                      for pathway, genes in zip(pathways['pathway'], pathways['common_genes'])]
                for ptw in pathways['pathway']:
                    if ptw not in ptw_dic:
                        ptw_dic[ptw] = []
                    ptw_dic[ptw].append(i)

            res.append({'cancer_pathways': pathways_and_genes, 'cancer_pathways_cnt': len(pathways),
                        'cancer_pathways_genes_cnt': len(cancer_ptw_genes),
                        'cancer_ptw_genes_ratio': len(cancer_ptw_genes) / len(com)})

        ptw_df = []
        for ptw in ptw_dic.keys():
            ptw_df.append({'pathway': ptw, 'coms': ptw_dic[ptw], 'coms_cnt': len(ptw_dic[ptw])})

        print('done')
        return pd.DataFrame(res), pd.DataFrame(ptw_df)


if __name__ == '__main__':
    cp = CancerPathways('../../../data/')
    sample = cp.cancer_pathway_genes['genes'].iloc[0]
    print(sample)
    print(type(sample))
