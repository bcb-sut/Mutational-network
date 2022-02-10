from biological_analysis.pubmed.pubmed_api import PubmedApi
from old4 import data_loader_old
import pandas as pd

if __name__ == '__main__':
    data_folder = ''

    network = data_loader_old.load_name_network(0.15)
    all_genes = list(network.nodes)

    pubmed_api = PubmedApi()

    cancer_df = []
    breast_cancer_df = []
    for gene in all_genes:
        cancer, breast_cancer = pubmed_api.is_gene_cancer_breast_cancer_related(gene)
        if len(cancer) > 0:
            cancer_df.append({'gene': gene, 'pmids': cancer})

        if len(breast_cancer) > 0:
            breast_cancer_df.append({'gene': gene, 'pmids': breast_cancer})

    pd.DataFrame(cancer_df).to_excel(data_folder + 'pubmed/cancer_genes.xlsx')
    pd.DataFrame(breast_cancer_df).to_excel(data_folder + 'pubmed/breast_cancer_genes.xlsx')
