import pandas as pd
import requests  ## python -m pip install requests
import json


class StringApi:
    def __init__(self):
        string_api_url = "https://version-11-5.string-db.org/api"
        output_format = "json"
        method = "enrichment"
        self.request_url = "/".join([string_api_url, output_format, method])

    def get_enrichment(self, genes, max_fdr=0.05):
        print(f'getting string enrichment for gene-set with size:{len(genes)}...')

        params = {
            "identifiers": "%0d".join(genes),  # your protein
            "species": 9606,  # species NCBI identifier
            "caller_identity": "www.awesome_app.org"  # your app name
        }
        response = requests.post(self.request_url, data=params)
        data = json.loads(response.text)

        try:
            df = pd.DataFrame(data)
            df['fdr'] = [float(fdr) for fdr in df['fdr']]
            df = df[df['fdr'] <= max_fdr]
            go = df[df['category'].isin(['Component', 'Process', 'Function'])]
            kegg = df[df['category'] == 'KEGG']

            print('done')
            return list(kegg['term']), list(go['term'])
        except:
            print(data)
            return [], []


if __name__ == '__main__':
    # genes = ['FBN2', 'RYR2', 'C2orf88', 'ZFHX4', 'GABRB2', 'RPS6KA2', 'KCNIP1', 'OTX2-AS1', 'RP11-457K10.1', 'MAD1L1', 'EGLN3', 'KLHL29', 'COL23A1', 'MRPS22']
    genes = ['FBN2', 'RYR2', 'C2orf88', 'ZFHX4', 'GABRB2', 'RPS6KA2', 'KCNIP1', 'OTX2-AS1', 'RP11-457K10.1', 'MAD1L1',
             'EGLN3', 'KLHL29', 'COL23A1', 'MRPS22']  # string_name = StringNames('../../../data/')
    # genes = [string_name.translate_gene(g) for g in genes]
    # print(genes)
    api = StringApi()
    kegg, go = api.get_enrichment(genes)
    print(len(kegg))
    print(kegg.head())
    # kegg.to_csv('kegg.csv')
    # go.to_csv('go.csv')
