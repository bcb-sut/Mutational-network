from pymed import PubMed


class PubmedApi:
    def __init__(self):
        self.pubmed = PubMed()

    def is_gene_cancer_breast_cancer_related(self, gene):
        results = self.pubmed.query(f"{gene} + cancer", max_results=100)
        cancer_pmids = []
        breast_cancer_pmids = []
        for res in results:
            title = res.title.lower()
            abstract = res.abstract.lower()
            keywords = [k.lower() for k in res.keywords]
            keywords = ', '.join(keywords)

            concat_all = title + abstract + keywords
            pmid = res.pubmed_id.split('\n')[0]

            if gene.lower() in concat_all:
                if 'cancer' in concat_all:
                    cancer_pmids.append(pmid)
                if 'breast cancer' in concat_all:
                    breast_cancer_pmids.append(pmid)

        return cancer_pmids, breast_cancer_pmids
