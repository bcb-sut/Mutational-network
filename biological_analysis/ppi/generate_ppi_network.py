import networkx as nx
import pickle as pkl
from biological_analysis.string.string_names_dic import StringNames


class PpiNetGenerator:
    def __init__(self, data_folder):
        self.data_folder = data_folder + 'stringdb/'
        self.link_file_path = self.data_folder + '9606.protein.links.v11.5.txt'

        self.string_names = StringNames(data_folder)

        self.low_confidence_score = 150
        self.medium_score = 400
        self.high_score = 700

    def get_score_name(self, score):
        if score <= self.low_confidence_score:
            return 'low'
        if score <= self.medium_score:
            return 'medium'
        return 'high'

    def generate_network(self, min_score):
        print(f'generating ppi network with min_score {min_score}...')
        f = open(self.link_file_path)
        line = f.readline()
        g = nx.Graph()
        while line:
            line = f.readline()
            try:
                p1, p2, score = line.split()[:3]
                p1, p2 = self.string_names.translate_gene(p1), self.string_names.translate_gene(p2)
                score = float(score)
                if score >= min_score:
                    g.add_edge(p1, p2, weight=score)
            except:
                print(f'cant split line: --{line}--')
        print(f'done. number of nodes: {len(g.nodes)}, number of edges: {len(g.edges)}')
        print('saving dump...')
        pkl.dump(g, open(self.data_folder + f'interaction_network_{self.get_score_name(min_score)}.pkl', 'wb'))
        print('done')


if __name__ == '__main__':
    netgen = PpiNetGenerator('../../../data/')
    netgen.generate_network(400)
