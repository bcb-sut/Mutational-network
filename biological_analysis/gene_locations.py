import pandas as pd


class GeneLoc:
    def __init__(self, data_folder):
        loc_df = pd.read_excel(data_folder + 'genes_list_hg19.xlsx')
        self.gene_loc_dic = {}
        print('generating gene-loc dic...')
        for _, row in loc_df.iterrows():
            gene = row['gene_symbol'][1:-1]
            chr = row['chr']
            start = int(row['start'])
            end = int(row['end'])
            strand = row['strand']
            gene_id = row['gene_name']

            dic_row = {'chr': chr, 'start': start, 'end': end,
                       'mid': (start + end) / 2, 'strand': strand}

            self.gene_loc_dic[gene] = dic_row
            self.gene_loc_dic[gene_id] = dic_row
        print('done')

    def distance(self, g1, g2):
        if g1 not in self.gene_loc_dic:
            print(f'{g1} not found in gene-loc')
            return -1
        if g2 not in self.gene_loc_dic:
            print(f'{g2} not found in gene-loc')
            return -1
        loc1 = self.gene_loc_dic[g1]['mid']
        loc2 = self.gene_loc_dic[g2]['mid']

        return abs(loc2 - loc1)

    def check_dic_contain(self, gene):
        if gene not in self.gene_loc_dic:
            print(f'{gene} not found in gene-loc dic')
            return False
        return True

    def analyse_cis_trans_distance_strand(self, com_net):
        cis = []
        trans = []
        distance_sum = 0
        same_strand = []
        all_chromes = set([])

        for g1, g2 in com_net.edges:
            if not (self.check_dic_contain(g1)):
                continue
            if not (self.check_dic_contain(g2)):
                continue

            chr1 = self.gene_loc_dic[g1]['chr']
            chr2 = self.gene_loc_dic[g2]['chr']
            all_chromes.add(chr1)
            all_chromes.add(chr2)
            weight = com_net[g1][g2]['weight']

            if chr1 == chr2:
                dist = self.distance(g1, g2)
                cis.append((g1, g2, chr1, dist, weight))
                distance_sum += dist
            else:
                trans.append((f'{g1}({chr1})', f'{g2}({chr2})', weight))

            if self.gene_loc_dic[g1]['strand'] == self.gene_loc_dic[g2]['strand']:
                same_strand.append((g1, g2))
        if len(cis) > 0:
            avg_dist = distance_sum / len(cis)
        else:
            avg_dist = 0
        return cis, trans, avg_dist, same_strand, all_chromes

    def analyze_coms(self, coms, network):
        cis_trans_dist_strand = {'cis': [], 'trans': [], 'dist': [], 'strand': []}
        res_df = pd.DataFrame([])
        all_chromes = []
        for com in coms:
            com_net = network.subgraph(com)
            res = self.analyse_cis_trans_distance_strand(com_net)
            cnt = 0
            for key in cis_trans_dist_strand.keys():
                cis_trans_dist_strand[key].append(res[cnt])
                cnt += 1
            all_chromes.append(res[4])
        res_df['avg_cis_distance'] = cis_trans_dist_strand['dist']
        for key in ['cis', 'trans', 'strand']:
            res_df[f'{key}_links'] = cis_trans_dist_strand[key]
            res_df[f'{key}_cnt'] = [len(arr) for arr in cis_trans_dist_strand[key]]
        res_df['all_chromes'] = all_chromes
        res_df['all_chromes_cnt'] = [len(chromes) for chromes in all_chromes]

        cis_rows = res_df[res_df['cis_cnt'] > 0]
        print(f'avg cis links distance: {cis_rows["avg_cis_distance"].sum() / len(cis_rows)}')
        return res_df
