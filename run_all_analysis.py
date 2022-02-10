import pandas as pd
from data_loader import DataLoader
from community_detection import community_detection
from community_analysis import CommunityAnalysis
from biological_analysis.ppi.ppi_analysis import PpiAnalysis
from biological_analysis.cancer_pathways import CancerPathways
from biological_analysis.gene_locations import GeneLoc
from biological_analysis.pubmed.cancer_genes_analysis import CancerGene
from biological_analysis.string.string_analysis import StringAnalysis

save_path = '../../results/new-run/'


def concat_and_save(all_analysis, df):
    all_analysis_df = pd.concat([all_analysis, df], axis=1)
    all_analysis_df.to_excel(save_path + 'all_analysis.xlsx')
    return all_analysis_df


if __name__ == '__main__':
    data_folder = '../../data/'
    data_loader = DataLoader(data_folder, save_path)

    network_thr = 0.15
    network = data_loader.load_network(network_thr)

    # all_analysis_df = pd.DataFrame([])
    all_analysis_df = pd.read_excel(save_path + 'all_analysis.xlsx')

    # ensemble_size = 100
    # coms = community_detection.hierarchical_ensemble(network, ensemble_size)
    # data_loader.save_com(0.15, coms)
    coms = data_loader.load_com(0.15)

    # all_analysis_df['coms'] = coms
    # all_analysis_df['node cnt'] = [len(com) for com in coms]
    # edges = [list(network.subgraph(com).edges) for com in coms]
    # all_analysis_df['links'] = edges
    # all_analysis_df['link cnt'] = [len(links) for links in edges]
    #
    # com_analysis = CommunityAnalysis(coms, save_path, data_loader)
    # coms_coverages = com_analysis.sample_coverages()
    # all_analysis_df['patient_coverages'] = coms_coverages
    # all_analysis_df.to_excel(save_path + 'all_analysis.xlsx')
    #
    # rand_cnt = 1000
    # ppi_analysis = PpiAnalysis(data_folder, save_path)
    # ppi_df = ppi_analysis.evaluate_coms_ppi_links(network, coms, rand_cnt)
    # all_analysis_df = concat_and_save(all_analysis_df, ppi_df)
    #
    # cancer_gene = CancerGene(data_folder)
    # cancer_gene_analysis = cancer_gene.cancer_genes_analysis(coms)
    # all_analysis_df = concat_and_save(all_analysis_df, cancer_gene_analysis)
    #
    # cancer_pathways = CancerPathways(data_loader)
    # coms_ptw_df, ptw_df = cancer_pathways.analyze_coms(coms)
    # all_analysis_df = concat_and_save(all_analysis_df, coms_ptw_df)
    # ptw_df.to_excel(save_path + 'cancer_pathways_analysis.xlsx')

    gene_loc = GeneLoc(data_folder)
    gene_loc_df = gene_loc.analyze_coms(coms, network)
    all_analysis_df = concat_and_save(all_analysis_df, gene_loc_df)

    # string_analysis = StringAnalysis()
    # string_res_df, go_df, kegg_df = string_analysis.analyze_coms(coms)
    # all_analysis_df = concat_and_save(all_analysis_df, string_res_df)
    # go_df.to_excel(save_path + 'go_analysis.xlsx')
    # kegg_df.to_excel(save_path + 'kegg_analysis.xlsx')
