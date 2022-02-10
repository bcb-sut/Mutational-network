from biological_analysis.string.string_api import StringApi
import pandas as pd


def set_col_arr_and_arr_cnt(df, col, col_values):
    df[col] = col_values
    df[f'{col}_cnt'] = [len(arr) for arr in col_values]
    return df


def reverse_coms_analysis(coms_results, term_key):
    all_terms = {}
    for i in range(len(coms_results)):
        for term in coms_results[i]:
            if term not in all_terms:
                all_terms[term] = []
            all_terms[term].append(i)
    res_df = []
    for term in all_terms.keys():
        res_df.append({term_key: term, 'coms': all_terms[term],
                       'coms_cnt': len(all_terms[term])})
    return pd.DataFrame(res_df)


class StringAnalysis:
    def __init__(self):
        self.string_api = StringApi()

    def analyze_coms(self, coms):
        all_kegg, all_go = [], []
        for com in coms:
            kegg, go = self.string_api.get_enrichment(com)
            all_kegg.append(kegg)
            all_go.append(go)
        res_df = pd.DataFrame([])
        res_df = set_col_arr_and_arr_cnt(res_df, 'kegg', all_kegg)
        res_df = set_col_arr_and_arr_cnt(res_df, 'go', all_go)

        go_res_analysis = reverse_coms_analysis(all_go, 'go')
        kegg_res_analysis = reverse_coms_analysis(all_kegg, 'kegg')

        return res_df, go_res_analysis, kegg_res_analysis
