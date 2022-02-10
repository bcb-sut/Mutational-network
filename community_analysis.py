import plotly_template
from vis_utils import VisUtils


def com_cover_patient(com, patient_genes):
    if len(set(com).intersection(set(patient_genes))) == len(com):
        return True
    return False


class CommunityAnalysis:
    def __init__(self, coms, save_path, data_loader):
        self.coms = coms
        self.patient_gene = data_loader.load_patient_gene()
        self.patients = list(self.patient_gene.keys())
        self.vis = VisUtils(save_path)

    def count_covering_coms(self, patient):
        cnt = 0
        for com in self.coms:
            if com_cover_patient(com, self.patient_gene[patient]):
                cnt += 1
        return cnt

    def count_covering_patients(self, com):
        cnt = 0
        for patient in self.patients:
            if com_cover_patient(com, self.patient_gene[patient]):
                cnt += 1
        return cnt

    def sample_coverages(self):
        print('count covering patients...')
        com_patient_cnt = [self.count_covering_patients(com) for com in self.coms]
        print('done')
        com_patient_coverage = [cnt / len(self.patients) for cnt in com_patient_cnt]
        cumulative_coverage = [com_patient_coverage[i - 1] + com_patient_coverage[i] for i in range(1, len(self.coms))]
        covered_com_cnt = [self.count_covering_coms(patient) for patient in self.patients]

        all_coms_patients = covered_com_cnt.count(len(self.coms))
        print(f'{all_coms_patients} number of patients are covered by all coms.')

        com_sizes = [len(com) for com in self.coms]
        com_ids = list(range(len(self.coms)))
        rows, cols = 2, 2

        plots = [plotly_template.scatter(com_sizes, com_patient_coverage),
                 plotly_template.scatter(com_ids, cumulative_coverage),
                 plotly_template.histogram(covered_com_cnt),
                 plotly_template.scatter(com_ids, covered_com_cnt)]
        plot_types = ['Scatter', 'Scatter', 'Histogram', 'Scatter']
        xlabels = ['Community size', 'Community size', 'Communities count', 'Patient id']
        ylabels = ['Patient coverage', 'Cumulative patient coverage',
                   'Covered patients count', 'Number of communities']
        titles = ['a', 'b', 'c', 'd']
        fig = plotly_template.make_multiple_plots(plots, rows, cols, plot_types, titles, xlabels, ylabels)
        self.vis.save_and_show_plotly(fig, 'coms-coverages')
        return com_patient_coverage
