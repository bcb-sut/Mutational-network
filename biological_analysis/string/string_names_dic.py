class StringNames:
    def __init__(self, data_folder):
        self.protein_info_path = data_folder + '/stringdb/9606.protein.info.v11.5.txt'
        self.gene_dic = self.generate_gene_dic()

    def generate_gene_dic(self):
        name_dic = {}
        print('generating id-name dic...')
        f = open(self.protein_info_path)
        line = f.readline()
        while line:
            line = f.readline()
            try:
                id_, name = line.split()[:2]
                name_dic[id_] = name
                name_dic[name] = id_
            except:
                print(f'cant split line: --{line}--')
        print(f'done. dic size: {len(name_dic)}')
        return name_dic

    def translate_gene(self, gene):
        try:
            return self.gene_dic[gene]
        except:
            print(f'{gene} not found in gene_dic')
            return gene
