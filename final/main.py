import pandas as pd
import numpy as np
from scipy.stats import kendalltau
import matplotlib.pyplot as plt
# from scikit-learn import ( ... )

# --------------------const

# get corr between EIIP and Bio Num
EII_FNDOM_CORR, EII_FNDOM_PVAL = kendalltau(AA_CHART[['Num']], AA_CHART[['EIIP']])
GRBP5 = 'IMVTESSDYSSY'
M6 = 'IMVTASSAYDDY'
PEP_PATH = 'NPS_Numerical_Codes.csv'

#--------------------func

class SequenceSimilarity(object):
    '''
    Class that takes in a path to a list of amino acid sequences as well
    as any number of peptide sequences explicitly that are known to have
    a certain set of properties. Generates metrics for similarity for each
    peptide in path and returns domains AA sequence with high similarity
    '''

    def __init__(self, seq_path, *binders):

        self.AA_COL = "AA_seq"   #--> Column title for AA sequences in .csv file
        self.AA_CHART = pd.read_csv('data/aa_chart.csv')   #--> Path of .csv with AA translations
        self.AA_EIIP = AA_CHART[['AA', 'EIIP']]
        self.AA_FNDOM = AA_CHART[['AA', 'Num', 'Function']]
        self.AA = list('LINGVEPHKAYWQMSCTFRD')
        self.binders = [binder for binder in binders]
        try:
            self.data = pd.Dataframe(seq_path, index=AA_COL)
        except:
            print("No .csv at that location.")

    def df_filter_subseq(self, data: pd.DataFrame, sub_seq: str, ind: int = None):
        '''
        Takes in a subsequence of equal or lesser length to
        peptides in class peptide dataframe and returns a dataframe
        containing only those peptides containing the sequence
        '''
        if not {*sub_seq}.issubset({*self.AA}):
            raise Exception('Invalid subsequence')
        if ind is None:
            return self.data.filter() #@TODO Implement the filter here

        data_filter = self.data[self.AA_COL].apply(
            lambda s: s[ind:len(sub_seq) == sub_seq])
        return data[data[self.AA_COL].filter(data_filter)]

    def get_binder_subseq(self):
        '''
        Generates all possible subsequences for binders
        provided in class constructor
        '''
        def gen_all_subseq(seq, sub_seq, i):
            if i == len(seq):
                if len(sub_seq) != 0:
                    yield(sub_seq)
                else:
                    gen_all_subseq(seq, sub_seq, sub_seq)
                gen_all_subseq(seq, sub_seq+[seq[i]], i+1)

        sseq = dict.fromkeys(self.binders)
        for binder in self.binders:
            sseq[binder] = [sseq for sseq in list(
                gen_all_subseq(binder, '', 0))]
        return sseq

    def get_PAM30_similarity(self):
        raise NotImplementedError

    def get_BLOSUM_similarity(self):
        raise NotImplementedError
    '''
    def get_RRM_SN_ratio(self):
        get_eiip_seq = lambda pep: list(map(lambda aa: self.AA_EIIP[aa], pep))
        get_dft_from_eiip = lambda eiip: np.fft.rfft(eiip)[1:]
    '''
def main():
    seq_sim = SequenceSimilarity(PEP_PATH, GRBP5, M6)

if __name__ == '__main__':
    main()

