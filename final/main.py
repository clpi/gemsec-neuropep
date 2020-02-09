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
        if not {*sub_seq}.issubset({*AA}):
            raise Exception('Invalid subsequence')
        if ind is None:
            return self.data.filter() #@TODO Implement the filter here

        data_filter = data[AA_COL].apply(lambda s: s[ind:len(dom) == dom])
        data_with_seq = data[data[AA_COL].filter(data_filter)]
        return data[data[AA_COL].filter(data_filter)]

    def get_PAM30_similarity(self):
        raise NotImplementedError

    def get_BLOSUM_similarity(self):
        raise NotImplementedError

    def get_RRM_SN_ratio(self):
        get_eiip_seq = lambda pep: list(map(lambda aa: AA_EIIP[aa], pep))
        get_dft_from_eiip = lambda eiip: np.fft.rfft(eiip)[1:]
        get_cross_spectrum = lambda p1, p2: [x1*x2 for x1, x2 in zip(p1, p2)]

def get_similarity_from_subseq(subseq: str):

"""
@author: Savvy Gupta
"""

domain = input("Enter a domain to select for: ")
index = int(input("Enter index at which domain starts: "))

def domain_selector(sequence):

	subSequence = sequence[index : len(domain)]

	if (subSequence == domain):
		return sequence
	else:
		return "-1"

'''
def main():

	df_filter = pd.read_csv('NPS_Numerical_Codes.csv', index_col=[0])
	df_filter['Sequences '] = df_filter['Sequences '].apply(domain_selector)
	nonDomain = df_filter[df_filter['Sequences '] == '-1'].index
	df_filter.drop(nonDomain, inplace = True)
	df_filter.to_csv('NPS_Filtered_Domain_' + str(domain) + '_Start_' + str(index) + '.csv', index=False)
'''

#-------------------main
def main():
    seq_sim = SequenceSimilarity(PEP_PATH, GRBP5, M6)

if __name__ == '__main__':
    main()

