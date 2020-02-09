import pandas as pd
from scipy.stats import kendalltau
# from scikit-learn import ( ... )

--------------------const

AA = list('LINGVEPHKAYWQMSCTFRD')
AA_COL = "AA_seq"   #--> Column title for AA sequences in .csv file
AA_CHART = pd.read_csv('data/aa_chart.csv')   #--> Path of .csv with AA translations
AA_EIIP = AA_CHART[['AA', 'EIIP']]
AA_FNDOM = AA_CHART[['AA', 'Num', 'Function']]
# get corr between EIIP and Bio Num
EII_FNDOM_CORR, EII_FNDOM_PVAL = kendalltau(AA_CHART[['Num']], AA_CHART[['EIIP']])
GRBP = 'IMVTESSDYSSY'
M6 = 'IMVTASSAYDDY'
PEP_PATH = 'NPS_Numerical_Codes.csv'

#--------------------func

class SequenceSimilarity(object):

    def __init__(self, seq_path, *binders):
        self.binders = [binder for binder in binders]
        try:
            self.seq = pd.Dataframe(seq_path, index=AA_COL)
        except:
            print("No .csv at that location.")


def df_filter_subseq(data: pd.DataFrame, sub_seq: str, ind: int = None):
    if not {*sub_seq}.issubset({*AA}):
        raise Exception('Invalid subsequence')
    if ind is None:
        return data.filter()

    def domain_filter(dom, seq, i):
        return seq[i:len(dom)] == dom

    data_filter = data[AA_COL].apply(lambda s: s[ind:len(dom) == dom])
    data_with_seq = data[data[AA_COL].filter(data_filter)]
    return data[data[AA_COL].filter(data_filter)]

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
    pass

if __name__ == '__main__':
    main()

