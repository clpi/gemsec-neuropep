'''
Title: main.py placeholder / skeleton program
Author: Chris Pecunies
Description: This is the program that will use the statistical PAM30 and cross-entropy similarity calculations (Savvy) and RRM signal-noise and correlation coefficient calculations (Aaron), as well as unsupervized statistical clustering to determine important sequence domains relative to both biological function mapping (provided in documentation) and EIIP mapping (also in docs). 
'''

import os
import pandas as pd
from scipy.stats import kendalltau
# from scikit-learn import ( ... )

#--------------------const
AA_CHART = pd.read_csv('data/aa_chart.csv')
AA_EIIP = AA_CHART[['AA', 'EIIP']]
AA_FNDOM = AA_CHART[['AA', 'Num', 'Function']]
# get corr between EIIP and Bio Num
EII_FNDOM_CORR, EII_FNDOM_PVAL = kendalltau(AA_CHART[['Num']], AA_CHART[['EIIP']])
GRBP = 'IMVTESSDYSSY'
M6 = 'IMVTASSAYDDY

#--------------------func
# @TODO: Implement cross correlation
# @TODO: Implement RRM (correct) w/ biolog. fn seq
# @TODO: Investigate kendall tau as a measure of seq sim.
def get_

def get_input():
    pass

def process_input(csv_in, *seqs):
    pass

#-------------------main
def main():
    pass

if __name__ == '__main__':
    main()

