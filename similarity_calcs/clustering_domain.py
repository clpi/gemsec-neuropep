"""
@author: Savvy Gupta
"""

import pandas as pd 
import numpy as np 

domain = input("Enter a domain to select for: ")
index = int(input("Enter index at which domain starts: "))

def domain_selector(sequence):

	subSequence = sequence[index : len(domain)]

	if (subSequence == domain):
		return sequence
	else:
		return "-1"


def main():

	df_filter = pd.read_csv('NPS_Numerical_Codes.csv', index_col=[0])
	df_filter['Sequences '] = df_filter['Sequences '].apply(domain_selector)
	nonDomain = df_filter[df_filter['Sequences '] == '-1'].index
	df_filter.drop(nonDomain, inplace = True)
	df_filter.to_csv('NPS_Filtered_Domain_' + str(domain) + '_Start_' + str(index) + '.csv', index=False)

if __name__ == "__main__":
	main()