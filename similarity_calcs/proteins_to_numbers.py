"""
@author: Savvy Gupta
"""

import pandas as pd 
import numpy as np

def numerical_conversion(sequence):
	clusterToNum = {
		'C': 7,
		'R': 6,
		'H': 6,
		'K': 6,
		'D': 5,
		'E': 5,
		'G': 4,
		'P': 3,
		'S': 2,
		'T': 2,
		'N': 2,
		'Q': 2,
		'A': 1,
		'V': 1,
		'I': 1,
		'L': 1,
		'M': 1,
		'F': 0,
		'Y': 0,
		'W': 0,
	}
	newNum = ""
	for letter in sequence:
		if letter.upper() in clusterToNum.keys():
			newNum += str(clusterToNum[letter.upper()])
		else:
			# still missing some letters: X, J, O, B
			# some had non-alphabetical entries: `, ., /, ", 1, which I removed
			newNum += ""
	return newNum

def main():
	
	df_proteins = pd.read_csv('Organized_NPS.csv')
	df_proteins['Sequences '] = df_proteins['Sequences '].apply(numerical_conversion)
	df_proteins.to_csv('NPS_Numerical_Codes.csv')

if __name__ == "__main__":
	main()



