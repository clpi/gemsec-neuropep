"""
@author: Savvy Gupta
full output csv of peptides and their AA sequences
additionally with the numerical sequences
with percent similarity calculated with PAM30 matrix
compared to both GrBP5 and M6.
(perfect match score - lowest match score) / PAM30 score
"""

import pandas as pd 
import numpy as np 

def PAM_similarity_score(code, reference):
	df_pam = pd.read_csv('PAM30.csv') 
	if (reference == 'IMVTESSDYSSY'):
		# 8 + 11 + 7 + 7 + 8 + 6 + 6 + 8 + 10 + 6 + 6 + 10
		perfect_match = 93
		# -14 + -13 + -15 + -13 + -17 + -8 + -8 + -15 + -14 + -8 + -8 + -14
		lowest_match = -147
	else:
		# 8 + 11 + 7 + 7 + 6 + 6 + 6 + 6 + 10 + 8 + 8 + 10 
		perfect_match = 93
		# -14 + -13 + -15 + -13 + -13 + -8 + -8 + -13 + -14 + -15 + -15 + -14
		lowest_match = -155

	if (len(code) == len(reference)):
		score = 0;
		i = 0;
		for element in reference:
			if (element.isalpha() and code[i].isalpha()):
				row = df_pam.columns.get_loc(element)
				col = df_pam.columns.get_loc(code[i])
				intersect = df_pam.iloc[row, col]
				score += intersect
				i += 1
		return float((perfect_match - lowest_match) / score)
	else:
		return ""

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
			newNum += ""
	return newNum

def main():
	GRBP5 = 'IMVTESSDYSSY'
	M6 = 'IMVTASSAYDDY'

	df_proteins = pd.read_csv('Organized_NPS.csv')

	# apply function of calculating score
	df_proteins['Similarity to GRBP5'] = df_proteins['Sequences '].apply(PAM_similarity_score, reference = GRBP5)
	df_proteins['Similarity to M6'] = df_proteins['Sequences '].apply(PAM_similarity_score, reference = M6)

	df_proteins['Numerical Codes'] = df_proteins['Sequences '].apply(numerical_conversion)

	df_proteins = df_proteins[df_proteins['Similarity to GRBP5'] != ""]

	df_proteins.to_csv('NPS_Similarity_Scores_PAM30.csv')

if __name__ == "__main__":
	main()

