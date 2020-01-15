# -*- coding: utf-8 -*-
"""
@author: Aaron Tsang
"""
import pywt 
import numpy as np
import matplotlib.pyplot as plt
import statistics as stat
import csv
import math

AA = ['L','I','N','G','V','E','P','H','K','A','Y','W','Q','M','S','C','T','F','R','D']
EIIP = [0,0,0.0036,0.005,0.0057,0.0058,0.0198,0.0242,0.0371,0.0373,0.0516,0.0548,0.0761,0.0823,0.0829,0.0829,0.0941,0.0946,0.0959,0.1263]

def pep_to_dft(peptide):
    eiip_value = []
    for aa in peptide:
        eiip_value.append(EIIP[AA.index(aa)])
    dft_value = np.fft.rfft(eiip_value)[1:]
    final_value = []
    for aa in dft_value:
        final_value.append(aa / max(dft_value))
    return final_value

        

GRBP = 'IMVTESSDYSSY'
subset = []
with open ('9_17mers.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter = ",")
    
    for row in csv_reader:
        peptide = "".join(row)
        rules = ['O' not in peptide, 'J' not in peptide,'/' not in peptide, len(peptide) == len(GRBP)]
        if all(rules):            
            subset.append(peptide)

dict = []
for peptides in range(0, len(subset)):
    dict.append(pep_to_dft(subset[peptides]))

grbp_dft = pep_to_dft(GRBP)


dict2 = []
for peptides in range(0, len(dict)):
    for aa in range(0, len(dict[peptides])):
        temp = dict[peptides][aa] * np.conj(grbp_dft[aa])
        cross_spec = []
        cross_spec.append((float)(math.sqrt( (temp.real)** 2 + (temp.imag)** 2)))
    dict2.append(cross_spec[0])
    
dict3 = []

# Calculate for S/N ratio as specified in the book
for peptides in range(0, len(dict2)):
    if(dict2[peptides]/(sum(dict2) / len(dict2)) > 20.0):
        toWrite = []
        toWrite.extend([dict2[peptides], subset[peptides]])
        dict3.append(toWrite)
        
# Saving CSV file
with open('set1.csv', 'w', newline = '') as write:
    csv_writer = csv.writer(write)
    
    csv_writer.writerow(['Peptide', 'Similarity score'])
    
    for row in dict3:
        csv_writer.writerow([row[1], row[0]])
    





    
  
