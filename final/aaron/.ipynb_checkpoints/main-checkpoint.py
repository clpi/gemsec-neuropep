# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 12:54:02 2019

@author: Admin
"""
import csv
import matplotlib.pyplot as plt
import statistics as stat
import numpy as np


AA = ['L','I','N','G','V','E','P','H','K','A','Y','W','Q','M','S','C','T','F','R','D']
EIIP = [0,0,0.0036,0.005,0.0057,0.0058,0.0198,0.0242,0.0371,0.0373,0.0516,0.0548,0.0761,0.0823,0.0829,0.0829,0.0941,0.0946,0.0959,0.1263]



seq = []
test = 'IMVTESSDYSSY'
testnum = []

    



with open('9_17mers.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter = ",")
    
    #change the spreadsheet to csv form with each element being strings
    for row in readCSV:
        string = "".join(row)
        rules = ['O' not in string, 'J' not in string,'/' not in string]
        if all(rules):            
            seq.append(string)
                     
newpep = []
less = []
more = []

for letters in test:
    newpep = []
    newpep.append(EIIP[AA.index(letters)])
    testnum.append(np.fft.fft(newpep))

# from EIIP -> fft 
for seqs in seq:
    newpep = []
    for letters in seqs:
        newpep.append(EIIP[AA.index(letters)]) 
    if len(newpep) <= len(test):
        less.append(np.fft.fft(newpep))
        
    else:
        more.append(np.fft.fft(newpep))

crosspec = []
temp = []
#        
## for sequences less than GRBP-5
#for seqs in less:
#    for i in range(0, len(testnum)):
#        newpep = []
#        if(i == len(seqs)):
#            print('AMINO:' + '0')
#            print('GRBP:' + testnum[i])
#        print('AMINO: ' + seqs[i])
#        print('GRBP:' + testnum[i])
#        
#        
#

        
        
        