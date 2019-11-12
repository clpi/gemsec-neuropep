# -*- coding: utf-8 -*-
"""
@author: Aaron Tsang
"""
import math
import numpy as np
import matplotlib.pyplot as plt
import statistics as stat

AA = ['L','I','N','G','V','E','P','H','K','A','Y','W','Q','M','S','C','T','F','R','D']
EIIP = [0,0,0.0036,0.005,0.0057,0.0058,0.0198,0.0242,0.0371,0.0373,0.0516,0.0548,0.0761,0.0823,0.0829,0.0829,0.0941,0.0946,0.0959,0.1263]

xaxis = []
cross = []

# amino sequence -> EIIP -> fft
seq = ['PALPEDGGSGAFPPGHFKDPKRLYCKNGGFFLRIHPDGRVDGVREKSDPHIKLQLQAEERGWSIKGVCANRYLAMKEDGRLLASKCVTDECFFFERLESNNYNTYRSRKYSSWYVALKRTGQYKLGPKTGPGQKAILFLPMSAKS',
       'FNLPLGNYKKPKLLYCSNGGYFLRILPDGTVDGTKDRSDQHIQLQLCAESIGEVYIKSTETGQFLAMDTDGLLYGSQTPNEECLFLERLEENHYNTYISKKHAEKHWFVGLKKNGRSKLGPRTHFGQKAILFLPLPVSSD']
newpep = []
result = {}

xaxis = []
for i in range(0, len(seq[0])):
    xaxis.append((i/145)/2)
    #xaxis.append(i)
    
newpep = []

# take the fourier transform, then the real part of the EIIP
for nums in seq[0]:
    newpep.append(EIIP[AA.index(nums)])
    result[seq[0]] = newpep
    
newpep = []
for nums in result[seq[0]]:
    newpep.append(nums * math.exp())
    


    

## find max
#max = 0
#for nums in result[seq[0]]:
#    if max < nums:
#        max = nums
#
## take the percentage out of 100% 
#newpep = []
#for nums in result[seq[0]]:
#    newpep.append(nums/ max)
#result[seq[0]] = newpep
#
## plot
#plt.plot(xaxis, result[seq[0]])
#
#
#
#
##
##newpep = []
##for nums in seq[1]:
##    newpep.append(EIIP[AA.index(nums)])
##for i in range (len(seq[1]), len(seq[0])):
##    newpep.append(0)
###result[seq[1]] = np.fft.fft(newpep)
##result[seq[1]] = newpep
##
##plt.plot(xaxis, result[seq[1]])
##    
##
#
#
#        
#        
