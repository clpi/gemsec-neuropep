# -*- coding: utf-8 -*-
"""
@author: Aaron Tsang
"""
import pywt 
import cmath as math
import numpy as np
import matplotlib.pyplot as plt
import statistics as stat

AA = ['L','I','N','G','V','E','P','H','K','A','Y','W','Q','M','S','C','T','F','R','D']
EIIP = [0,0,0.0036,0.005,0.0057,0.0058,0.0198,0.0242,0.0371,0.0373,0.0516,0.0548,0.0761,0.0823,0.0829,0.0829,0.0941,0.0946,0.0959,0.1263]

xaxis = []
cross = []

# amino sequence -> EIIP -> fft
seq = ['PALPEDGGSGAFPPGHFKDPKRLYCKNGGFFLRIHPDGRVDGVREKSDPHIKLQLQAEERGWSIKGVCANRYLAMKEDGRLLASKCVTDECFFFERLESNNYNTYRSRKYSSWYVALKRTGQYKLGPKTGPGQKAILFLPMSAKS']
      # 'FNLPLGNYKKPKLLYCSNGGYFLRILPDGTVDGTKDRSDQHIQLQLCAESIGEVYIKSTETGQFLAMDTDGLLYGSQTPNEECLFLERLEENHYNTYISKKHAEKHWFVGLKKNGRSKLGPRTHFGQKAILFLPLPVSSD']
newpep = []
result = {}

xaxis = []
for i in range(0, int((len(seq[0])/2) + 1)):
    xaxis.append((i/len(seq[0])/2)/2)
    #xaxis.append(i)
    
newpep = []

# from amino acids to EIIP value
for nums in seq[0]:
    newpep.append(EIIP[AA.index(nums)])
    result[seq[0]] = pywt.dwt(newpep, 'db1')


#
#pep2 = []
#holder = 0
#for nums in range(0, len(seq[0])):
#   if(nums % 2 == 0):
#        pep2.append(abs(result[seq[0]][nums] * math.exp((1j ** -1) * ((2 * math.pi)/len(seq[0])) * holder * holder).real))
#        holder = holder + 1
#        



final = []
for amino in result[seq[0]]:
    final.append(amino/max(result[seq[0]]))
result[seq[0]] = final
    

#
#
plt.plot(xaxis, result[seq[0]])
#plt.title("with norm = ortho")
#plt.savefig("test.png")



    


    

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
