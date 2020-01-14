# -*- coding: utf-8 -*-
"""
@author: Aaron Tsang
"""
# import pywt 
import cmath as cmath
import math as math
import numpy as np
import matplotlib.pyplot as plt
import statistics as stat

AA = ['L','I','N','G','V','E','P','H','K','A','Y','W','Q','M','S','C','T','F','R','D']
EIIP = [0,0,0.0036,0.005,0.0057,0.0058,0.0198,0.0242,0.0371,0.0373,0.0516,0.0548,0.0761,0.0823,0.0829,0.0829,0.0941,0.0946,0.0959,0.1263]

xaxis = []
cross = []

# FGF amino acid sequence
seq = ['PALPEDGGSGAFPPGHFKDPKRLYCKNGGFFLRIHPDGRVDGVREKSDPHIKLQLQAEERGWSIKGVCANRYLAMKEDGRLLASKCVTDECFFFERLESNNYNTYRSRKYSSWYVALKRTGQYKLGPKTGPGQKAILFLPMSAKS']
      # 'FNLPLGNYKKPKLLYCSNGGYFLRILPDGTVDGTKDRSDQHIQLQLCAESIGEVYIKSTETGQFLAMDTDGLLYGSQTPNEECLFLERLEENHYNTYISKKHAEKHWFVGLKKNGRSKLGPRTHFGQKAILFLPLPVSSD']
result = {}

    
newpep = []

# from amino acids to EIIP value
for nums in range (0, int((len(seq[0]))/2)):
    newpep.append(EIIP[AA.index(seq[0][nums])])
result[seq[0]] = (np.fft.rfft(newpep)[1:])

xaxis = range(0, len(newpep))

# manually computing the DFT 
pep2 = []
holder = 1
for nums in range (0, len(newpep)):
        pep2.append(abs(newpep[nums] * cmath.exp((1j ** -1) * ((2 * math.pi)/len(seq[0])) * nums * holder).real ))
        holder = holder + 1
        value = 0
newpep = []        
for nums in range(0, len(pep2)):
    newpep.append(pep2[nums] / max(pep2))
pep2 = newpep

plt.plot(xaxis, pep2)
plt.ylim(0,1)
plt.savefig("final.png")
