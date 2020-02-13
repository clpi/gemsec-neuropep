# Resonant Recognition Model implementation
# Author: Chris Pecunies

import numpy as np
import matplotlib.pyplot as plt

# set up eiip dict
AA = list('LINGVEPHKAYWQMSCTFRD')
EIIP = [0,0,0.0036,0.005,0.0057,0.0058,0.0198,0.0242,0.0371,0.0373,0.0516,0.0548,0.0761,0.0823,0.0829,0.0829,0.0941,0.0946,0.0959,0.1263]
AA_EIIP = dict(zip(AA, EIIP))

# helper functions
get_eiip_seq = lambda pep: list(map(lambda aa: AA_EIIP[aa], pep))
get_dft_from_eiip = lambda eiip: np.fft.rfft(eiip)[1:]
get_cross_spectrum = lambda p1, p2: [x1*x2 for x1, x2 in zip(p1, p2)]

if __name__ == '__main__':
  seq1 = list('PALPEDGGbSGAFPPGHFKDPKRLYCKNGGFFLRIHPDGRVDGVREKSDPHIKLQLQAEERGWSIKGVCANRYLAMKEDGRLLASKCVTDECFFFERLESNNYNTYRSRKYSSWYVALKRTGQYKLGPKTGPGQKAILFLPMSAKS')
  seq2 = list('FNLPLGNYKKPKLLYCSNGGYFLRILPDGTVDGTKDRSDQHIQLQLCAESIGEVYIKSTETGQFLAMDTDGLLYGSQTPNEECLFLERLEENHYNTYISKKHAEKHWFVGLKKNGRSKLGPRTHFGQKAILFLPLPVSSD')
  seq3 = list('VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR') #alpha hemoglobin
  seq4 = list('VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH') #beta hemoglobin
  eiip1 = get_eiip_seq(seq1)
  eiip2 = get_eiip_seq(seq2)
  dft1 = get_dft_from_eiip(eiip1)
  dft2 = get_dft_from_eiip(eiip2)
  cross_spectrum = get_cross_spectrum(dft1, dft2)

  fig = plt.figure(figsize=(10, 8))
  ax = fig.add_subplot(111)
  ax.plot(cross_spectrum)
  plt.show()
  fig.savefig('graph.png')
