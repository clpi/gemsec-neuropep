import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from typing import Set, Tuple, Dict, List
from scipy import interpolate, stats, fftpack, signal
#import sci-kit learn
import textdistance as td
#import biopython as bio


def main():
# Consider using Organized NPS from Savvy's folder?
    DATA_PATHS = {
        "BLOSUM":"src_data/BLOSUM.csv",
        "PAM30":"src_data/pam30.csv",
    }
    SEQS = [
        ('GRBP5','IMVTESSDYSSY'),
        ('M6','IMVTASSAYDDY')
    ]
    AA_COL = 'Sequences'
    PEP_PATH = 'src_data/Sequence_data.csv' 
    dat1 = {
        'grbp5_sim' : SequenceSimilarity(SEQS[0], DATA_PATHS, PEP_PATH, AA_COL, True, False),
        'm6_sim' : SequenceSimilarity(SEQS[1], DATA_PATHS, PEP_PATH, AA_COL, True, False),
        'grbp5_sim_match' : SequenceSimilarity(SEQS[0], DATA_PATHS, PEP_PATH, AA_COL, True, True),
        'm6_sim_match' : SequenceSimilarity(SEQS[1], DATA_PATHS, PEP_PATH, AA_COL, True, True),
    }
    dat2 = {
        'both_sep_sim' : dat1['grbp5_sim'].merge_data(other=dat1['m6_sim'], sep_cols=True),
        'both_avg_sim' : dat1['grbp5_sim'].merge_data(other=dat1['m6_sim'], sep_cols=False),
        'both_sep_sim_match' : dat1['grbp5_sim_match'].merge_data(other=dat1['m6_sim_match'], sep_cols=True),
        'both_avg_sim_match' : dat1['grbp5_sim_match'].merge_data(other=dat1['m6_sim_match'], sep_cols=False),
    }

    d1 = list(dat1.values())
    d1d = [d.pep_data for d in d1]
    d2d = list(dat2.values())
    print([len(d) for d in d1d] + [len(d) for d in d2d])

    for name in dat1.keys():
        path = "../out/v2/{}_samelen_allmetrics_v2.csv".format(name)
        dat1[name].pep_data.to_csv(path)
    for i, name in enumerate(dat2.keys()):
        if i == 1 or 1 == 3:
            cols="sep_cols"
            path = "../out/v2/{}_samelen_allmetrics_v2_{}.csv".format(name, cols)
            dat2[name].to_csv(path)
        if i == 0 or i == 2:
            cols="avg cols"
            path2 = "../out/v2/{}_samelen_allmetrics_combined_m6_grbp5_{}.csv".format(name,cols)
            dat2[name].to_csv(path2)

    new1 = dat2['both_avg_sim_match'].groupby('FNS_Seq').mean()
    new2 = dat2['both_avg_sim_match'].groupby('Sequences')

    dat2['both_avg_sim_match'].describe()

    import seaborn
    import scipy.signal

    pd.options.display.max_rows = 999
    pd.set_option('mode.sim_interactive', True)
    pd.set_option('expand_frame_repr', True)
    pd.set_option('large_repr', 'truncate')

    for name in dat1.keys():
        d = dat1[name].pep_data #.plot()
        #print(d.head())
        #df = d.groupby('NUM_Seq').plot()
        #ds = d.sort_values(by="weighted_matches").plot()
        #df.head(10)
        #dat1[name].pep_data.hist()
    #dat2['both_avg_sim_match'][[*distances.keys()] + ["NUM_Seq"]].groupby("NUM_Seq").mean().plot()
    d1 = list(dat1.keys())
    d2 = list(dat2.keys())

    dat2[d2[0]][['sseq_matches_GRBP5', 'sseq_matches_M6']]

    def get_RRM():
        data = similarity.pep_data
        matrix = similarity.data['PA]
        print(data.head(), len(data))

    get_RRM()

    SEQS = [
        ('GRBP5','IMVTESSDYSSY'), #->, -----22---220
        ('M6','IMVTASSAYDDY')    #->  -----22---550
    ]
    s1 = np.fft.rfft(list(dat1['grbp5_sim'].pep_data['EIIP_Seq'].iloc[1]))
    s1 = np.fft.rfft(list(dat1['grbp5_sim'].pep_data['EIIP_Seq'].iloc[2]))
    s1 = np.fft.rfft(list(dat1['grbp5_sim'].pep_data['EIIP_Seq'].iloc[3]))
    plt.plot(np.real(s1))
    plt.show()

    sorted_grbp5_sim_match = data['grbp5_sim_match'].pep_data.sort_values(by="weighted_matches")
                                 
if __name__ == '__main__':
    main()