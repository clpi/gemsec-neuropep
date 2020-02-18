'''
Class to calculate several simialrity metrics for an input list of peptides given a sequence string to compare it to.
USAGE: 
1. Create SequenceSimilarityObject({sequence string}, {dictionary of data paths} (to be deprecated soon -- use
   internal class data), {path to peptides csv}, {column title of sequence values of aforementioned peptides csv})
   for any number of sequences you want to be compared to the list of peptides
2. For each object, call object.generate_similarity() to fill out the Dataframe with similarity metrics
3. To whittle down this similarity matrix to list only those peptides with pattern matching of a minimum length
   at a matching index in the rerence peptide (henceforth the binder), call object.get_df_with_binder_subseqs(min_length={#})
@ Author: Chris Pecunies, with help from Savvy Gupta and Aaron Tsang
@ Date: February 12, 2020
'''
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from typing import Set, Tuple, Dict, List
from scipy import stats, signal, fft
#import sci-kit learn
import textdistance as td
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist

class SequenceSimilarity:
    '''
    Class that takes in a path to a list of amino acid sequences as well
    as any number of peptide sequences explicitly that are known to have
    a certain set of properties. Generates metrics for similarity for each
    peptide in path and returns domains AA sequence with high similarity
    '''
    AA = list('FYWAVILMSTNQPGDEKHRC')
    
    def __init__(self, binder: Tuple[str, str],
                 data: SeqData,
                 p_path: str,       
                 aa_col: str,
                 dists: List = [],   #set dists = None for no dists
                 only_match: bool = False,
                 min_length: int = 0): 
        
        # ---- setting data ---------
        self.d = {d: vars(data)[d] for d in list(vars(data).keys())}
        self.conv = list(self.d['conv'].keys())
        self.aa_map = {conv:dict(zip(self.AA, self.d['conv'][conv])) for conv in self.conv}
        self.conv_pep = lambda t, p: [self.aa_map[t][AA] if AA in self.AA else 0 for AA in p]
        self.sim_sum = lambda p1, p2, m: sum([self.d['matr'][m][a1][a2] for a1, a2 in zip(p1, p2)])
        
        self.seq = aa_col
        self.sim_cols = ['PAM30', 'BLOSUM45', 'RRM_SN', 'RRM_Corr', 'weighted_matches']
        if dists is not None:
            self.sim_cols += self.d['dist'].keys() if not dists else dists
        self.conv_cols = [conv+'_Seq' for conv in self.conv]
        self.cols = [aa_col] + self.sim_cols + self.conv_cols
        self.bname, self.b = binder
        self.bsseq = [(self.b[i:j], i) for i in range(len(self.b)) for j in range(i+1, len(self.b)+1)]
    
        # ---- helper lambda functions
        
        self.p_og = pd.read_csv(p_path)
        self.p_og.columns = [aa_col]
        self.p_og[aa_col].drop_duplicates()
        self.p_og = self.p_og[~self.p_og[aa_col].str.contains("O")]
        self.p_sl = self.p_og[self.p_og[aa_col].str.len() == len(self.b)]
        self.p_ls = self.p_sl[aa_col].tolist()
        self.p_og = self.p_og[aa_col].tolist()
        if len(self.p_ls) == 0:  
            raise Exception("No peptides of same length as binder found")
        self.p = pd.DataFrame(columns=self.cols)
        self.p[self.seq] = self.p_ls
        # @TODO P STILL HAS DUPLICATES ????
        self.update_similarities(dists)
        
        self.p_match = self.filter_by_bsseq(min_length)
        self.p_match_ls = self.p_match[self.seq].tolist()
        if only_match: 
            self.p, self.p_ls = self.p_match, self.p_match_ls
        
    #----------------SET UP FUNCTIONS (void)-------------------------------#     

    def _update_RRM_similarity(self, match_len = True) -> None:
        """
        Uses the Resonant Recognition Model as described by Irena Cosic to 
        """
        rrm = pd.DataFrame(index=self.p.index)
        eiip_seq = self.p['EIIP_Seq'].tolist()
        rrm = rrm.assign(
            seq = self.p_ls,
            eiip = eiip_seq,
            dft = [np.fft.rfft(eiip)/len(eiip) for eiip in eiip_seq],
            dft2 = [2*np.abs(np.fft.rfft(eiip)/len(eiip)) for eiip in eiip_seq],
            freq = [np.abs(np.fft.fftfreq(len(eiip))[0:int(len(eiip)/2+1)]) for eiip in eiip_seq],
            power = [np.abs(np.fft.rfft(eiip)/len(eiip))**2 for eiip in eiip_seq],
        )
        bind_eiip = self.conv_pep('EIIP', self.b)
        self.bind_eiip_dft = 2*np.abs(np.fft.rfft(bind_eiip)/len(bind_eiip))
        self.bind_freq = np.abs(np.fft.fftfreq(len(bind_eiip))[0:int(len(bind_eiip)/2+1)])
        self.bind_peak = signal.find_peaks(self.bind_eiip_dft)
        peaks = [signal.find_peaks(d) for d in rrm.dft2]
        rrm['peak_loc'] = [p[0] for p in peaks]
        rrm['peak_val'] = [p[1] for p in peaks]
        rrm['peak_dist'] = [np.abs(p[0]-self.bind_peak[0]) for p in peaks]
        rrm['correlate'] = [signal.correlate(d, self.bind_eiip_dft) for d in rrm.dft2]
        rrm['convolve'] = [signal.convolve(d, self.bind_eiip_dft) for d in rrm.dft2]
        self.rrm = rrm
        
        def get_max_seqs(num = 10, typ='correlate'):
            maxes = [signal.find_peaks(rrm.iloc[n].correlate) for n in range(len(rrm))]
            maxheight = [rrm[typ].iloc[i][m[0][0]] for i, m in enumerate(maxes)]
            top_idx = np.argsort(maxheight)[-num:]
            top_values = [maxheight[i] for i in top_idx]
            top_seq = [top_seq.append(rex['seq'].iloc[i]) for i in top_idx]
            top_data = self.p[self.p[self.seq].isin(top_seq)]
            return top_data
            
        def plot(col = 'column'): plt.plot(rrm.[col])
    
    
    # NOTE! Adds columns "Matching_sseqs" and "Num_matching" to output
    # Might be too unwieldy / unhelpful for output similarity data
    # if so, just comment out _update_matching_sseqs()
    def _update_matching_sseqs(self, w1: float = 1, w2: float = 1) -> None:
        """
        Returns a number as a new column representing the number of "matches" a peptide
        has for all possible subsequences for the binder inputted at a given index. For
        weighting=1, all matches are treated equally ('Y' at position 3 is treated equal
        to IMV at position 0) but lowering weighting lowers smaller-length matches
        """
        # @TODO Remove "duplicates" which occur at different matching indexes of binder
        # but are part of a larger pattern already recorded at an earlier index
        self.p.sseq_matches, self.p.weighted_matches = None, None
        self.bsseq.sort(key = lambda ss: len(ss[0]), reverse=True)
        score: float = lambda s: (w1 * 1) + (len(s)**w2)
        matches: List = list(); nmatches = list()
        for i, seq in enumerate(self.p_ls):
            matches.append(list()); nmatches.append(int())
            trigger = False
            for j in range(len(seq)):
                if trigger: break
                trigger = False
                for (sseq, bin_i) in self.bsseq:
                    in_seq = seq[j:len(sseq)+j]
                    if (bin_i == j) and (in_seq == sseq):
                        if len(matches[i]) > 0 and matches[i][-1][0].find(sseq)>=0:
                            trigger = True
                            break
                        matches[i].append((sseq, bin_i))
                        nmatches[i] += score(sseq)
                        break
#         matches = [pairwise2.align.localxx(s, self.b) for s in self.p_ls]
        self.p.sseq_matches, self.p.weighted_matches = matches, nmatches
                
    def update_similarities(self, metrics: List = []) -> None:
        '''
        Updates the similarity values whenever called (for now should be only once right
        after creating the object, ecept possibly if the Binding peptide is updated
        (should be handled automatically)
        '''
        matrices = list(self.d['matr'].keys())
        cdata = [[self.conv_pep(c, p) for c in self.conv] for p in self.p_ls]
        mdata = [[self.sim_sum(p, self.b, m) for m in matrices] for p in self.p_ls]
        self.p[self.conv_cols] = cdata
        self.p[matrices] = mdata
        self._update_matching_sseqs()
        self._update_RRM_similarity()
        if metrics is not None:
            dists = metrics if not len(metrics) == 0 else list(self.d['dist'].keys())
            self.p[dists] = [[self.d['dist'][d](p, self.b) for d in dists] for p in self.p_ls]
        # OPTIONAL
        # self._unpack_num_encoding()
        
        
        
    #----------MAIN CLASS FUNCTIONS (returns data) ------------------------#
    
    def filter_by_sseq(self, sseq: str, ind: int):
        return self.p[self.p[self.seq].str.find(sseq) == ind]
        
    def filter_by_bsseq(self, min_len: int = 0) -> pd.DataFrame:
        bsseq_dfs = [self.filter_by_sseq(ss,i) for (ss,i) in self.bsseq if len(ss)>=min_len]
        return pd.concat(bsseq_dfs)
    
    def get_distances(self, seqs1: list, seqs2: list) -> List[float]:
        pass
    
    def merge_data(self, other, sep_cols = False) -> pd.DataFrame:
        # !!! IMPORTANT: "other" must also be SequenceSimilarity object (couldnt compile)
        """
        Returns a merged Dataframe of self.p and another SequenceSimilarity's pep_data.
        If sep_cols=True, then the other SequenceSimilarity's columns will simply be appended
        to the returned DataFrame (self.p is unchanged). If False, results will be averaged.
        @NOTE: This is a super naive implementatoin -- expand this to make it more configurable
        @TODO: Take in *others as a list of arbitrarily many other SequenceSimilarities to compare
        """
        # must be same length binders -> so same peptides of interest
        this_data = self.p.copy()
        other_data = other.p.copy()
        sim_cols = self.sim_cols.copy()
        print(sim_cols)
        if sep_cols:
            suf = ("_"+self.bname, "_"+other.bname)
            out = this_data.merge(right=other_data, on=self.seq, suffixes=(suf))
            out = out.drop_duplicates()
            return out

        new_cols = ['{}_{}_{}'.format(col, self.bname, other.bname) for col in sim_cols]
        out_data = pd.DataFrame(index=this_data.index, columns=seq_cols + new_cols)
        out_data[seq_cols] = this_data[seq_cols]
        both = pd.concat([this_data[non_seq_cols],other_data[non_seq_cols]])
        out_data[new_cols] = both.groupby(both.index).mean()
        if 'sseq_matches' in self.cols or 'sseq_matches' in other.cols:
            both_match = self.p['sseq_matches'].append(other.pdata['sseq_matches'])
            both['sseq_matches'] = both_match
            out_data.join(both_match)
        return out_data
                    
        #@TODO Finish
                
                
    #-------------------miscellaneous methods------------------------------#
    
#     def get_kendalltau_corr_map(self) -> Tuple:
#         return stats.kendalltau(self.data['AA_MAP'][['Num']], self.data['AA_MAP'][['EIIP']])

from Bio.SubsMat import MatrixInfo as matlist

class SeqData:

    CONV = {
        'NUM': [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 4, 5, 5, 6, 6, 6, 7],
        'EIIP' : [0.0946, 0.0516, 0.0548, 0.0373, 0.0057, 0.0, 0.0, 0.0823,\
                 0.0829, 0.0941, 0.0036, 0.0761, 0.0198, 0.005, 0.1263, 0.0058, \
                 0.0371, 0.0242, 0.0959, 0.0829],
        'FNS' : ['Aromatic', 'Aromatic', 'Aromatic', 'Hydrophobic', 'Hydrophobic', \
                'Hydrophobic', 'Hydrophobic', 'Hydrophobic', 'Polar', 'Polar', 'Polar', \
                'Polar', 'Proline', 'Glycine', 'Charge (-)', 'Charge (-)', 'Charge (+)', \
                'Charge (+)', 'Charge (+)', 'Excluded'],
    }
    # @TODO: Import Biopython dicts of matrices, not read from .csvs
    MATR = {
        'PAM30': pd.read_csv('./src_data/pam30.csv', index_col=0).to_dict(),
        'BLOSUM': pd.read_csv('./src_data/BLOSUM.csv', index_col=0).to_dict(),
    }
    DIST = {
        'jaro_winkler': (lambda p1, p2: td.jaro_winkler.normalized_similarity(p1, p2)),
        'needleman_wunsch': (lambda p1, p2: td.needleman_wunsch.normalized_similarity(p1, p2)),
        'smith_waterman': (lambda p1, p2: td.smith_waterman.normalized_similarity(p1, p2)),
        'levenshtein': (lambda p1, p2: td.levenshtein.normalized_similarity(p1, p2))
    }
    
    def __init__(self, conv=CONV.keys(), matr=MATR.keys(), dist=DIST.keys()):
        
        self.conv = {cnv:self.CONV[cnv] for cnv in self.CONV.keys() if cnv in conv}
        self.matr = {mtr:self.MATR[mtr] for mtr in self.MATR.keys() if mtr in matr}
        self.dist = {dst:self.DIST[dst] for dst in self.DIST.keys() if dst in dist}
    