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
from scipy import interpolate, stats, fftpack, signal
#import sci-kit learn
import textdistance as td
#import biopython as bio

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
                 min_length: int = 0,
                 dists: List = [],
                 only_match: bool = False):       
        
        # ---- setting data ---------
        
        #self.conv = ["NUM", "EIIP", "FNS"]
        self.d = {d: vars(data)[d] for d in list(vars(data).keys())}
        self.conv = list(sl['conv'].keys())
        self.aa_col = [aa_col]
        self.sim_cols = ['PAM30', 'BLOSUM', 'RRM_SN', 'RRM_Corr', 'weighted_matches']
        self.match_cols = ['sseq_matches', ]
        self.conv_cols = [conv+'_Seq' for conv in self.conv]
        self.cols = self.aa_col + self.sim_cols + self.conv_cols + self.match_cols
        self.AA_map = {cnv:dict(zip(self.AA, self.d['conv'][cnv])) for cnv in self.conv}
        self.bname, self.b = binder
        # finds all possible subsequences and locations in binder in tuple (sseq, ind)
        self.bsseq = [(self.b[i:j], i) for i, _ in enumerate(self.b) for j in range(i+1, len(self.b)+1)]
        
        # ---- helper lambda functions
        self.AA_conv = lambda typ, p: tuple(self.AA_map[typ][AA] if AA in self.AA else 0 for AA in p)
        self.sim_sum = lambda p1, p2, t: sum([self.d['matr'][t][a1][a2] for a1 in p1 for a2 in p2])
        
        self.p_og = pd.read_csv(p_path)
        self.__set_peps(aa_col, min_length, only_match)
        self.__update_similarities(dists)
        
        
        
    #----------------SET UP FUNCTIONS (void)-------------------------------#
    
    def __set_peps(self, aa_col:str, minlen: int = 0, matchonly: bool = False):
        self.p_og = self.p_og.drop_duplicates()
        self.p_og = self.p_og[~self.p_og[aa_col].str.contains("O")]
        self.p = self.p_og[self.p_og[aa_col].str.len() == len(self.b)]
        if not list(self.p):  raise Exception("No peptides of same length as binder found")
        self.pdata = pd.DataFrame(index=self.p.index, columns=self.cols)
        self.pdata_match = pd.concat([self.df_filter_subseq(ss,i) for (ss,i) in self.bsseq if len(ss) >= minlen])
        self.p_match = self.pdata_match[aa_col]
        if matchonly: self.p, self.pdata = self.pmatch, self.pdata_match

            
    def _update_AA_conversion(self) -> None:
        """if seq.signaltonoise(cross) < max(sn):
        Adds to the initially empty column values for conv_type (possible choices
        'EIIP' and 'Num' for now) the conversion of the AA sequence in the self.aa_col
        column a list representing its conversion
        """
        for conv in enumerate(self.conv):
            if conv == 'NUM':
                self.pdata[conv+"_Seq"] = [str([str(n) for n in self.AA_conv(conv, p)]) for p in self.p]
            else:
                self.pdata[conv+"_Seq"] = [self.AA_conv(conv, p) for p in self.p]
    
    def _update_matrix_similarity(self) -> None:
        """
        Just updates the similarity columns for the output similarity dataframe.
        Uses lambda helper function self.sim_sum in __init__
        """
        for m in list(self.d['matr'].keys()):
            sim = [self.sim_sum(p, self.b, m) for p in self.p]
            self.pdata[m] = np.interp(sim, (min(sim), max(sim)), (0,1))
        

    def _update_RRM_similarity(self) -> None:
        """
        Uses the Resonant Recognition Model as described by Irena Cosic to 
        """
        
        bnd_eiip = self.AA_conv('EIIP', self.b)
        bnd_dft = get_dft_from_eiip(bnd_eiip)
        sn = []
        do = []
        best_sn = None
        best_dot = None
        
        def signaltonoise(a, axis=0, ddof=0):
            a = np.asanyarray(a)
            m = a.mean(axis)
            sd = a.std(axis=axis, ddof=ddof)
            return np.where(sd == 0, 0, m/sd)
        
        for pep in self.p:
            seq_eiip =self.AA_conv('EIIP', pep)
            seq_dft = np.fft.rfft(seq_eiip)
            
            cross = signal.correlation(seq_dft, bnd_dft)
            dot = np.dot(seq_dft, bnd_dft)
            SN = np.mean(np.real(signaltonoise(cross, axis=None)))
            do.append(dot)
            sn.append(SN)
            if not sn and SN >= max(sn):
                best_sn = (pep, seq_dft)
            if not do or dot >= max(do):
                best_dot = (pep, seq_dft)
            
        sn_out = np.interp(sn, (min(sn), max(sn)), (0,1))
        dot_out = np.interp(do, (min(do), max(do)), (0,1))
        self.pdata['RRM_Corr'] = dot_out
        self.pdata['RRM_SN'] = sn_out
        
    def get_spectrums() -> pd.DataFrame:
        pass
    
    # NOTE! Adds columns "Matching_sseqs" and "Num_matching" to output
    # Might be too unwieldy / unhelpful for output similarity data
    # if so, just comment out _update_matching_sseqs()
    def _update_matching_sseqs(self, single_match_weight: float = 1, weight: float = 1) -> None:
        """
        Returns a number as a new column representing the number of "matches" a peptide
        has for all possible subsequences for the binder inputted at a given index. For
        weighting=1, all matches are treated equally ('Y' at position 3 is treated equal
        to IMV at position 0) but lowering weighting lowers smaller-length matches
        """
        # @TODO Remove "duplicates" which occur at different matching indexes of binder
        # but are part of a larger pattern already recorded at an earlier index
        self.pdata['sseq_matches'] = None
        self.pdata['weighted_matches'] = None
        
        score: float = lambda s: (single_match_weight * 1) + (len(s)**weight)
        matches: List = list(); nmatches = list()
        for i, seq in enumerate(self.p):
            matches.append(list()); nmatches.append(int())
            for j, AA in enumerate(seq): 
                lmatch: Tuple = None
                for k, (sseq, bin_i) in enumerate(self.bsseq):
                    in_seq = seq[j:len(sseq)+j]
                    if (bin_i == j) and (in_seq == sseq): 
                        if lmatch is not None:
                            if seq[lmatch[0]:lmatch[0]+lmatch[1]].find(in_seq) >= 0:
                                continue
                        if self.bsseq[k-1][1] == bin_i:
                            nmatches[i] -= score(matches[i].pop()[0])
                            lmatch = (bin_i, len(sseq))
                        matches[i].append((sseq, bin_i))
                        nmatches[i] += score(sseq)
        
        self.pdata['sseq_matches'] = matches
        self.pdata['weighted_matches'] = np.interp(nmatches, (min(nmatches), max(nmatches)), (0,1))
        self.sim_cols += ['weighted_matches']
        self.cols += ['sseq_matches', 'weighted_matches']
        
                    
    def _remove_matching_sseqs_column(self) -> None:
        """
        Just removes the binder sseq pattern matches list column and number
        of matching (with/without) weighting if they exist
        """
        if self.pdata.columns.contains(self.match_cols):
            self.cols.remove(self.match_cols)
            self.pdata = self.pdata.drop(columns=self.match_cols)
            
    def _update_distances(self, metrics: List = []) -> None:
        """
        Adds Hamming, Levenstein, etc. distance metrics for sequences in peptide list
        Metrics can be specified by name string in parameter
        """
        dists = metrics if metrics else list(self.d['dist'].keys())
        dist_vals = [[d['dist'][d](p, self.b) for p in self.p] for d in dists]
        self.pdata[dists] = dist_vals
        self.cols += dists
        self.sim_cols += dists
                
                
    def update_similarities(self, metrics: List = []) -> None:
        '''
        Updates the similarity values whenever called (for now should be only once right
        after creating the object, ecept possibly if the Binding peptide is updated
        (should be handled automatically)
        '''
        self._update_AA_conversion()
        self._update_matrix_similarity()
        self._update_RRM_similarity()
        self._update_matching_sseqs()
        if metrics is not None: 
            self._update_distances(metrics) if metrics else self._update_distances()
        # OPTIONAL
        # self._unpack_num_encoding()
        
    #----------MAIN CLASS FUNCTIONS (returns data) ------------------------#

    def df_filter_subseq(self, sub_seq: str, ind: int = None) -> pd.DataFrame:
        '''
        Takes in a subsequence of equal or lesser length to
        peptides in class peptide dataframe and returns a dataframe
        containing only those peptides containing the sequence
        '''
        if not {*sub_seq}.issubset({*self.AA}):
            raise Exception('Invalid subsequence')
        if ind is None:
            return self.pdata[self.p.str.contains(sub_seq)]
        return self.pdata[self.p.str.find(sub_seq) == ind]

    def sim_sum_matrix(self, seq) -> pd.DataFrame:
        return self.data.filter
        
    
    def merge_data(self, other, sep_cols = False) -> pd.DataFrame:
        # !!! IMPORTANT: "other" must also be SequenceSimilarity object (couldnt compile)
        """
        Returns a merged Dataframe of self.pdata and another SequenceSimilarity's pep_data.
        If sep_cols=True, then the other SequenceSimilarity's columns will simply be appended
        to the returned DataFrame (self.pdata is unchanged). If False, results will be averaged.
        @NOTE: This is a super naive implementatoin -- expand this to make it more configurable
        @TODO: Take in *others as a list of arbitrarily many other SequenceSimilarities to compare
        """
        # must be same length binders -> so same peptides of interest
        this_data = self.pdata.copy()
        non_seq_cols = self.sim_cols.copy()
        seq_cols = self.aa_col + self.conv_cols
        other_data = other.pep_data.copy()
        if sep_cols:
            if len(list(this_data.cols)) == len(list(other_data.cols)):
                d = this_data.merge(right=other_data, on=self.aa_col, suffixes=("_"+self.bname, "_"+other.bname))
                return d.drop_duplicates()
            else:
                raise Exception("Mismatched columns")

        new_cols = ['{}_{}_{}'.format(col, self.bname, other.bname) for col in self.sim_cols]
        out_data = pd.DataFrame(index=this_data.index, columns=seq_cols + new_cols)
        out_data[seq_cols] = this_data[seq_cols]
        both = pd.concat([this_data[non_seq_cols],other_data[non_seq_cols]])
        out_data[new_cols] = both.groupby(both.index).mean()
        if 'sseq_matches' in self.cols or 'sseq_matches' in other.cols:
            both_match = self.pdata['sseq_matches'].append(other.pdata['sseq_matches'])
            both['sseq_matches'] = both_match
            out_data.join(both_match)
        return out_data
                    
        #@TODO Finish
                
                
    #-------------------miscellaneous methods------------------------------#
    
#     def get_kendalltau_corr_map(self) -> Tuple:
#         return stats.kendalltau(self.data['AA_MAP'][['Num']], self.data['AA_MAP'][['EIIP']])
    
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
        'BLOSUM45': pd.read_csv('./src_data/BLOSUM.csv', index_col=0).to_dict(),
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
    
'''
import pandas as pd
from scipy.spatial.distance import euclidean, pdist, squareform


def similarity_func(u, v):
    return 1/(1+euclidean(u,v))

DF_var = pd.DataFrame.from_dict({"s1":[1.2,3.4,10.2],"s2":[1.4,3.1,10.7],"s3":[2.1,3.7,11.3],"s4":[1.5,3.2,10.9]})
DF_var.index = ["g1","g2","g3"]

dists = pdist(DF_var, similarity_func)
DF_euclid = pd.DataFrame(squareform(dists), columns=DF_var.index, index=DF_var.index)
'''

# @TODO Dynamically filter peptide set based on length(s) of input sequences of binders
#       i.e. 2 binders, one 11 AA long, one 13 AA long, each gets their own "subset" of the
#       full peptide lilst that can be compared to it. For any number of input sequences

# @TODO implement method for both similarities to M6 and GrBP5 to interrelate and act as their own feature set:
# i.e. if a peptide matches both peptides in temrs of sequence at some index, that should be important rather than
# having it equal to one matching only one

# @TODO: Maybe as originally planned implement binders inputted as dictionary of multiple values, each with separate
#      dataframes within the same class --> for big similarity calculations. Or just create multiple SequenceSimilarity instances
# @TODO: Implement method to apply relevant similarity scoring algorithms (distance metrics, ex) for
#       non-same-length peptides. Could even be stored alongside same length peptides
# @TODO Simplify some of the df filtering going on with Dataframe.where() or Dataframe.mask or Dataframe.query or isin
# @TODO Implement some patterns as being DEFINING features of binders -- ie. endi with letter 1 letter 2 letter 2 letter 1 or starting with IMVT always