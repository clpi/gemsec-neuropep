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
from scipy.stats import kendalltau
import matplotlib.pyplot as plt
from typing import Set, Tuple, Dict, List
from scipy import interpolate
from scipy import stats
import textdistance as td

class SequenceSimilarity:
    '''
    Class that takes in a path to a list of amino acid sequences as well
    as any number of peptide sequences explicitly that are known to have
    a certain set of properties. Generates metrics for similarity for each
    peptide in path and returns domains AA sequence with high similarity
    '''

    def __init__(self, binder: str,
                 binder_name: str,
                 data_paths: Dict,     #@TODO Make another .py file containing object version of necessary sim matrices/conversions, etc.
                 peps_path: str,       #      to reduce reliance on outside data
                 aa_col: str,
                 dists: bool = False,
                 only_matching: bool = False):         #-> The column in peps_path csv where sequences are held
        
        # ---- setting data ---------
        self.AA = list('FYWAVILMSTNQPGDEKHRC')
        self.conv = ["NUM", "EIIP", "FNS"]
        self.conv_cols = [conv+'_Seq' for conv in self.conv]
        self.con_vals = dict.fromkeys(self.conv)
        self.con_vals['NUM'] = [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 4, 5, 5, 6, 6, 6, 7]
        self.con_vals['EIIP'] = [0.0946, 0.0516, 0.0548, 0.0373, 0.0057, 0.0, 0.0, 0.0823,\
                     0.0829, 0.0941, 0.0036, 0.0761, 0.0198, 0.005, 0.1263, 0.0058, \
                     0.0371, 0.0242, 0.0959, 0.0829]
        self.con_vals['FNS'] = ['Aromatic', 'Aromatic', 'Aromatic', 'Hydrophobic', 'Hydrophobic', \
                    'Hydrophobic', 'Hydrophobic', 'Hydrophobic', 'Polar', 'Polar', 'Polar', \
                    'Polar', 'Proline', 'Glycine', 'Charge (-)', 'Charge (-)', 'Charge (+)', \
                    'Charge (+)', 'Charge (+)', 'Excluded']
        self.AA_map = dict.fromkeys(self.conv)
        for conv in self.conv:
            self.AA_map[conv] = dict(zip(self.AA, self.con_vals[conv]))
        self.binder = binder   # to get binder_len, just use len(self.binders)
        self.bname = binder_name
        self.__read_similarity_data(data_paths) # !TO BE DEPRECATED! Just store data in class
        
        # ---- helper lambda functions
        self.AA_conv = lambda typ, pep: tuple(self.AA_map[typ][AA] if AA in self.AA else 0 for AA in pep)
        self.get_sim = lambda p1, p2, t: sum([self.data[t][a1][a2] for a1 in p1 for a2 in p2])
        
        # ---- setting up dataframes
        # @TODO: It looks like some sequences show up as duplicates in the overall sequence list
        #        OR some sequences are being overwritten / returned multiple times? Find out why
        #        this is and debug
        self.aa_col = aa_col
        self.sim_cols = ['PAM30', 'BLOSUM', 'RRM_SN', 'RRM_Corr']
        self.columns = [*self.sim_cols] + [*self.conv_cols]
        ## @TODO: Add "# of matching sseqs, cross entropy AA, cross entropy Num columns ?"
        
        self.peps = pd.read_csv(peps_path)
        self.peps.columns = [aa_col]
        self.peps = self.peps.drop_duplicates()
        self.peps = self.peps[~self.peps[self.aa_col].str.contains("O")]
        self.peps_same_len = self.peps[self.peps[aa_col].str.len() == len(binder)]
        if len(self.peps_same_len) == 0:
            raise Exception("No peptides of same length as binder found")
        self.pep_data = self.peps_same_len.copy()
        for col in self.columns:
            self.pep_data[col] = None
        self.update_similarities(use_distance=True)
        self.pep_match = self.get_df_with_binder_subseqs()
        if only_matching: self.pep_data = self.pep_match
        
    #----------------SET UP FUNCTIONS (void)-------------------------------#

    def __read_similarity_data(self, data_path_dict) -> None:
        """
        Private method to store the paths of any data needed
        for similarity calcs and create Dataframes from them
        """
        self.data = dict.fromkeys(data_path_dict.keys())
        for data in self.data.keys():
            self.data[data] = pd.read_csv(data_path_dict[data], index_col=0)
            self.data[data] = self.data[data].to_dict()
            
    def _update_AA_conversion(self) -> None:
        """if seq.signaltonoise(cross) < max(sn):
        Adds to the initially empty column values for conv_type (possible choices
        'EIIP' and 'Num' for now) the conversion of the AA sequence in the self.aa_col
        column a list representing its conversion
        """
        for conv in self.conv:
            if conv == 'NUM':
                out = []
                for p in self.pep_data[self.aa_col]:
                    num_list = self.AA_conv(conv, p)
                    s = ''
                    for n in num_list:
                        s += str(n)
                    out.append(s)
                self.pep_data[conv+"_Seq"] = out
            else:
                self.pep_data[conv+"_Seq"] = [self.AA_conv(conv, p) for p in self.pep_data[self.aa_col]]
    
    def _update_matrix_similarity(self) -> None:
        """
        Just updates the similarity columns for the output similarity dataframe.
        Uses lambda helper function self.get_sim in __init__
        """
        for data in self.data.keys():
            sim = [self.get_sim(p, self.binder, data) for p in self.pep_data[self.aa_col]]
            self.pep_data[data] = np.interp(sim, (min(sim), max(sim)), (0,1))
        

    def _update_RRM_similarity(self) -> None:
        """
        Uses the Resonant Recognition Model as described by Irena Cosic to 
        """
        get_dft_from_eiip = lambda ls: np.real(np.fft.rfft(ls))
        get_cross_spectrum = lambda p1, p2: [x1*x2 for x1, x2 in zip(p1, p2)]
        bnd_eiip = self.AA_conv('EIIP', self.binder)
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
        
        for pep in self.pep_data[self.aa_col]:
            seq_eiip =self.AA_conv('EIIP', pep)
            seq_dft = get_dft_from_eiip(seq_eiip)
            cross = get_cross_spectrum(seq_dft, bnd_dft)
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
        self.pep_data['RRM_Corr'] = dot_out
        self.pep_data['RRM_SN'] = sn_out
    
    # NOTE! Adds columns "Matching_sseqs" and "Num_matching" to output
    # Might be too unwieldy / unhelpful for output similarity data
    # if so, just comment out _update_matching_sseqs()
    def _update_matching_sseqs(self, single_match_weight:float = 1, weight:float = 1.5) -> None:
        """
        Returns a number as a new column representing the number of "matches" a peptide
        has for all possible subsequences for the binder inputted at a given index. For
        weighting=1, all matches are treated equally ('Y' at position 3 is treated equal
        to IMV at position 0) but lowering weighting lowers smaller-length matches
        """
        # @TODO Remove "duplicates" which occur at different matching indexes of binder
        # but are part of a larger pattern already recorded at an earlier index
        self.pep_data['sseq_matches'] = None
        self.pep_data['weighted_matches'] = None
        
        score = lambda s: (single_match_weight * 1) + (len(s)**weight)
        all_sseqs = self.get_binder_subseq()
        matches = list(); num_matches = list()
        for i, seq in enumerate(self.pep_data[self.aa_col]):
            matches.append([]); num_matches.append(0)
            longest_match = None
            for j, AA in enumerate(seq): #-----------goes thru each sequence
                prev_match: Tuple = None
                for k, (sseq, bin_i) in enumerate(all_sseqs): #---------for each AA in seq, then go through all sseq and indexes
                    if (bin_i == j) and (seq[j:len(sseq)+j] == sseq): #--- if AA matches sseq and index of pep == index of sseq:
                        if longest_match is not None:
                            if seq[longest_match[1]:longest_match[1]+longest_match[2]].find(seq[j:len(sseq)+j]) >= 0:
                                continue
                        prev_binder_i: int = all_sseqs[k-1][1]
                        prev_binder_sseq: str = all_sseqs[k-1][0]
                        if prev_binder_i == bin_i:
                            num_matches[i] -= score(matches[i].pop()[0])
                            longest_match = (sseq, bin_i, len(sseq))
                        matches[i].append((sseq, bin_i))   
                        num_matches[i] += score(sseq)
                        prev_match = (sseq, bin_i)
        
        self.pep_data['sseq_matches'] = [*matches]
        self.pep_data['weighted_matches'] = np.interp(num_matches, (min(num_matches), max(num_matches)), (0,1))
        self.sim_cols += ['weighted_matches']
        self.columns += ['sseq_matches', 'weighted_matches']
        
                    
    def _remove_matching_sseqs_column(self) -> None:
        """
        Just removes the binder sseq pattern matches list column and number
        of matching (with/without) weighting if they exist
        """
        for col in ['sseq_matches', 'weighted_matches']:
            if self.pep_data.columns.contains(['sseq_matches', 'weighted_matches']):
                self.pep_data = self.pep_data.drop(columns=['sseq_matches', 'weighted_matches'])
                self.columns.remove(col)
            self.sim_cols.remove('weighted_matches')
            
    def _update_distances(self, metrics: List = []) -> None:
        """
        Adds Hamming, Levenstein, etc. distance metrics for sequences in peptide list
        Metrics can be specified by name string in parameter
        """
        distances = {
            'jaro_winkler': (lambda p1, p2: td.jaro_winkler.normalized_similarity(p1, p2)),
            'needleman_wunsch': (lambda p1, p2: td.needleman_wunsch.normalized_similarity(p1, p2)),
            'smith_waterman': (lambda p1, p2: td.smith_waterman.normalized_similarity(p1, p2)),
            'levenshtein': (lambda p1, p2: td.levenshtein.normalized_similarity(p1, p2)),
        }
        dists = list(distances.keys()) if len(metrics)==0 else metrics
        data = pd.DataFrame(index=self.pep_data.index, columns=dists)
        for dist in dists:
            self.pep_data[dist] = [distances[dist](pep, self.binder) for pep in self.pep_data[self.aa_col]]
        self.columns += dists
        self.sim_cols += dists
        
    def _unpack_num_encoding(self) -> None:
        for i, num_tup in enumerate(list(self.pep_data['NUM_Seq'])):
            for j, num in enumerate(num_tup):
                self.pep_data['Nm_'+str(j)].iloc[i]
                
                
    def update_similarities(self, use_distance=False, metrics: List = []) -> None:
        '''
        Updates the similarity values whenever called (for now should be only once right
        after creating the object, ecept possibly if the Binding peptide is updated
        (should be handled automatically)
        '''
        self._update_AA_conversion()
        self._update_matrix_similarity()
        self._update_RRM_similarity()
        self._update_matching_sseqs()
        if use_distance:
            if len(metrics) == 0:
                self._update_distances()
            else:
                self._update_distances(metrics)
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
            return self.pep_data[self.pep_data[self.aa_col].str.contains(sub_seq)]
        return self.pep_data[self.pep_data[self.aa_col].str.find(sub_seq) == ind]

    def get_sim_matrix(self, seq) -> pd.DataFrame:
        return self.data.filter

    def get_binder_subseq(self) -> pd.DataFrame: #Each binder associated with list of (sseq, index)
        '''
        Generates all possible subsequences for binders. Returns ia list of tuples, where
        each entry of the list is of the form ({Sub Seq}, {Ind}) of binder
        '''
        self.sseqs = [(self.binder[i:j], i) for i in range(len(self.binder)) for j in range(i+1, len(self.binder)+1)]
        return self.sseqs


    def get_df_with_binder_subseqs(self, min_length: int = 0) -> Dict[str, pd.DataFrame]:
        '''
        Returns a filtered version of self.peps_same_len DataFrame containing only
        those rows with sequences which contain subsequences (of min length specified in parameter) 
        of the two binder sequences in the locations where they occur in the binders
        '''
        # @TODO Check that the output df is actually right
        sseq = self.get_binder_subseq()
        self.pep_match = pd.concat([self.df_filter_subseq(ss,i) for (ss,i) in sseq if len(ss) >= min_length])
        return self.pep_match
        
    
    def merge_data(self, other, sep_cols = False) -> pd.DataFrame:
        # !!! IMPORTANT: "other" must also be SequenceSimilarity object (couldnt compile)
        """
        Returns a merged Dataframe of self.pep_data and another SequenceSimilarity's pep_data.
        If sep_cols=True, then the other SequenceSimilarity's columns will simply be appended
        to the returned DataFrame (self.pep_data is unchanged). If False, results will be averaged.
        @NOTE: This is a super naive implementatoin -- expand this to make it more configurable
        @TODO: Take in *others as a list of arbitrarily many other SequenceSimilarities to compare
        """
        # must be same length binders -> so same peptides of interest
        this_data = self.pep_data.copy()
        non_seq_cols = self.sim_cols.copy()
        seq_cols = [self.aa_col] + [*self.conv_cols]
        other_data = other.pep_data.copy()
        if sep_cols:
            if len(list(this_data.columns)) == len(list(other_data.columns)):
                return this_data.merge(right=other_data, on=self.aa_col, suffixes=("_"+self.bname, "_"+other.bname))        
            else:
                raise Exception("Mismatched columns")

        new_cols = ['{}_{}_{}'.format(col, self.bname, other.bname) for col in self.sim_cols]
        out_data = pd.DataFrame(index=this_data.index, columns=seq_cols + new_cols)
        out_data[seq_cols] = this_data[seq_cols]
        both = pd.concat([this_data[non_seq_cols],other_data[non_seq_cols]])
        out_data[new_cols] = both.groupby(both.index).mean()
        if 'sseq_matches' in self.columns or 'sseq_matches' in other.columns:
            both_match = self.pep_data['sseq_matches'].append(other.pep_data['sseq_matches'])
            both['sseq_matches'] = both_match
            out_data.join(both_match)
        self.both = both
        return out_data
                    
        #@TODO Finish
                
                
    #-------------------miscellaneous methods------------------------------#

    def get_kendalltau_corr_map(self) -> Tuple:
        return kendalltau(self.data['AA_MAP'][['Num']], self.data['AA_MAP'][['EIIP']])
    
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