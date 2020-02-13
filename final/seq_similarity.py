# @TODO Dynamically filter peptide set based on length(s) of input sequences of binders
#       i.e. 2 binders, one 11 AA long, one 13 AA long, each gets their own "subset" of the
#       full peptide lilst that can be compared to it. For any number of input sequences

# @TODO implement method for both similarities to M6 and GrBP5 to interrelate and act as their own feature set:
# i.e. if a peptide matches both peptides in temrs of sequence at some index, that should be important rather than
# having it equal to one matching only one

# @TODO: Maybe as originally planned implement binders inputted as dictionary of multiple values, each with separate
#      dataframes within the same class --> for big similarity calculations. Or just create multiple SequenceSimilarity instances

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

class SequenceSimilarity:
    '''
    Class that takes in a path to a list of amino acid sequences as well
    as any number of peptide sequences explicitly that are known to have
    a certain set of properties. Generates metrics for similarity for each
    peptide in path and returns domains AA sequence with high similarity
    '''

    def __init__(self, binder: str,   
                 data_paths: Dict,     #@TODO Make another .py file containing object version of necessary sim matrices/conversions, etc.
                 peps_path: str,       #      to reduce reliance on outside data
                 aa_col: str):         #-> The column in peps_path csv where sequences are held
        
        self.AA = list('FYWAVILMSTNQPGDEKHRC')
        self.conv = ["NUM", "EIIP", "FNS"]
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
        self.__read_similarity_data(data_paths) # !TO BE DEPRECATED! Just store data in class

        self.aa_col = aa_col
        self.columns = ['EIIP_Seq', 'NUM_Seq', 'PAM30', 'BLOSUM', 'RRM_SN', 'RRM_Corr', 'FNS_Seq']
        ## @TODO: Add "# of matching sseqs, cross entropy AA, cross entropy Num columns ?"
        self.AA_conv = lambda typ, pep: [self.AA_map[typ][AA] if AA in self.AA else 0 for AA in pep]
        self.get_sim = lambda p1, p2, t: sum([self.data[t][a1][a2] for a1 in p1 for a2 in p2])
        self.peps = pd.read_csv(peps_path)
        self.peps = self.peps[~self.peps[self.aa_col].str.contains("O")]
        self.peps.columns = [aa_col]
        self.peps_same_len = self.peps[self.peps[aa_col].str.len() == len(binder)]
        if len(self.peps_same_len) == 0:
            raise Exception("No peptides of same length as binder found")
        self.pep_data = self.peps_same_len.copy()
        for col in self.columns:
            self.pep_data[col] = None
        
        #self.get_similarities()

    def __read_similarity_data(self, data_path_dict) -> None:
        """
        Private method to store the paths of any data needed
        for similarity calcs and create Dataframes from them
        """
        self.data = dict.fromkeys(data_path_dict.keys())
        for data in self.data.keys():
            self.data[data] = pd.read_csv(data_path_dict[data], index_col=0)
            self.data[data] = self.data[data].to_dict()

            
    def _update_AA_conversion(self, conv_type: str = None) -> None:
        """
        Adds to the initially empty column values for conv_type (possible choices
        'EIIP' and 'Num' for now) the conversion of the AA sequence in the self.aa_col
        column a list representing its conversion
        """
        seqs = self.pep_data[self.aa_col]
        self.pep_data[conv_type+"_Seq"] = [self.AA_conv(conv_type, p) for p in seqs]

    def _update_matrix_similarity(self) -> None:
        """
        Just updates the similarity columns for the output similarity dataframe.
        Uses lambda helper function self.get_sim in __init__
        """
        for data in self.data.keys():
            sim = [self.get_sim(p, self.binder, data) for p in self.pep_data[self.aa_col]]
            self.pep_data[data] = np.interp(sim, (min(sim), max(sim)), (0,100))

    def _update_RRM_similarity(self) -> None:
        get_dft_from_eiip = lambda eiip: np.fft.rfft(eiip)[1:]
        get_cross_spectrum = lambda p1, p2: [x1*x2 for x1, x2 in zip(p1, p2)]
        pass
        #@TODO Finish
        
    # NOTE! Adds columns "Matching_sseqs" and "Num_matching" to output
    # Might be too unwieldy / unhelpful for output similarity data
    # if so, just comment out _update_matching_sseqs()
    def _update_matching_sseqs(self, single_match_weight:float = 0.5, weight:float = 2) -> None:
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
                            matches[i].pop()
                            num_matches[i] -= score(prev_match[0])
                            longest_match = (sseq, bin_i, len(sseq))
                        matches[i].append((sseq, bin_i))   
                        num_matches[i] += score(sseq)
                        prev_match = (sseq, bin_i)
        
        self.pep_data['sseq_matches'] = matches
        self.pep_data['weighted_matches'] = np.interp(num_matches, (min(num_matches), max(num_matches)), (0,100))
        
                    
    def remove_matching_sseqs_column(self) -> None:
        """
        Just removes the binder sseq pattern matches list column and number
        of matching (with/without) weighting if they exist
        """
        if self.pep_data.columns.contains(['sseq_matches', 'num_matches']):
            self.pep_data = self.pep_data.drop(columns=['sseq_matches', 'num_matches'])
    
    def update_similarities(self) -> None:
        '''
        Updates the similarity values whenever called (for now should be only once right
        after creating the object, ecept possibly if the Binding peptide is updated
        (should be handled automatically)
        '''
        self._update_AA_conversion(conv_type='EIIP')
        self._update_AA_conversion(conv_type='NUM')
        self._update_matrix_similarity()
        self._update_RRM_similarity()
        self._update_matching_sseqs()

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
        sseq = self.get_binder_subseq()
        filtered_data = [self.df_filter_subseq(ss,i) for (ss,i) in sseq if len(ss) >= min_length]
        return pd.concat(filtered_data)

    def get_bio_functions_by_loc(self):
        pass
    
    #------------------helper methods ------------------------------------#
                
                
    #-------------------miscellaneous methods------------------------------#

    def get_kendalltau_corr_map(self) -> Tuple:
        return kendalltau(self.data['AA_MAP'][['Num']], self.data['AA_MAP'][['EIIP']])

    #==========================notes===========================================#
    #TODO implement method to charcterize different locations on each peptide according to their
    #funtion from the number map. I.e. a (head -- that is, first 3 amino acids) of num seq
    # 112 can be considered "Primarily 1 (whatever that is)". The tail can be characteriszed in the
    # same way