# @TODO Dynamically filter peptide set based on length(s) of input sequences of binders
#       i.e. 2 binders, one 11 AA long, one 13 AA long, each gets their own "subset" of the
#       full peptide lilst that can be compared to it. For any number of input sequences

import pandas as pd
import numpy as np
from scipy.stats import kendalltau
import matplotlib.pyplot as plt
from typing import Set, Tuple, Dict

class SequenceSimilarity:
    '''
    Class that takes in a path to a list of amino acid sequences as well
    as any number of peptide sequences explicitly that are known to have
    a certain set of properties. Generates metrics for similarity for each
    peptide in path and returns domains AA sequence with high similarity
    '''

    def __init__(self, binders: Dict, data_paths: Dict, peps_path: str, aa_col: str):
        
        self.binder_dict = binders
        self.data_paths_dict = data_paths
        self.__update_binders()
        self.__update_similarity_data()

        self.peps = pd.read_csv(peps_path)
        self.peps.columns = [aa_col]
        self.peps_same_len = self.peps[self.peps[aa_col].str.len() == self.binder_len]
        self.aa_col = aa_col

        self.AA = set('LINGVEPHKAYWQMSCTFRD')
        self.sseq_set: Set[Tuple[int, str]] # full set of binder subseqs
        self.top_sseq: Set[Tuple[int, str]] # set of sub sequences w/ high simil
        self.top_seq: Set[str]              # set of peptides with high simil
    
    def __update_binders(self) -> None:
        '''
        Private method to set the binders stored by
        class and check to make sure they are of the
        same length as one another
        '''
        print(self.binder_dict)
        self.binders = [binder for binder in self.binder_dict.values()]
        try:
            self.binder_len = len(self.binders[0])
        except:
            print("Need at least one input binder")
        
        for binder in self.binders:
            if len(binder) != self.binder_len:
                print("All binders must be of same length")
                # @TODO Handle multiple lengths of binders, each differently
                #       lengthed binder is compared with different parts of the
                #       full peptide set of the same length

    def __update_similarity_data(self) -> None:
        """
        Private method to store the paths of any data needed
        for similarity calcs and create Dataframes from them
        """
        self.data_paths: Dict = self.data_paths_dict
        self.data = {data:pd.read_csv(self.data_paths[data], header=None) \
                     for data in self.data_paths.keys()}


    def df_filter_subseq(self, sub_seq: str, ind: int = None):
        '''
        Takes in a subsequence of equal or lesser length to
        peptides in class peptide dataframe and returns a dataframe
        containing only those peptides containing the sequence
        '''
        if not {*sub_seq}.issubset({*self.AA}):
            raise Exception('Invalid subsequence')
        if ind is None:
            return self.peps_same_len[peps_same_len[self.aa_col].str.contains(sub_seq)]

        data_filter = self.peps_same_len[self.aa_col].apply(
            lambda s: s[ind:len(sub_seq) == sub_seq])
        return self.peps_same_len.filter(data_filter)

    def get_sim_matrix(self, seq) -> pd.DataFrame:
        return self.data.filter


    def get_binder_subseq(self) -> pd.DataFrame:
        '''
        Generates all possible subsequences for binders
        provided in class constructor
        '''
        def gen_all_subseq(seq, sub_seq, i):
            if i == len(seq):
                if len(sub_seq) != 0:
                    yield(sub_seq)
                else:
                    gen_all_subseq(seq, sub_seq, sub_seq)
                gen_all_subseq(seq, sub_seq+[seq[i]], i+1)

        sseq = dict.fromkeys(self.binders)
        for binder in self.binders:
            sseq[binder] = [sseq for sseq in list(
                gen_all_subseq(binder, '', 0))]
        return sseq

    def get_PAM30_similarity(self) -> pd.DataFrame:
        '''
        Returns the PAM30 similarity of peptides to
        specified binder sequences
        @TODO: Automatically get perfect match and lowest match
        @TODO: Generalize for all binder sequences inputted
        '''
        raise NotImplementedError

    def get_BLOSUM_similarity(self) -> pd.DataFrame:
        raise NotImplementedError
    '''
    def get_RRM_SN_ratio(self):
        get_eiip_seq = lambda pep: list(map(lambda aa: self.AA_EIIP[aa], pep))
        get_dft_from_eiip = lambda eiip: np.fft.rfft(eiip)[1:]
    '''

    def get_kendalltau_corr_map(self) -> Tuple:
        return kendalltau(self.data['AA_MAP'][['Num']], self.data['AA_MAP'][['EIIP']])