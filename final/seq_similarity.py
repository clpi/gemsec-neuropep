import pandas as pd
import numpy as np
from scipy.stats import kendalltau
import matplotlib.pyplot as plt
from typing import Set, Tuple, Dict, List

class SequenceSimilarity:
    '''
    Class that takes in a path to a list of amino acid sequences as well
    as any number of peptide sequences explicitly that are known to have
    a certain set of properties. Generates metrics for similarity for each
    peptide in path and returns domains AA sequence with high similarity
    '''

    def __init__(self, binders: Dict, data_paths: Dict, peps_path: str, aa_col: str):
        
        self.AA = set('LINGVEPHKAYWQMSCTFRD')
        self.binder_dict = binders
        self.data_paths_dict = data_paths
        self.__update_binders()
        self.__update_similarity_data()

        self.aa_col = aa_col
        self.columns = ['EIIP_Seq', 'Num_Seq', 'PAM30', 'BLOSUM', 'RRM_SN', 'RRM_Corr']

        self.peps = pd.read_csv(peps_path)
        self.peps.columns = [aa_col]
        self.peps_same_len = self.peps[self.peps[aa_col].str.len() == self.binder_len]

        self.peps_sl_sim = self.peps_same_len.copy()
        for col in self.columns:
            self.peps_sl_sim[col] = None
        self.__get_similarities()

        self.sseq_set: Set[Tuple[int, str]] # full set of binder subseqs
        self.top_sseq: Set[Tuple[int, str]] # set of sub sequences w/ high simil
        self.top_seq: Set[str]              # set of peptides with high simil
    
    def __update_binders(self) -> None:
        '''
        Private method to set the binders stored by class and check to 
        make sure they are of the same length as one another
        @TODO Allow for different lengthed binders, creates different lengthd
        peptide dataframes for analysis
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
        self.data = dict.fromkeys(self.data_paths.keys())
        for data in self.data_paths.keys():
            if data == "AA_MAP":
                self.data[data] = pd.read_csv(self.data_paths[data])
            else:
                self.data[data] = pd.read_csv(self.data_paths[data], header=None)

    def __get_similarities(self) -> None:
        self._get_AA_conversion(conv_type='EIIP')
        self._get_AA_conversion(conv_type='Num')
        self._get_BLOSUM_similarity()
        self._get_PAM30_similarity()
        self._get_RRM_SN_similarity()
        self._get_RRM_corr_similarity()

    def _get_AA_conversion(self, conv_type: str = None) -> None:
        AA_map = self.data["AA_MAP"][['AA', conv_type ]]
        def get_aa_conv(pep):
            new_seq = list()
            for AA in pep:
                if AA not in self.AA:
                    val = 0
                else:
                    val = AA_map.loc[AA_map['AA']==AA][conv_type].values[0]
                new_seq.append(val)
            return new_seq
        # get_aa_conv = lambda pep: map(lambda aa: AA_map.loc[AA_map['AA']==aa][conv_type].values[0], pep)
        self.peps_sl_sim[conv_type+"_Seq"] = self.peps_sl_sim[self.aa_col].apply(get_aa_conv)

    def _get_PAM30_similarity(self) -> pd.DataFrame:
        pass

    def _get_BLOSUM_similarity(self) -> pd.DataFrame:
        pass

    def _get_RRM_SN_similarity(self) -> None:
        pass

    def _get_RRM_corr_similarity(self) -> None:
        pass

    def df_filter_subseq(self, sub_seq: str, ind: int = None) -> pd.DataFrame:
        '''
        Takes in a subsequence of equal or lesser length to
        peptides in class peptide dataframe and returns a dataframe
        containing only those peptides containing the sequence
        '''
        if not {*sub_seq}.issubset({*self.AA}):
            raise Exception('Invalid subsequence')
        if ind is None:
            return self.peps_sl_sim[self.peps_sl_sim[self.aa_col].str.contains(sub_seq)]
        return self.peps_sl_sim[self.peps_sl_sim[self.aa_col].str.find(sub_seq) == ind]

    def get_sim_matrix(self, seq) -> pd.DataFrame:
        return self.data.filter

    def get_binder_subseq(self) -> Dict[str, List[Tuple[str, int]]]: #Each binder associated with list of (sseq, index)
        '''
        Generates all possible subsequences for binders. Returns in dict, where
        each binder corresponds to a dictionary of the index where it occurs,
        and the subsequence that occurs
        '''
        all_sseq = lambda s: [(s[i:j], i) for i in range(len(s)) for j in range(i+1, len(s)+1)]
        sseq = dict.fromkeys(self.binders)
        for binder in sseq.keys():
            sseq[binder] = all_sseq(binder)
        return sseq


    def get_df_with_binder_subseqs(self, min_length: int = 0) -> Dict[str, pd.DataFrame]:
        '''
        Returns a filtered version of self.peps_same_len DataFrame containing only
        those rows with sequences which contain subsequences (of min length specified in parameter) 
        of the two binder sequences in the locations where they occur in the binders
        '''
        data: Dict[str, List[pd.DataFrame]] = dict.fromkeys(self.binders)
        sseq = self.get_binder_subseq()
        for binder in self.binders:
            filtered_data = []
            for (ss, i) in sseq[binder]:
                if len(ss) >= min_length:
                    filtered_data.append(self.df_filter_subseq(ss, i))
            data[binder] = pd.concat(filtered_data)
        return data

    def get_kendalltau_corr_map(self) -> Tuple:
        return kendalltau(self.data['AA_MAP'][['Num']], self.data['AA_MAP'][['EIIP']])