from seq_similarity import SequenceSimilarity

def main() -> None:
    DATA_PATHS = {
        "BLOSUM":"src_data/BLOSUM.csv",
        "PAM30":"src_data/pam30.csv",
        "AA_MAP":"src_data/aa_chart.csv",
    }
    SEQS = {
        'GRBP5':'IMVTESSDYSSY',
        'M6':'IMVTASSAYDDY'
    }
    AA_COL = 'Sequences'
    PEP_PATH = 'src_data/Sequence_data.csv' 
    similarity = SequenceSimilarity(SEQS, DATA_PATHS, PEP_PATH, AA_COL)

    # --------------------------- debug
    # Check to make sure df filter works
    print(similarity.df_filter_subseq('SVP', ind=0)) # it works!
    # ----------------------------

if __name__ == '__main__':
    main()