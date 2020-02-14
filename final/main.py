from seq_similarity import SequenceSimilarity

def main() -> None:
     
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
    grbp5_sim = SequenceSimilarity(SEQS[0][1], SEQS[0][0], DATA_PATHS, PEP_PATH, AA_COL, True)
    m6_sim = SequenceSimilarity(SEQS[1][1], SEQS[1][0], DATA_PATHS, PEP_PATH, AA_COL, True)
    both_sep_sim = grbp5_sim.merge_data(other=m6_sim, sep_cols=True)
    both_avg_sim = grbp5_sim.merge_data(other=m6_sim, sep_cols=False)
    grbp5_sim.pep_data.to_csv('../out/grbp5_data_v2.csv')
    m6_sim.pep_data.to_csv('../out/m6_data_v2.csv')
    both_sep_sim.to_csv('../out/grbp5_m6_sep_v2.csv')
    both_avg_sim.to_csv('../out/grbp5_m6_avg_v2.csv')
    both_avg_sim.groupby()

    print(grbp5_sim.df_filter_subseq('SVP', ind=0)) # it works!
    print(grbp5_sim[['PAM30', 'BLOSUM']])
    # ----------------------------

if __name__ == '__main__':
    main()