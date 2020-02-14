from seq_similarity import SequenceSimilarity

def main() -> None:
     

    # --------------------------- debug
    # Check to make sure df filter works
    print(grbp5_sim.df_filter_subseq('SVP', ind=0)) # it works!
    print(grbp5_sim[['PAM30', 'BLOSUM']])
    # ----------------------------

if __name__ == '__main__':
    main()