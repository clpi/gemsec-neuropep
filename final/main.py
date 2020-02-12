from seq_similarity import SequenceData, SequenceSimilarity

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
    AA_COL = 'AA_seq'
    PEP_PATH = 'src_data/Organized_NPS.csv'
    seq_data = SequenceData(SEQS, )

if __name__ == '__main__':
    main()