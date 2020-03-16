from seq_similarity import SeqData, SeqSimilarity

def main():
    DATA = SeqData()
    SEQS = [
        ('GRBP5','IMVTESSDYSSY'),
        ('M6','IMVTASSAYDDY')
    ]
    AA_COL = 'Sequences'
    PEP_PATH = './src_data/Sequence_data.csv' 
    data = {
        'grbp5_sim' : SeqSimilarity(SEQS[0], DATA, PEP_PATH, AA_COL),
        'm6_sim' : SeqSimilarity(SEQS[1], DATA, PEP_PATH, AA_COL)
    }
    dval = list(data.values())

if __name__ == '__main__':
    main()
