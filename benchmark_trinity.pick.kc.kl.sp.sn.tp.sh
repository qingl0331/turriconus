#! /bin/bash


set -x
LEFT_READS=$1
RIGHT_READS=$2
REF_BLAST_DB_NAME=$3
#TRINITY_BLAST_DB_NAME=$4
REFERENCE_FASTA_FILE_NAME=$4
KMER_COVERAGE=$5
KMER_LENGTH=$6
#SP_SN_OUT=$6

#-z means zero length
if [[ -z $LEFT_READS || -z $RIGHT_READS || -z $REF_BLAST_DB_NAME || -z $REFERENCE_FASTA_FILE_NAME || -z $KMER_COVERAGE || -z $KMER_LENGTH ]]; then
    echo "Usage: $0 LEFT_READS  RIGHT_READS  REF_BLAST_DB_NAME  REFERENCE_FASTA_FILE_NAME KMER_COVERAGE KMER_LENGTH"
    exit 1
#fi closes the if statement
fi
echo -e  "kmer_coverage\tsensitivity\tspecificity"
#make ref blastdb ->run Trinity ->make Trinity blastdb -> run blastn ( ref2trinity & trinity2ref)-> benchmark 
			makeblastdb -in $REFERENCE_FASTA_FILE_NAME -out $REF_BLAST_DB_NAME -dbtype nucl
                        rm -rf kc.$KMER_COVERAGE.kl.$KMER_LENGTH.Trinity
                	rm -rf kc.$KMER_COVERAGE.kl.$KMER_LENGTH.Trinity.Trinity.fasta
			Trinity --seqType fq --max_memory 30G --left $LEFT_READS --right $RIGHT_READS --CPU 23 --bflyHeapSpaceMax 5G --bflyCPU 12 --SS_lib_type RF --KMER_SIZE $KMER_LENGTH --min_kmer_cov $KMER_COVERAGE --min_contig_length 75 --output kc.$KMER_COVERAGE.kl.$KMER_LENGTH.Trinity --full_cleanup
# output: kc.$KMER_COVERAGE.kl.$KMER_LENGTH.Trinity.Trinity.fasta
			makeblastdb -in kc.$KMER_COVERAGE.kl.$KMER_LENGTH.Trinity.Trinity.fasta -out kc${KMER_COVERAGE}.kl${KMER_LENGTH} -dbtype nucl
			blastn -db kc${KMER_COVERAGE}.kl${KMER_LENGTH} -query $REFERENCE_FASTA_FILE_NAME -num_threads 22 -outfmt 6 -evalue 1e-3 -word_size 30 -reward 1 -penalty -1 -gapopen 1 -gapextend 2 -soft_masking true -dust yes -max_target_seqs 100 -out ref_vs_kc${KMER_COVERAGE}.kl${KMER_LENGTH}.bn
			blastn -db $REF_BLAST_DB_NAME -query kc.$KMER_COVERAGE.kl.$KMER_LENGTH.Trinity.Trinity.fasta -num_threads 22 -outfmt 6 -evalue 1e-3 -word_size 30 -reward 1 -penalty -1 -gapopen 1 -gapextend 2 -soft_masking true -dust yes -max_target_seqs 100 -out kc${KMER_COVERAGE}.kl${KMER_LENGTH}_vs_ref.bn
			benchmark_trinity_blastnPlus_rbh_cov.pl ref_vs_kc${KMER_COVERAGE}.kl${KMER_LENGTH}.bn kc${KMER_COVERAGE}.kl${KMER_LENGTH}_vs_ref.bn $REFERENCE_FASTA_FILE_NAME kc.$KMER_COVERAGE.kl.$KMER_LENGTH.Trinity.Trinity.fasta kc${KMER_COVERAGE}.kl${KMER_LENGTH}
