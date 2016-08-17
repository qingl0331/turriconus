#! /bin/bash


set -x
LEFT_READS=$1
RIGHT_READS=$2
REF_BLAST_DB_NAME=$3
#TRINITY_BLAST_DB_NAME=$4
REFERENCE_FASTA_FILE_NAME=$4
#SP_SN_OUT=$6

#-z means zero length
if [[ -z $LEFT_READS || -z $RIGHT_READS || -z $REF_BLAST_DB_NAME || -z $REFERENCE_FASTA_FILE_NAME ]]; then
    echo "Usage: $0 LEFT_READS  RIGHT_READS  REF_BLAST_DB_NAME  REFERENCE_FASTA_FILE_NAME"
    exit 1
#fi closes the if statement
fi
echo -e  "kmer_coverage\tsensitivity\tspecificity"
#make ref blastdb ->run Trinity ->make Trinity blastdb -> run blastn ( ref2trinity & trinity2ref)-> benchmark 
/data1/qli/tool/ncbi-blast-2.2.28+/bin/makeblastdb -in $REFERENCE_FASTA_FILE_NAME -out $REF_BLAST_DB_NAME -dbtype nucl
for i in {1..20}
do	
	for j in {10..32}
	do
                        rm -rf kc.$i.kl.$j.Trinity
                	rm -rf kc.$i.kl.$j.Trinity.Trinity.fasta
			/home/qli/tool/trinityrnaseq-2.0.5/Trinity --seqType fq --max_memory 30G --left $LEFT_READS --right $RIGHT_READS --CPU 23 --bflyHeapSpaceMax 5G --bflyCPU 12 --KMER_SIZE $j --min_kmer_cov $i --min_contig_length 75 --output kc.$i.kl.$j.Trinity --full_cleanup
# output: kc.$i.kl.$j.Trinity.Trinity.fasta
		/data1/qli/tool/ncbi-blast-2.2.28+/bin/makeblastdb -in kc.$i.kl.$j.Trinity.Trinity.fasta -out kc${i}.kl${j} -dbtype nucl
		/data1/qli/tool/ncbi-blast-2.2.28+/bin/blastn -db kc${i}.kl${j} -query $REFERENCE_FASTA_FILE_NAME -num_threads 22 -outfmt 6 -evalue 1e-3 -word_size 30 -reward 1 -penalty -1 -gapopen 1 -gapextend 2 -soft_masking true -dust yes -max_target_seqs 100 -out ref_vs_kc${i}.kl${j}.bn
		/data1/qli/tool/ncbi-blast-2.2.28+/bin/blastn -db $REF_BLAST_DB_NAME -query kc.$i.kl.$j.Trinity.Trinity.fasta -num_threads 22 -outfmt 6 -evalue 1e-3 -word_size 30 -reward 1 -penalty -1 -gapopen 1 -gapextend 2 -soft_masking true -dust yes -max_target_seqs 100 -out kc${i}.kl${j}_vs_ref.bn
		perl ~/bin/transcriptome_data_mining/benchmark_trinity_blastnPlus_rbh_cov.pl ref_vs_kc${i}.kl${j}.bn kc${i}.kl${j}_vs_ref.bn $REFERENCE_FASTA_FILE_NAME kc.$i.kl.$j.Trinity.Trinity.fasta kc${i}.kl${j}
	done
done
