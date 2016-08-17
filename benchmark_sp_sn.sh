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
echo -e  "kmer_length\tsensitivity\tspecificity"
#make ref blastdb ->run Trinity ->make Trinity blastdb -> run blastn ( ref2trinity & trinity2ref)-> benchmark 
/data1/qli/tool/ncbi-blast-2.2.28+/bin/makeblastdb -in $REFERENCE_FASTA_FILE_NAME -out $REF_BLAST_DB_NAME -dbtype nucl
for i in {10..32}
	do
                rm -rf $i.Trinity
                rm -rf $i.Trinity.fasta
		/home/qli/tool/trinityrnaseq-2.0.5/Trinity --seqType fa --max_memory 30G --left $LEFT_READS --right $RIGHT_READS --CPU 26 --bflyHeapSpaceMax 5G --bflyCPU 19 --KMER_SIZE $i --min_contig_length 75 --min_per_id_same_path 99 --max_diffs_same_path 1 --max_internal_gap_same_path 3 --output $i.Trinity --full_cleanup
# output: $i.Trinity.Trinity.fasta
		/data1/qli/tool/ncbi-blast-2.2.28+/bin/makeblastdb -in /data1/qli/conus_transcriptome/trinity_param_twist/bfly_seri/$i.Trinity.Trinity.fasta -out k$i -dbtype nucl
		/data1/qli/tool/ncbi-blast-2.2.28+/bin/blastn -db k$i -query $REFERENCE_FASTA_FILE_NAME -num_threads 22 -outfmt 6 -evalue 1e-3 -word_size 30 -reward 1 -penalty -1 -gapopen 1 -gapextend 2 -soft_masking true -dust yes -max_target_seqs 100 -out ref_vs_k${i}.bn
		/data1/qli/tool/ncbi-blast-2.2.28+/bin/blastn -db $REF_BLAST_DB_NAME -query $i.Trinity.Trinity.fasta -num_threads 22 -outfmt 6 -evalue 1e-3 -word_size 30 -reward 1 -penalty -1 -gapopen 1 -gapextend 2 -soft_masking true -dust yes -max_target_seqs 100 -out k${i}_vs_ref.bn
		perl ~/bin/transcriptome_data_mining/benchmark_trinity_blastnPlus_rbh_cov.pl /data1/qli/conus_transcriptome/trinity_param_twist/bfly_seri/ref_vs_k${i}.bn /data1/qli/conus_transcriptome/trinity_param_twist/bfly_seri/k${i}_vs_ref.bn $REFERENCE_FASTA_FILE_NAME /data1/qli/conus_transcriptome/trinity_param_twist/bfly_seri/$i.Trinity.Trinity.fasta $i
	done

