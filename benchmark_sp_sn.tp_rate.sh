#! /bin/bash


set -x
REFERENCE_FASTA_FILE_NAME=$1

#-z means zero length
if [-z $REFERENCE_FASTA_FILE_NAME ]; then
    echo "Usage: $0 REFERENCE_FASTA_FILE_NAME"
    exit 1
#fi closes the if statement
fi
echo -e  "kmer_length\tsensitivity\tspecificity\ttrue_positive_rate"
for i in {10..32}
	do
		perl ~/bin/transcriptome_data_mining/benchmark_trinity_blastnPlus_rbh_cov.pl /data1/qli/conus_transcriptome/trinity_param_twist/bfly_seri/ref_vs_k${i}.bn /data1/qli/conus_transcriptome/trinity_param_twist/bfly_seri/k${i}_vs_ref.bn $REFERENCE_FASTA_FILE_NAME /data1/qli/conus_transcriptome/trinity_param_twist/bfly_seri/$i.Trinity.Trinity.fasta $i
	done


