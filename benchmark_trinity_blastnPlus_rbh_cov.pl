##!/usr/bin/env perl
#!/usr/bin/perl -w 
use strict;
use warnings;

#use FindBin;
use lib ("trinityrnaseq-2.0.5/PerlLib");
use Fasta_reader;
use Data::Dumper;





# run blast+ reciprocal: A->B, B->A
#blastn -query Trinity.fasta -db refTranscripts -out blastn.fmt6.txt  -evalue 1e-20 -dust no -task megablast -num_threads 2 -max_target_seqs 1 -outfmt 6

# analyze results
#analyze_blastPlus_topHit_coverage.pl blastn.fmt6.txt refTranscripts.fasta Trinity.fasta


#--------------------------------main---------------------------------------------------------------------------------------------------------

  


my $usage = "usage: $0 O2T.outfmt6.bn T2O.outfmt6.bn ori.fasta trinity.fasta [output_prefix=NameOfBlastFileHere] [verbose=0]\n
mainly for trinity run discovery mode for generating more splice variants while still want to get genuine assembly\n";

my $o2t_blast_out = $ARGV[0] or die $usage;
my $t2o_blast_out = $ARGV[1] or die $usage;
my $fasta_file_o = $ARGV[2] or die $usage;
my $fasta_file_t = $ARGV[3] or die $usage; 
my $output_prefix = $ARGV[4] || "$o2t_blast_out";
my $verbose = $ARGV[5] || 0;



my $best_hit_o= best_hits($o2t_blast_out); # hash ref of blast output for original fasta as query .
my $best_hit_t= best_hits($t2o_blast_out); # hash ref of blast output for trinity assembly fasta as query .
my %rbh_o2t;#rbh pair original fasta id as key.	
my %rbh_t2o;#rbh pair trinity assembly fasta id as key.	
for my $o (keys %{$best_hit_o}){
	my $best_o=$best_hit_o->{$o}{db_id};
	my $best_t=$best_hit_t->{$best_o}{db_id};
	next unless ($o && $best_t);
	if ($o eq $best_t){
		$rbh_o2t{$o}=$best_o;
		$rbh_t2o{$best_o}=$o;
	}
}


## get sequence length info, coverage=query(ori seq) aligned length / query seq length, TP,FN,FP for sp and sn calculation
           my ($specif,$sensit,$ori_count); 
	   my $tp=0;
	   my $fp=0;
	   my $fn=0;
            my $fasta_reader = new Fasta_reader($fasta_file_o);
            open (my $ofh, ">$output_prefix.id_cov_rbh_qlength") or die $!;
	    print $ofh join("\t", "#qseqid", "sseqid", "pident","pcover","qlen") . "\n";	
            while (my $seq_obj = $fasta_reader->next()) {
                my $acc = $seq_obj->get_accession();
		$ori_count++;
                if (exists $rbh_o2t{$acc}) {
                    my $sequence = $seq_obj->get_sequence();
                    my $seq_length = length($sequence);
               	    my $q_coverage= ($best_hit_o->{$acc}{query_match_len})/$seq_length;	 
		    print $ofh join("\t",$acc,$best_hit_o->{$acc}{db_id},$best_hit_o->{$acc}{percent_id},$q_coverage,$seq_length)."\n";	
		    $tp=$tp+ $q_coverage;
		    my $test;	
                }else{
			$fn++;
		}
            }

            my $fasta_reader_B = new Fasta_reader($fasta_file_t);
    	    while (my $seq_obj_B=$fasta_reader_B->next()){
		my $acc_B=$seq_obj_B->get_accession();
		$fp++ unless exists $rbh_t2o{$acc_B};	
	    }	
    ##  calculate sensitivity and specificity;
        
$specif=$tp/($tp+$fp);
$sensit=$tp/($tp+$fn);
my$tp_rate=$tp/$ori_count;
open ( FH,">$output_prefix.sn.sp.tpr.txt"); 
 print FH "$output_prefix\t$sensit\t$specif\t$tp_rate\n";
        
close $ofh;
close FH;    
    

#-------------------------------subroutines---------------------------------------------------------------------------------------------------
    # outfmt6:
    # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
sub best_hits{    
    print STDERR "-parsing blast output\n" if $verbose;
    my $file= shift;
    open (my $fh, $file) or die "Error, cannot open file $file";
    my %best_hits;
    while (<$fh>) {
        chomp;
        my @x = split(/\t/, $_);
        my $query_id = $x[0];
        my $db_id = $x[1];
        my $percent_id = $x[2];
        my $query_start = $x[6];
        my $query_end = $x[7];
        my $db_start = $x[8];
        my $db_end = $x[9];

        my $evalue = $x[10];
        my $bitscore = $x[11];
        
        if ( (! exists $best_hits{$query_id}) || ($evalue < $best_hits{$query_id}->{evalue}) ) {
            
            $best_hits{$query_id} = {	     db_id => $db_id,
                                             percent_id => $percent_id,
                                             query_start => $query_start,
                                             query_end => $query_end,
                                             db_start => $db_start,
                                             db_end => $db_end,
                                             evalue => $evalue,
                                             bitscore => $bitscore,
                                         
                                             query_match_len => abs($query_end - $query_start) + 1,
                                             db_match_len => abs($db_end - $db_start) + 1,
                                             
                                         };
        }



    }
    close $fh;
    return \%best_hits;	
}
        
