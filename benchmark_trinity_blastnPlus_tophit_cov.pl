##!/usr/bin/env perl
#!/usr/bin/perl -w 
use strict;
use warnings;

#use FindBin;
use lib ("/home/qli/tool/trinityrnaseq-2.0.5/PerlLib");
use Fasta_reader;
use Data::Dumper;




# make blastable
#makeblastdb  -in refTranscripts.fasta -out refTranscripts -dbtype nucl

# run blast+
#blastn -query Trinity.fasta -db refTranscripts -out blastn.fmt6.txt  -evalue 1e-20 -dust no -task megablast -num_threads 2 -max_target_seqs 1 -outfmt 6

# analyze results
#analyze_blastPlus_topHit_coverage.pl blastn.fmt6.txt refTranscripts.fasta Trinity.fasta



  


my $usage = "usage: $0 blast+.outfmt6.txt query.fasta db.fasta.file [output_prefix=NameOfBlastFileHere] [verbose=0]\n
mainly for trinity run for generating more splice variants while still want to get full length of transcripts from tophits\n
original ref seq as query, trinity assembly as subject\n";

my $blast_out = $ARGV[0] or die $usage;
my $fasta_file_A = $ARGV[1] or die $usage;
my $fasta_file_B = $ARGV[2] or die $usage; 
my $output_prefix = $ARGV[3] || "$blast_out";
my $verbose = $ARGV[4] || 0;

main: {



    my %query_to_top_hit; # only storing the hit with the greatest blast score.
    my %sid_key_top_hit;	
    # outfmt6:
    # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    
    print STDERR "-parsing blast output: $blast_out\n" if $verbose;
    open (my $fh, $blast_out) or die "Error, cannot open file $blast_out";
    while (<$fh>) {
        chomp;
        my $line = $_;
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
        
        if ( (! exists $query_to_top_hit{$query_id}) || ($bitscore > $query_to_top_hit{$query_id}->{bitscore}) ) {
            
            $query_to_top_hit{$query_id} = { query_id => $query_id,
                                             db_id => $db_id,
                                             percent_id => $percent_id,
                                             query_start => $query_start,
                                             query_end => $query_end,
                                             db_start => $db_start,
                                             db_end => $db_end,
                                             evalue => $evalue,
                                             bitscore => $bitscore,
                                         
                                             query_match_len => abs($query_end - $query_start) + 1,
                                             db_match_len => abs($db_end - $db_start) + 1,
                                             
                                             line => $line,
                                         };
		$sid_key_top_hit{$db_id}=$query_id if ! exists $sid_key_top_hit{$db_id};			
        }



    }
    close $fh;

        
        ## get sequence length info, coverage=query(ori seq) aligned length / query seq length, TP,FN,FP for sp and sn calculation
           my ($fn,$fp,$specif,$sensit); 
	   my $tp=0;
            my $fasta_reader = new Fasta_reader($fasta_file_A);
            open (my $ofh, ">$output_prefix.id_cov_hit_length") or die $!;
	    print $ofh join("\t", "#qseqid", "sseqid", "pident","pcover","qlen") . "\n";	
            while (my $seq_obj = $fasta_reader->next()) {
                my $acc = $seq_obj->get_accession();
                if (exists $query_to_top_hit{$acc}) {
                    my $sequence = $seq_obj->get_sequence();
                    my $seq_length = length($sequence);
               	    my $q_coverage= ($query_to_top_hit{$acc}->{query_match_len})/$seq_length;	 
		    print $ofh join("\t",$query_to_top_hit{$acc}->{query_id},$query_to_top_hit{$acc}->{db_id},$query_to_top_hit{$acc}->{percent_id},$q_coverage,$seq_length)."\n";	
		    $tp=$tp+ $q_coverage;
		    my $test;	
                }else{
			$fn++;
		}
            }
            my $fasta_reader_B = new Fasta_reader($fasta_file_B);
    	    while (my $seq_obj_B=$fasta_reader_B->next()){
		my $acc_B=$seq_obj_B->get_accession();
		$fp++ unless exists $sid_key_top_hit{$acc_B};	
	    }	
    ##  calculate sensitivity and specificity;
        
$specif=$tp/($tp+$fp);
$sensit=$tp/($tp+$fn);
 print "sensitivity: $sensit\nspecificity: $specif\n";
        
close $ofh;
    
    
    exit(0);
    

}
