#!/usr/bin/perl -w
use Getopt::Long;




sub usage {
	
    print STDERR "\nUsage:\n";
    print STDERR "perl demultiplex_ecker_scWGBS.pl [Options] outputPrefix index_file.txt sample.1st_end.fastq.bz2 sample.2nd_end.fastq.bz2\n\n";
	
	exit(1);
}


my $output_prefix = $ARGV[0];
my $index_file = $ARGV[1];
my $first_end = $ARGV[2];
my $second_end = $ARGV[3];

usage() if ( scalar(@ARGV) < 4 );



my $linecount_within = 1;
my $linecount_global = 0;
my $linecount_find_match_global = 0;

my $seq_within_1st_line1 = "";
my $seq_within_2nd_line1 = "";

my %index_dict_r1;
my %index_dict_r2;
open(INDEX,"<$index_file") or die "can't open read index file $index_file: $!\n";
while(my $in=<INDEX>){
	chomp($in);
	$index_dict_r1{$in}="$output_prefix.$in.R1";
	$index_dict_r2{$in}="$output_prefix.$in.R2";
	open($index_dict_r1{$in}, ">$output_prefix.$in.R1.fastq") or die "can't open invert dups file $output_prefix.$in.R1.fastq: $!\n";
	open($index_dict_r2{$in}, ">$output_prefix.$in.R2.fastq") or die "can't open invert dups file $output_prefix.$in.R2.fastq: $!\n";

}
close(INDEX);

if($first_end=~/\.bz2/){
	open(FH1,"bzcat $first_end |") or die "can't open read file $first_end: $!\n";
	open(FH2,"bzcat $second_end |") or die "can't open read file $second_end: $!\n";
	
}elsif($first_end=~/\.gz/){
	open(FH1,"zcat $first_end |") or die "can't open read file $first_end: $!\n";
	open(FH2,"zcat $second_end |") or die "can't open read file $second_end: $!\n";
	
}else{
	open(FH1,"<$first_end") or die "can't open read file $first_end: $!\n";
	open(FH2,"<$second_end") or die "can't open read file $second_end: $!\n";
	
}


my $find_match=0;
my $key="";
while(my $seq1 = <FH1>, my $seq2 = <FH2>){
    
	chomp($seq1);
	chomp($seq2);
		if ($linecount_within == 1)
        {
            $seq_within_1st_line1 = "$seq1\n";
        	$seq_within_2nd_line1 = "$seq2\n";

        }
        elsif ($linecount_within == 2)
        {
			my @bases1 = split "",$seq1;
			
			for(my $i=0;$i<6;$i++){
				$key.=$bases1[$i];
			}
			if(exists $index_dict_r1{$key}){
				$find_match=1;
				my $trimmed="";
				for(my $i=6;$i<$#bases1;$i++){
					$trimmed.=$bases1[$i];
				}
				 $seq_within_1st_line1 .= "$trimmed\n";
				 $seq_within_2nd_line1 .= "$seq2\n";
			}else{
        		#print STDERR "$seq1\n";
        	}
        }
        elsif ($linecount_within == 3)
        {
        	$seq_within_1st_line1 .= "$seq1\n";
        	$seq_within_2nd_line1 .= "$seq2\n";
        }
        elsif ($linecount_within == 4)
        {
        	if($find_match==1){
        		my @bases1 = split "",$seq1;
        		my $trimmed="";
				for(my $i=6;$i<$#bases1;$i++){
					$trimmed.=$bases1[$i];
					
				}
				$seq_within_1st_line1 .= "$trimmed\n";
				 $seq_within_2nd_line1 .= "$seq2\n";
				 my $fh1=$index_dict_r1{$key};
				 my $fh2=$index_dict_r2{$key};
				 print $fh1 "$seq_within_1st_line1";
				 print $fh2 "$seq_within_2nd_line1";
				$linecount_find_match_global++; 
        	}
        	$linecount_global++;
        	$linecount_within = 0;
            $find_match=0;
            $key="";
            $seq_within_1st_line1 = "";
			$seq_within_2nd_line1 = "";
        }
        
   
     
	$linecount_within++;
	
}
close(FH1);
close(FH2);

foreach my $key(keys %index_dict_r1){
	close($index_dict_r1{$key});
	close($index_dict_r2{$key});
	`gzip -f $output_prefix.$key.R1.fastq`;
	`gzip -f $output_prefix.$key.R2.fastq`;
}

print STDERR "Find matched reads: $linecount_find_match_global\n";
print STDERR "Total reads: $linecount_global\n";

