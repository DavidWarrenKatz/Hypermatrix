#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

# Configuration
my $java_path = "/home/dwk681/workspace/softwareFiles/java/jdk1.8.0_281/bin/java";
my $java_opts = "-Xmx18G";
my $java_classpath = "/home/dwk681/workspace/softwareFiles/java/FinaleMe/lib/dnaaseUtils-0.14-jar-with-dependencies.jar:/home/dwk681/workspace/softwareFiles/java/FinaleMe/lib/java-genomics-io.jar:/home/dwk681/workspace/softwareFiles/java/FinaleMe/lib/igv.jar";
my $main_class = "main.java.edu.mit.compbio.utils.AlignMultiWigInsideBed";
my $bed_file = "/home/dwk681/workspace/softwareFiles/hg19/b37.common_chr.500bp_interval.autosome.bed";
my $output_bed = "b37.autosome.500bp_interval.add_value.methy.GM_IMR90_153_samples.bed.gz";

# Open the output file for writing command strings
open(my $out_fh, '>', 'methy_summary.cmd.txt') or die "Cannot open methy_summary.cmd.txt: $!";

# List all .cov.b37.bw files in the specified directory
my @files = glob('/home/dwk681/workspace/cluster_cells_from_GSE189158_NOMe_HiC/filesFromCluster/bam/methylation/filter_low_qual/*.cov.b37.bw');

# Process each file to create the command strings
foreach my $cov (@files) {
    chomp($cov);
    my $m = $cov;
    $m =~ s/cov/methy_count/;
    print $out_fh " -bigWig $m -useMean0 0 -regionMode 0 -bigWig $cov -useMean0 0 -regionMode 0";
}

close($out_fh);

# Read the command strings from the file
my $cmd = `cat methy_summary.cmd.txt`;
chomp($cmd);

# Construct the full Java command
my $full_cmd = "$java_path $java_opts -cp \"$java_classpath\" $main_class $bed_file $output_bed $cmd";

# Print the command to be executed
print STDERR "Executing command: $full_cmd\n";

# Execute the Java command directly and capture error output
my $result = system("$full_cmd 2>java_error.log");

if ($result != 0) {
    # Print the error log content
    open(my $err_fh, '<', 'java_error.log') or die "Cannot open java_error.log: $!";
    while (my $line = <$err_fh>) {
        print STDERR $line;
    }
    close($err_fh);

    die "Failed to execute Java command\n";
}

print STDERR "Computation completed successfully\n";

