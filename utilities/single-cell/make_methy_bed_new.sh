#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

sub load_config {
    my %config;
    my $python_output = `python3 config_and_print.py`;
    die "Failed to run config_and_print.py: $!" if $?;
    
    for my $line (split /\n/, $python_output) {
        if ($line =~ /(\w+)=\'(.+)\'/) {
            $config{$1} = $2;
        } elsif ($line =~ /(\w+)=([\w:,]+)/) {
            $config{$1} = [split(/,/, $2)];
        }
    }
    return %config;
}

# Example usage in the main script
my %config = load_config();

# Read the filtered BAM list
my $filtered_list = $config{filtered_list};
open(my $fh, '<', $filtered_list) or die "Cannot open $filtered_list: $!";
my @prefixes;
while (my $line = <$fh>) {
    chomp($line);
    push @prefixes, $line;
}
close($fh);

# Configuration
my $java_path = "$config{software_directory}/java";
my $java_opts = "-Xmx18G";
my $java_classpath = join(":", (
    "$config{software_directory}/dnaaseUtils-0.14-jar-with-dependencies.jar",
    "$config{software_directory}/java-genomics-io.jar",
    "$config{software_directory}/igv.jar"
));
my $main_class = "main.java.edu.mit.compbio.utils.AlignMultiWigInsideBed";
my $bed_file = "$config{software_directory}/b37.common_chr.1Mb_interval.autosome.bed";
my $output_bed_base = "$config{output_directory}/b37.autosome.1Mb_interval.add_value.methy";
my @resolutions = @{$config{resolutions}};

# Open the output file for writing command strings
open(my $out_fh, '>', 'methy_summary.cmd.txt') or die "Cannot open methy_summary.cmd.txt: $!";

# List all .cov.b37.bw files in the specified directory
my @files = glob("$config{methy_directory}/*.cov.b37.bw");

# Create a hash for quick lookup of available files
my %available_files = map { basename($_, ".cov.b37.bw") => $_ } @files;

# Process each resolution
foreach my $res (@resolutions) {
    my ($resolution, $label) = split(/:/, $res);
    my $output_bed = "$output_bed_base.$label.bed.gz";
    
    # Process each prefix in the order of the filtered list
    foreach my $prefix (@prefixes) {
        if (exists $available_files{$prefix}) {
            my $cov = $available_files{$prefix};
            my $m = $cov;
            $m =~ s/cov/methy_count/;
            print $out_fh " -bigWig $m -useMean0 0 -regionMode 0 -bigWig $cov -useMean0 0 -regionMode 0 -res $resolution";
        }
    }
    
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

    print STDERR "Computation completed successfully for resolution $res\n";
}

close($out_fh);




























