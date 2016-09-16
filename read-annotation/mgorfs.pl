#!/usr/bin/env perl
#$ -S /usr/bin/perl
# ===========================================================================
# Author: Gene Tyson 08/29/07
#
# Modified: John Eppley 10/1/07
#
# Description:
#  runs metagene on a FASTA database and produces FASTA db of called ORFs
# ===========================================================================

use strict;
use warnings;
use lib qw(/RemotePerl/5.8.6);
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
my $modules = ["Bio::Tools::CodonTable"];

my $description = "runs metagene on a FASTA database and produces FASTA db of called ORFs\n";

### Process command line options ###
my ($opt_help, $opt_aa, $opt_nuc, $opt_mg, $opt_debug, $opt_unann, $fasta_file, $opt_prefix, $optabout);
my $version='mga';
my $old_mg=0;
Getopt::Long::Configure("pass_through", "auto_abbrev", "no_ignore_case");
my $result = GetOptions(
                        "input_file=s"      =>  \$fasta_file,
                        "output_prefix=s"      =>  \$opt_prefix,
                        "Amino_acids:s"      =>  \$opt_aa,
			"nucleotides:s"      =>  \$opt_nuc,
			"metagene:s"      =>  \$opt_mg,
			"unannotated:s"  => \$opt_unann,
			"Version=s"    => \$version,
			"verbose"    => \$opt_debug,
                        "help"      =>  \$opt_help,
			"about" => \$optabout,
	  );

if (defined $optabout) {
  print $description;
  exit(0);
}

# load modules after help, so we can get usage outside of perfect environment
if (not $opt_help) {
  foreach my $modname (@$modules) {
    print "loading $modname\n";
    eval "use $modname";
    die $@ if $@;
  }
}

die "Failed to parse command line options, use -h to see usage.\n" unless $result;
pod2usage({-exitval => 0, -verbose => 2}) if $opt_help;
pod2usage({-exitval => 1, -verbose => 2}) unless defined($fasta_file);

my $debug = defined($opt_debug);
print STDERR "Debugging on\n" if $debug;

my $output_prefix = (defined($opt_prefix)) ? $opt_prefix : "output";

## chose metagene version
my $command;
if ($version =~ m/^[oO]/ or $version eq 'metagene') {
    print STDERR "Using original metagene\n" if $debug;
    $command = 'metagene';
    $old_mg = 1;
} elsif ($version =~ m/^[mM]/) {
    print STDERR "Using mga -m\n" if $debug;
    $command = 'mga -m';
} elsif ($version =~ m/^[sS]/) {
    print STDERR "Using mga -s\n" if $debug;
    $command = 'mga -s';
} else {
    die "Unkown version string: $version\n Use -h option to see usage.\n";
}

### Prepare output streams ###
my %handles;
if (defined $opt_mg) {
    my $file = ($opt_mg ne "") ? $opt_mg : 
	$old_mg ? "$output_prefix.mge" : "$output_prefix.mga";
    open (MGOUT, ">$file");
    $handles{'mg'} = \*MGOUT;
    $opt_mg = 1;
}
if (defined $opt_nuc) {
    my $file = ($opt_nuc eq "") ? "$output_prefix.ffn" : $opt_nuc;
    $opt_nuc = 1;
    open (FFN, ">$file");
    $handles{'ffn'} = \*FFN;
}
if (defined $opt_aa) {
    my $file = ($opt_aa eq "") ? "$output_prefix.faa" : $opt_aa;
    $opt_aa = 1;
    open (FAA, ">$file");
    $handles{'faa'} = \*FAA;
}
if (defined $opt_unann) {
    my $file = ($opt_unann eq "") ? "$output_prefix.noPredGenes.fna" : $opt_unann;
    $opt_unann = 1;
    open (UNANN, ">$file");
    $handles{'unann'} = \*UNANN;
}

# bail if no output to be created
die "ERROR: Must choose at least one output method!\n(metagene output, predicted genes, or predicted proteins)\n\nUse -h to see usage.\n" unless ($opt_aa or $opt_nuc or $opt_mg or $opt_unann);

# open fasta file if we'll need it
if ($opt_nuc || $opt_aa || $opt_unann) {
    open (FASTA, $fasta_file) || die "Cannot open $fasta_file\n";
    $handles{'fasta'} = \*FASTA;
}


### Run metagene ###
$command .= " $fasta_file |";
print STDERR "Runing: '$command'\n" if $debug;
open (MG, $command) or die "Error running metagene!\n";

### metagene parser- stores gene and read information ###
my $read;
my $name;
my $count;
my %reads;

my %orfInfo;
while (my $a =<MG>) {
    $handles{'mg'}->print("$a") if $opt_mg;
    next unless ($opt_nuc or $opt_aa or $opt_unann);

    chomp $a;
    next unless $a;

    if ($a =~ m/^# /) {
    	# this is a header line
	if ($a =~ m/^# gc = (.+)/) {
	    # gc content line...ignore
	    next;
	} elsif ($a =~ m/^# (archaea|bacteria)\s*$/) {
	    # metagene domain line...ignore
	    next;
	} elsif ($a =~ m/^# self:\s+(\S+)/) {
	    # mga domain line...ignore
	    next;
	} elsif ($a =~ m/^# ([^\s]+)(.*)/) {
       	    #  header line...get data
	    # but first process last record
	    print STDERR "\n" if $debug;
            processRead(\%reads,$read,\%orfInfo,\%handles,$old_mg) if defined $read;

	    # now get data
            $name = $1; # full header line
	    $read = $1 . $2; # first word (before space) in header line
       	    print STDERR "reading MG data for $name\n" if $debug;
            $count = 0;
	    %orfInfo = ();
	}
    } else {
	# count orfs
	$count++;
	my $gene_name = $name."_GENE_".$count;
	$orfInfo{$gene_name} = $a;
       	print STDERR "... $gene_name " if $debug;
    }
}
print STDERR "\n" if $debug;

# do last read
processRead(\%reads,$read,\%orfInfo,\%handles,$old_mg) if defined $read;

# take in mg read info
# check skipped reads for passed read
# flip through input Fasta till read is found
# print genes
# store any skipped genes to reads hash
sub processRead {
    my $reads = shift;
    my $read = shift;
    my $orfInfo = shift;
    my $handles = shift;
    my $old_mg = shift;

    print STDERR " getting orfs for $read\n" if $debug;
    
    # check stored reads
    if (defined($reads->{'skipped'}{$read})) {
	print STDERR " matched stored read\n" if $debug;

	printOrfs($read,$reads->{'skipped'}{$read},$orfInfo,$handles,$old_mg);
	delete $reads->{'skipped'}{$read};

	return;
    }

    # look for read in fasta file
    my $thisRead;
    $thisRead = $reads->{'nextRead'} if defined($reads->{'nextRead'});
    #print STDERR "This read: $thisRead\n" if $debug;
    my $sequence = "";
    my $fasta = ${$handles}{'fasta'};
    while (my $line = <$fasta>) {
	chomp $line;
	next unless $line;

	if ($line =~ m/^>(.+)$/) {
	    #print $line . "\n";
	    # did we get the read we wanted?
	    if (defined($thisRead)) {
		if ($thisRead eq $read) {
		    printOrfs($read,$sequence,$orfInfo,$handles,$old_mg);
		    
		    # store read name for next time
		    $reads->{'nextRead'} = $1;

		    # done
		    return;
		} else {
		    print STDERR "  skipping $thisRead\n" if $debug;
		    $reads->{'skipped'}{$thisRead} = $sequence;
		}
	    }

	    # store read name for next time
	    $thisRead = $1;
	    $sequence = "";
	} else {
	    $sequence .= $line;
	}
    }

    # if we get here, we hit the end of the file...check last record
    if (defined($thisRead)) {
	if ($thisRead eq $read) {
	    printOrfs($read,$sequence,$orfInfo,$handles,$old_mg);
	} else {
	    print STDERR "  skipping last read: $thisRead\n" if $debug;
	    $reads->{'skipped'}{$thisRead} = $sequence;

	    print "WARNING: No match in fasta for metagene output record! ($read)\n";
	}
    } else {
	print "ERROR: Unexpected result! Is the fasta file corrupt or empty?";
    }

    # reached the end of fasta file
    delete $reads->{'nextRead'};
}	

if ($opt_nuc or $opt_aa or $opt_unann) {
# check for missed reads at end of fasta file
my $thisRead;
$thisRead = $reads{'nextRead'} if defined($reads{'nextRead'});
#print STDERR "This read: $thisRead\n" if $debug;
my $sequence = "";
my $fasta = $handles{'fasta'};
while (my $line = <$fasta>) {
    chomp $line;
    next unless $line;

    if ($line =~ m/^>(.+)$/) {
	# save previous record
	if (defined($thisRead)) {
	    print STDERR "  skipping $thisRead\n" if $debug;
	    $reads{'skipped'}{$thisRead} = $sequence;
	}
	# store read name for next time
	$thisRead = $1;
	$sequence = "";
    } else {
	$sequence .= $line;
    }
}

# at the end of the file...check last record
if (defined($thisRead)) {
    print STDERR "  skipping last read: $thisRead\n" if $debug;
    $reads{'skipped'}{$thisRead} = $sequence;
}

my $skippedReads = scalar keys %{$reads{skipped}};
if ($skippedReads > 0) {
    # Reads that were too short for metagene to consider
    print STDERR "WARNING: $skippedReads reads not found in MG output! (Too short?)\n" if $debug;

    if ($opt_unann) {
	#print skipped reads to unann (this should have been done in printorfs)
	foreach my $read (keys(%{$reads{skipped}})) {
	    $handles{'unann'}->print(">$read\n");
	    my $seq = $reads{skipped}{$read};
	    $handles{'unann'}->print("$seq\n");
	}
    }
}
}

# close Handles
close MG;
close MGOUT if $opt_mg;
close FAA if $opt_aa;
close FFN if $opt_nuc;
close UNANN if $opt_unann;
close FASTA;


sub printOrfs {
    my $read = shift;
    my $read_sequence = shift;
    my $orfInfo = shift;
    my $handles = shift;
    my $old_mg = shift;

    if ((scalar keys %{$orfInfo}) == 0) {	
	# no genes!
	if ($opt_unann) {
	    print STDERR "saving sequence of ORF-less read: $read\n" if $debug;
	    # print to unann file is user asked for it
	    $handles->{'unann'}->print(">$read\n");
	    $handles->{'unann'}->print("$read_sequence\n");
	}
	# nothing else to do
	return;
    }

    print STDERR "printing orfs for $read\n" if $debug;

    foreach my $gene (sort keys %{$orfInfo}) {	
	my @orf = split /\t/, $orfInfo{$gene};
# mge orf: 'start','end','strand','frame','score','wholeness (complete or partial)'	
# mga orf: 'name','start','end','strand','frame','wholeness (complete or partial)','score'
	# new version (mga) has extra field at start
	print STDERR "Shifting" if $debug && !$old_mg; 
	shift(@orf) unless $old_mg;
	my $start = $orf[0];
	my $stop = $orf[1];
	my $length = $stop - $start + 1;
	my $sense = $orf[2];
	my $frame = $orf[3];
	# we could adjust for mg versions, but it's not critical
	my $header = ">$gene start:$start end:$stop $orf[4] $orf[5]";	
	print STDERR "$header\n" if $debug;

	if ($sense eq "-") {
	    my $nt_seq = substr($read_sequence, $start - 1, $length - $frame);
	    my $rev_nt_seq = reverse $nt_seq;
	    $rev_nt_seq =~ tr/ATGCatgc/TACGtacg/;
	    $nt_seq = $rev_nt_seq;
	    ${$handles}{'ffn'}->print("$header\n$nt_seq\n") if $opt_nuc;
	    if ($opt_aa) {
		my $aa_seq = translate($nt_seq);
		$handles->{'faa'}->print("$header\n$aa_seq\n");
	    }
	} else {
	    my $nt_seq = substr($read_sequence, $start - 1 + $frame, $length - $frame);
	    $handles->{'ffn'}->print("$header\n$nt_seq\n") if $opt_nuc;
	    if ($opt_aa) {
		my $aa_seq = translate($nt_seq);
		$handles->{'faa'}->print("$header\n$aa_seq\n");
	    }
	}
    }
}

sub translate {
    my $myCodonTable = Bio::Tools::CodonTable->new();
    my $input = shift;
    my $length = length ($input);
    my $aa_seq;
    my $a;
    for ($a = 0; $a <= $length; $a = $a+3) {
	my $codon = substr($input, $a, 3);
	my $aa = $myCodonTable->translate($codon);
	$aa_seq .= $aa;
    }
    return($aa_seq);
}
exit;

__END__

=head1 NAME

B<mgorfs.pl> - call ORFs wih metagene and produce list of ORF sequences

=head1 SYNOPSIS

mgorfs.pl -i FASTA_FILE [options]

=head1 ARGUMENTS

=head1 OPTIONS

At least one output type must be specified. Options can be abbreviated with first letter.

=over 2

=item B<--input_file>

Filename of a FASTA list of sequences to pass through metagene

=item B<--output_prefix>

A filename prefix (can include path information) to write output files to. Defaults to "output". 
I.E. if no output prefix set, amino acid output is written to "output.faa" in the working directory.

=item B<--Amino-acids>

Output ORFs database as sequences of amino acids (> prefix.faa), a value given with this option (e.g. mgorfs -i file.fastsa -a file.faa) will override the output prefix and write the amino acid file to the indicated value.

=item B<--nucleotides>

Output ORFs database as sequences of nucleotides (> prefix.ffn)

=item B<--unann>

Output FASTA list of reads with no predicted genes (> prefix.noPredGenes.fna)

=item B<--metagene>

Write raw metagene output to file (> prefix.mge)

=item B<--Version>

         Use one of three different Versions/Modes (default: m): 
            'o', 'original', or 'metagene' will run the original version
            'm', 'mga', or 'multi' will run the new version in default mode
            's', or 'single' will run the new version in single species mode


=item B<--help>

Prints this message.

=back

=head1 DESCRIPTION

This script will submit the given FASTA list of sequences to metagene and capture the output.
The output will be used to produce a list of ORF sequences. These can be written out as
amino acid sequences or nucleotide sequences (or both).

=head1 EXIT CODES

This script returns 0 on success and a non-zero value on errors.

=head1 BUGS

Please report them to <jmeppley@mit.edu>

=head1 COPYRIGHT

2007 MIT (?)

=cut
