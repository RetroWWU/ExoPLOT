#!/usr/bin/env perl

#
# ==============================================================================
# This program counts the hits of a specific chromosome and coordinate
#
# Version 2
# ==============================================================================
#

use strict;

use vars qw/ %opt /;
use Getopt::Std;

$| = 1;

my $debug = 0;

# The maximum range hit in the chromosome files (actually it is about 350000)
my $range = 300000;

# Map sorted chromosome file -> offsets
my %alloffsets;

#
# ------------------------------------------------------------------------------
# Load the offsets of a sorted file
# ------------------------------------------------------------------------------
#
sub getOffsets {
	my ($file) = @_;

	if (!$alloffsets{$file}) {
		print "Load offsets for $file...\n" if ($debug);
		my $offsetfile = "$file.offset";
		open(F, "<", $offsetfile) or die "Cannot load offsets for $file!";
		my @offsets;
		while (my $line = <F>) {
			chomp($line);
			my @a = split(/\s+/, $line);
			push(@offsets, \@a);
		}
		close(F);
		$alloffsets{$file} = \@offsets;
		return \@offsets;
	}

	return $alloffsets{$file};
}

#
# ------------------------------------------------------------------------------
# Set the offset for start position
# ------------------------------------------------------------------------------
#
sub setOffset {
	my ($fh, $offsets, $start) = @_;

	$start = 0 if ($start < 0);
	my $offset = 0;
	for (my $i = 0; $i < @$offsets; $i++) {
		last if ($offsets->[$i][0] > $start);
		$offset = $offsets->[$i][1];
	}
	print "Offset: $offset\n" if ($debug);
	seek($fh, $offset, 0);
}

#
# ------------------------------------------------------------------------------
# Run the input on libpath
# ------------------------------------------------------------------------------
#
sub run {
	my ($input, $output, $libpath) = @_;

	my $starttime = time;
	my $lineno = 0;

	my ($name, $hitname, $nonhitname);
	my (%counthits, %countnonhits);
	
	# Loop over the input lines
	open(IN, "<", $input);
	while (my $line = <IN>) {
		chomp($line);

		$lineno++;
		print "==> Run request $lineno...\n";

		# Get the input for one line
		my ($rloc, $rchr, $rcoord3, $rcoord5);
		if ($line =~ m/(.*?)\s+(.*?)\s+(\d+)\s+(\d+)/) {
			($rloc, $rchr, $rcoord3, $rcoord5) = ($1, $2, $3, $4);
		}

		# Loop over all library folders
		my @dirs = grep {-d $_} <$libpath/*>;
		for my $dir (@dirs) {
			my ($count3, $count5, $countnonhit) = (0, 0, 0);

			my $file = "$dir/$rchr.txt";
			my $sortfile = "$file.sort";

			my $offsets = getOffsets($file);

			my $fh;
			open($fh, "<", $sortfile) or die "Cannot open sorted file for $file!";
			setOffset($fh, $offsets, $rcoord3 - $range);

			my $p = 0;
			my $count = 0;

			# Loop in the appropriate chromosome file
			while (my $line = <$fh>) {
				$count++;
				chomp($line);
				my @a = $line =~ /\s(\d+)-(\d+)/g;
				next if ($rcoord3 > $a[0] + $range * 2);
				last if ($rcoord3 < $a[0] - $range * 2);
				$p = $a[1];
				for my $i (1 .. (@a / 2)) {
					my $x_1 = $a[2 * $i - 1];
					my $x_2 = $a[2 * $i - 2];
					$count3++ if ($rcoord3 >= $x_2 && $rcoord3 <= $x_1);
					$count5++ if ($rcoord5 >= $x_2 && $rcoord5 <= $x_1);
					$countnonhit++ if ($rcoord3 >= $p && $rcoord5 <= $x_2);
					$p = $x_1;
				}
			}

			close($fh);

			print "Loops $count\n" if ($debug);;
			my $counthit = $count3 >= $count5 ? $count3 : $count5;

			$dir =~ m/E-MTAB-6814\.(\d+)/;
			$name = $1;
			$hitname = "${rloc}_hit";
			$nonhitname = "${rloc}_nonhit";
			$counthits{$name} += $counthit;
			$countnonhits{$name} += $countnonhit;
		}
	}

	my @keys = sort keys %counthits;
	open(OUT, ">", $output);
	print OUT "lib";
	for my $key (@keys) {
		print OUT "\t", $key;
	}
	print OUT "\n";
	print OUT "$hitname";
	for my $key (@keys) {
		print OUT "\t", $counthits{$key};
	}
	print OUT "\n";
	print OUT "$nonhitname";
	for my $key (@keys) {
		print OUT "\t", $countnonhits{$key};
	}
	print OUT "\n";
	close(OUT);

	close(IN);

	print "Duration ", (time - $starttime), "seconds.\n" if ($debug);;
}

#
# ------------------------------------------------------------------------------
# Print usage info
# ------------------------------------------------------------------------------
#
sub usage {
	print qq(USAGE
    count.pl -h
    count.pl -i input -o output -l libpath [-r range] [-d]

WHERE
    -h      - show help message
    -d      - debug messages
    input   - the input file
    output  - the output file
    libpath - the path to the library / database
    range   - the maximum difference in the library files between a
              start and end position in a line
              (default: $range)
);
	exit(1);
}

#
# ------------------------------------------------------------------------------
# main
# ------------------------------------------------------------------------------
#
sub main {
	getopts('hdi:o:l:r:', \%opt) or usage();
	usage() if ($opt{h});
	my $input   = $opt{i} or usage();
	my $output  = $opt{o} or usage();
	my $libpath = $opt{l} or usage();
	$range = $opt{r} if ($opt{r});
	$debug = $opt{d} if ($opt{d});
	run($input, $output, $libpath);
}

main();
