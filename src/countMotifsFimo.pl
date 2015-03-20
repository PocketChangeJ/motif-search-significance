#! /usr/bin/env perl


use strict;
use warnings;
use Getopt::Long;



my $file;
my $out;
my $num_arg  = scalar @ARGV;

my $result = GetOptions(
	"fimo=s"      => \$file,
	"out=s"	      => \$out,
);

my $usage = <<END;

The script takes the fimo.txt file and count the number of motifs per sequences.
It outputs a matrix where each row is a sequence and each column is a motif.

Usage: $0 --fimo=FILENAME --out=FILENAME

  --fimo=FILENAME (mandatory - basically it is the fimo.txt file)
  --out=FILENAME (mandatory - output file name)

END

die $usage if (@ARGV);
die $usage if ( $num_arg == 0 );

my %data; ## Contain the number of motifs for each sequence
my %motifList; ## Contain the list of all possible motifs

open(FIMO, "<".$file) or die "Cannot open file $file: $!";
while(<FIMO>){
	next if(/^#/);
	my @tab = split "\t";
	$data{$tab[1]}{$tab[0]}++;
	$motifList{$tab[0]}=1;
}
close (FIMO);

open(OUT, ">".$out) or die "Cannot open output file $out: $!"; 
## Printing out headers
print OUT "Sequence Name";
foreach my $motif (keys %motifList){
	print OUT "\t".$motif;
}
print OUT "\n";

## printing the matrix :
## for each sequence, printing the number of motif occurence
foreach my $seq (keys %data){
	print OUT $seq;
	foreach my $motif (keys %motifList){
		if($data{$seq}{$motif}){
			print OUT "\t".$data{$seq}{$motif};
		}
		else
		{	
			print OUT "\t0";
		}
	}
	print OUT "\n";
}

close(OUT);





