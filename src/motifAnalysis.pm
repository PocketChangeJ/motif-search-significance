package motifAnalysis;

use strict;
use Exporter;
use File::Copy;
our @ISA= qw(Exporter);
our @EXPORT= qw(&getFastafromBED &getRandomSeq &runSearchMotif &removeDuplicateLine &countMotifs &countCoOccurence &annotMotifName);

use config;
use vars qw(%config );
*config=\%config::config;

###### This fonction takes a tabular file and turn the names of the transfac matrices
###### into their protein names.
sub annotMotifName{
	my ($motifFile, $input) = @_;

	### reading the file with the correspondance between the MA***** and the protein associated
	my %convert;
	open(MOTIF, "<".$motifFile) or die "Cannot open file $motifFile: $!";
	while(<MOTIF>){
		next unless(/^MOTIF/);
		chomp;
		my ($id, $motif) = $_ =~ /\s(MA.*)\s(.*$)/;
		if(defined $id and defined $motif){
			$convert{$id} = $motif;
		}
		
	}
	close(MOTIF);

	### Changing the name of the transfac matrix
	open(INPUT, "<".$input) or die "Cannot open file: $input: $!";
	open(OUT, ">".$input.".tmp") or die "Cannot create file: $input.tmp: $!";
	while(<INPUT>){
		chomp;
		my @tab = split "\t";
		my @modif = map { if(defined $convert{$_} ){ $_ =~ s/$_/$convert{$_}/}; $_  } @tab;
		print OUT join "\t", @modif;
		print OUT "\n";
	}
	close(OUT);
	close(INPUT);

	move($input.".tmp", $input);
}

sub countCoOccurence{
	my ($input, $annot, $out) = @_;

	my %data; ## Contain the number of motifs for each sequence
	my %motifList = &getMotifList($annot); ## Contain the list of all possible motifs

	open(FIMO, "<".$input) or die "Cannot open file $input: $!";
	while(<FIMO>){
		next if(/^#/);
		my @tab = split "\t";
		$data{$tab[1]}{$tab[0]}++;
	}
	close (FIMO);

	## Counting motif Coocurrence
	## A motif cooccure with another one if they are seen in the same sequence
	my %motifOcc;
	foreach my $seq (keys %data){
		foreach my $motif (keys %motifList){
			foreach my $key (keys %motifList){
				if($data{$seq}{$key} && $data{$seq}{$motif} && $key ne $motif){
					$motifOcc{$motif}{$key}++;
				}
				if($data{$seq}{$key} && $data{$seq}{$motif} && $key eq $motif){
					$motifOcc{$motif}{$key}=$motifOcc{$motif}{$key}+($data{$seq}{$key}-1) if($data{$seq}{$key} > 1);
				}				
			}
		}
	}

	open(OUT, ">".$out) or die "Cannot open output file $out: $!"; 
	## Printing out headers
	print OUT "Motifs";
	foreach my $motif (keys %motifList){
		print OUT "\t".$motif;
	}
	print OUT "\n";

	## Outputing matrix
	foreach my $motif (keys %motifList){
		print OUT $motif;
		foreach my $key (keys %motifList){
			if ($motifOcc{$motif}{$key}){
				print OUT "\t".$motifOcc{$motif}{$key};
			}
			else{
				print OUT "\t0";
			}
		}
		print OUT "\n";
	}
	close(OUT);
}


sub countMotifs{
	my ($input, $annot, $out) = @_;

	my %data; ## Contain the number of motifs for each sequence
	my %motifList = &getMotifList($annot); ## Contain the list of all possible motifs

	open(FIMO, "<".$input) or die "Cannot open file $input: $!";
	while(<FIMO>){
		next if(/^#/);
		my @tab = split "\t";
		$data{$tab[1]}{$tab[0]}++;
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

}


sub getFastafromBED{
	my ($bed, $genome_fa, $out_fa, $log) = @_;

	## Build up the command line
	my $cmdline = "$config{'bedtools'}/fastaFromBed -fi $config{$genome_fa} -bed $bed -fo - | \
		perl -ne 's/[:-]/_/g; print \$_' | gzip > $out_fa.gz";

	my $title = "Extracting nucleotide sequence out of the randomly selected coordinates";

	&runCmdLine($cmdline, $log, $title);

	return 1;

}

###### This function is built up to get the total list of motifs
###### from the *.meme files
###### The function return a hash table with all the matrix names as keys
sub getMotifList{
	my ($input) = @_;

	my %motifList;
	open(INPUT, "<".$input) or die "Cannot open file $input: $!";
	while(<INPUT>){
		next unless(/^MOTIF/);		
		chomp;
		## case of a jaspar matrix name
		my ($motif) = $_ =~ /MOTIF\s(MA.*)\s/;
		## case if there is only the name of the protein
		if($motif eq ""){
			($motif) = $_ =~ /MOTIF\s(.*)$/;
		}
		$motifList{$motif}=1 if(defined $motif);
	}
	close INPUT;

	return %motifList;

}

sub getRandomSeq {
	my ($bed, $genome_fa, $out_bed, $log) = @_;

	## Build up the command line
	my $cmdline = "$config{'bedtools'}/shuffleBed -i $bed -g $config{$genome_fa}.fai | cut -f1-3 | sort -k 1,1 -k 2,2n > $out_bed";

	my $title = "Selecting randomly selected sequences";

	&runCmdLine($cmdline, $log, $title);

	return 1;

}

sub removeDuplicateLine {
	my ($out_dir, $output, $log) = @_;

	my $cmdline = "head -1 $out_dir/fimo.txt > $out_dir/$output ; tail -n +2 $out_dir/fimo.txt | sort -u >> $out_dir/$output";

    my $title = "Removing duplicate lines from fimo results.";
	&runCmdLine($cmdline, $log, $title);
}

sub runCmdLine {
	my ($cmdline, $log, $title) = @_;

	## Run the command line and output into the log file
	print $log "\n##############\n" ;
	print $log "## Start analysis: ".`date`."\n";
	print $log "## $title. \n\n";
	print $log "$cmdline\n";
	print $log "#<--- Output: --------------------------------------------------\n";

	my $result = `$cmdline 2>&1`;
	print $log $result."\n";
	
	print $log "\n";

	print $log "#--------------------------------------------------------------->\n";
	print $log "## End of analysis: ".`date`."\n";

	return 1;
}

sub runSearchMotif{
	my ($fasta, $motifs, $out_dir, $log) = @_;
	
	## Build up the command line for fimo
	my $cmdline = "";
	if ($fasta =~ /.gz$/){
		$cmdline = "gunzip $fasta; "; 
		$fasta =~ s/\.gz$//;
	}
	$cmdline .= "$config{'fimo'} --thresh 1e-6 --oc $out_dir $motifs $fasta 2> $out_dir.nohup; gzip $fasta;";

	my $title = "Running the motif search analysis";
	&runCmdLine($cmdline, $log, $title);

	return 1;

}


1;
