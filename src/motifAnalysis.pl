#! /usr/bin/env perl

BEGIN {
	use File::Basename;
	push( @INC, dirname(__FILE__) );
}

use strict;
use warnings;
use threads;
use Getopt::Long;
use lib qw(.);
use motifAnalysis;

use config;
use vars qw(%config );
*config=\%config::config;

## Objectives : 


my $bed;
my $workingDir;
my $iterate = 1;
my $num_arg  = scalar @ARGV;

my $result = GetOptions(
	"bed=s"      => \$bed,
	"workingDir=s" => \$workingDir,
	"iterate=i" => \$iterate,
);

my $usage = <<END;

Usage: $0 --workingDir=DIRNAME --bed=FILENAME --iterate=INT 

  --workingDir=DIRNAME (mandatory - )
  --bed=FILENAME (mandatory - )
  --iterate=INT (mandatory - the number of iteration to get reliable control)

END

die $usage if (@ARGV);
die $usage if (! $result);
die $usage if ( $num_arg < 3 );

###### Check if working Dir exists and is writable
if( !-d $workingDir or !-w $workingDir){
	warn "Output directory doesn't exist or not writable. Check your parameters.\n";
	exit 0;
}

###### Getting the prefix of output file
my $prefix=$bed;
$prefix =~ s/\.(.*)//g;

##################################################################
################################# Creating a log file
## Creating the LOG directory if it doesn't exist
if(!-d "$workingDir/LOG"){
	mkdir "$workingDir/LOG";
}

## Getting the log file name
if ( -f "$workingDir/".basename($0).".log"){
	unlink "$workingDir/".basename($0).".log";
}

my $i = 1;
my $logName = "$workingDir/LOG/".basename($0).".$i.log";

while ( -f $logName ){
	$i++;
	$logName = "$workingDir/LOG/".basename($0).".$i.log";
}

symlink $logName, "$workingDir/".basename($0).".log";

open(my $log, ">".$logName) or die "Cannot create log file : $!";

##################################################################
################################# Creating output Dir
my $mainDir= "$workingDir/Main";
if(!-d $mainDir){
	mkdir $mainDir;
}

my $controlDir= "$workingDir/Controls";
if(!-d $controlDir){
	mkdir $controlDir
}

##################################################################
################################# MAIN
print $log "################ Starting Analysis\n";
print $log "##".`date`;
print $log "\n";


###### Running analysis on main file
my $thr;
$thr = threads->create(\&runMainAnalysis,"$workingDir/$bed", $mainDir, 'human_fa', $log);

###### Running analysis on controls
foreach my $i (1..$iterate){

	## checking if the number of threads is ok
	&waitThreads($config{'threads'}-1);

	$thr = threads->create(\&runSigAnalysis,"$workingDir/$bed", $controlDir, "shuff_$prefix\_$i", 'human_fa', $log);

}

###### Wainting for the threads to end
&waitThreads("0");

print $log "################ End of Analysis\n";
print $log "##".`date`;
print $log "\n";

##################################################################
################################# Functions

###### This function is used to :
## - wait for a thread to end (to run a new thread)
## - join a thread when run is over
## if the number of threads is higher than 0 and lower than the number of autorized threads
sub waitThreads {

	my ($nbThreads) = @_;

	my @threads;
	my @joinable;

	@threads = threads->list(threads::running);

	if($#threads > -1 && $#threads >= $nbThreads ){

		while( $#threads >= $nbThreads ){

			#print "Nb of threads:".$#threads."\n";
			@threads = threads->list(threads::running);
			sleep 1;
		}

		@joinable = threads->list(threads::joinable);
		if($#joinable > -1){
			foreach my $thr (@joinable){
				$thr->join();
				#print "one kill\n";
			}
		}
	}
}

###### This function is used to :
## - get fasta sequences out of chromosomal coordinates
## - run FIMO analysis
sub runMainAnalysis{

	my ($input, $outputDir, $genome, $log) = @_;

	my $prefix= basename($input);
	$prefix =~ s/\..*//g;

	###### Getting fasta sequences out of the bed files (Controls and input)
	if( ! -f "$outputDir/$prefix.fa.gz" ){
		&getFastafromBED($input, $genome, "$outputDir/$prefix.fa", $log);
	}

	###### Running motif analysis
	foreach my $motif ( keys %{$config{'motif'} } ){
		my $motifPrefix= basename( $motif );

		my $out_dir="$outputDir/FIMO-$motifPrefix-$prefix";

		if ( ! -f "$out_dir/fimo.txt" and -f "$outputDir/$prefix.fa.gz"){
			runSearchMotif("$outputDir/$prefix.fa.gz", $config{'motif'}{$motif}, "$out_dir", $log);
		}

		###### removing duplicate lines in FIMO results
		my $outputFimo="fimo_uniq.txt";
		if(! -f "$out_dir/$outputFimo"){
			&removeDuplicateLine($out_dir, $outputFimo, $log);
		}
	
	
		###### Count motifs
		my $output="$out_dir/CountOccurence-$prefix.tsv";
		if(! -f "$output"){
			&countMotifs("$out_dir/$outputFimo", $config{'motif'}{$motif}, $output);
			&annotMotifName($config{'motif'}{$motif}, $output);
		}
		
		###### Count motif Co-occurence
		$output="$out_dir/CountCoOccurence-$prefix.tsv";
		if(! -f $output){
			&countCoOccurence("$out_dir/$outputFimo", $config{'motif'}{$motif}, $output);
			&annotMotifName($config{'motif'}{$motif}, $output);
		}

		###### Cleaning unused data
    		unlink "$outputDir/$prefix.fa.gz";
		unlink "$outputDir/FIMO-$motifPrefix-$prefix/cisml.xml";
		unlink "$outputDir/FIMO-$motifPrefix-$prefix/fimo.gff";
		unlink "$outputDir/FIMO-$motifPrefix-$prefix/fimo.txt";
		`gzip $outputDir/FIMO-$motifPrefix-$prefix/$outputFimo`;

	}
	

	return 1;

}

###### This function is used to :
## do the same analysis as in runMainAnalysis function
## but first extract random coordinates first
sub runSigAnalysis{

	my ($input, $outputDir, $outputPrefix, $genome, $log) = @_;

	## Getting randomly selected sequences
	if( ! -f "$outputDir/$outputPrefix.bed" ){
		&getRandomSeq($input, $genome, "$outputDir/$outputPrefix.bed", $log);
	}

	&runMainAnalysis("$outputDir/$outputPrefix.bed", $outputDir, $genome, $log);

	return 1;

}




