package config;

use strict;

use vars qw(%config);

our %config = (

	##### Tools 
	"bedtools" => "", # <- Enter path for bedtools software
	"fimo" => "", # <- Enter path for fimo software

	"countMotifsFimo" => "countMotifsFimo.pl",

	##### Genome data
	"mouse_genome" => "/PATH/TO/Genomes/Mus_musculus/mm9", # <- enter path to the genome directory
	"mouse_fa" => "/PATH/TO/Genomes/Mus_musculus/mm9/Fasta/Unmasked/mm9.fa", # <- enter path to fasta file of genome

	"human_genome" => "/PATH/TO/Genomes/Homo_sapiens/hg19", # <- enter path to the genome directory
	"human_fa" => "/PATH/TO/Genomes/Homo_sapiens/hg19/Fasta/Unmasked/hg19.fa", # <- enter path to fasta file of genome

	##### Motif files
	"motif" => {
		'vertebrate' => '/PATH/TO/Motifs/JASPAR_CORE_2014_vertebrates.meme'
	},

	##### Number of threads
	"threads" => 6,

);

1;
