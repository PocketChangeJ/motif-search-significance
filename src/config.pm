package config;

use strict;

use vars qw(%config);

our %config = (

	##### Tools - Rufus
	#"bedtools" => "/ifs/illumina/share/software/BEDTools/bedtools-2.17.0/bin",
	#"fimo" => "/ifs/illumina/share/software/meme490/bin/fimo",

	##### Tools - Furious
	"bedtools" => "/ifs/illumina/share/Utilities/softwareSL/BEDTools/bedtools-2.17.0/bin",
	"fimo" => "/ifs/illumina/share/Utilities/softwareSL/MEME/meme491/bin/fimo",

	"countMotifsFimo" => "~/Scripts/Motifs/countMotifsFimo.pl",

	##### Genome data
	"mouse_genome" => "/ifs/illumina/share/Genomes/Mus_musculus/mm9",
	"mouse_fa" => "/ifs/illumina/share/Genomes/Mus_musculus/mm9/Fasta/Unmasked/mm9.fa",

	"human_genome" => "/ifs/illumina/share/Genomes/Homo_sapiens/hg19",
	"human_fa" => "/ifs/illumina/share/Genomes/Homo_sapiens/hg19/Fasta/Unmasked/hg19.fa",

	##### Motif files
	"motif" => {
		#'selec' => "/ifs/illumina/slegras/S13051_Merienne/MotifFinding/Motifs/interestingMotifs.meme",
		'vertebrate' => '/ifs/illumina/slegras/S14098_Laurette/MotifAnalysis/Motifs/JASPAR_CORE_2014_vertebrates.meme'
	},

	##### Number of threads
	"threads" => 6,

);

1;
