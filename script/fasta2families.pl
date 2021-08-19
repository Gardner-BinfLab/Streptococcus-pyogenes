#!/usr/bin/perl

#Input: fasta file of proteins, database file (may be the same as the original) 
#Output: stockholm+hmm files, at most 1 per sequence, homologous sequences within the fasta are grouped into "families"
#Method:
#-foreach sequence, if not grouped into a family, run jackhmmer, build alignment and hmm for it and homologous found by jackhmmer

use strict;
use warnings;
use Getopt::Long;

my ($help, $faFile, $dbFile, $verb);
my $threshold      = 0.00001;
my $outDir         = '.'; 
my $iterations     = 3; 
my $nuc            = 0;
my $minClusterSize = 1;
&GetOptions( 
    "i|fa=s"       => \$faFile,
    "d|db=s"       => \$dbFile,
    "t|thresh=s"   => \$threshold,
    "o|outdir=s"   => \$outDir,
    "N|iterations" => \$iterations,
    "n|nuc=s"      => \$nuc,
    "m|minClust=s" => \$minClusterSize,
    "v|verbose"    => \$verb,
    "h|help"       => \$help 
    );

if( $help ) {
    &help();
    exit(1);
}
elsif ( not defined($faFile) ){
    #ERROR!!!
    print "ERROR: fasta file input is required!\n";
    &help();
    exit(1);
}
elsif ( not defined($dbFile) ){
    $dbFile = $faFile;    
}

#0. Index files if needed 
system("esl-sfetch --index $faFile")  if(not -s $faFile . '.ssi');
system("esl-sfetch --index $dbFile")  if(not -s $dbFile . '.ssi');

my %faIDs;
my $totInSeqs=0;
#1. fetch IDs for each input sequence 
open(FIN, "< $faFile") or die "Error: unable to open output file $faFile ($!)";
while (my $inf = <FIN>){
    if($inf =~ /^>\s*(\S+)/){
	$faIDs{$1} = 0;
	$totInSeqs++;
    }    
}
close(FIN);
#2. jackhmmer loop! 
my ($pid, $clusterNum) = ($$, 0, 0);
my @clusterCounts;
my %overlaps;
foreach my $faID (keys %faIDs){
    next if($faIDs{$faID} > 0); #seq is already in a cluster

    my $outTbl    = "";
    my $outRoot   = "$outDir/$pid-cluster$clusterNum";
    my $sFetch    = "esl-sfetch $faFile \47$faID\47 > $outDir/$pid\.fa";
    my $jackHmmer = "jackhmmer -N $iterations -o $outRoot\.jackhmmer  -A $outRoot\.stk --domtblout $outRoot\.domtblout -E $threshold --incE $threshold $outDir/$pid\.fa  $dbFile";
    my $hmmBuild  = "hmmbuild -o /dev/null  -n \42$outRoot\42 $outRoot\.hmm $outRoot\.stk";  
    system($sFetch   ) and die "FATAL: failed to run [$sFetch]!   \n[$!]\n";
    if($nuc==1){
	jacknhmmer($iterations, $outRoot, $threshold, "$outDir/$pid\.fa", $dbFile);
	$outTbl = "$outRoot\.tblout";
    }
    else {
	system($jackHmmer) and die "FATAL: failed to run [$jackHmmer]!\n[$!]\n";    
	system($hmmBuild ) and die "FATAL: failed to run [$hmmBuild]! \n[$!]\n";
	$outTbl = "$outRoot\.domtblout";
    }

    print "$outTbl\t" if ($verb);
    $clusterCounts[$clusterNum]=0;
    open(CIN, "< $outTbl");
    while (my $cin = <CIN>){
	next if ($cin =~ /^#/);
	if($cin =~ /^(\S+)\s+/){
	    if (defined($faIDs{$1}) && $faIDs{$1}==0){
		$faIDs{$1} = $clusterNum;
	    }
	    elsif(defined($faIDs{$1}) && $faIDs{$1}>0 && $faIDs{$1} != $clusterNum){ #recording overlapping families:
		my $oid = "$faIDs{$1}:$clusterNum";
		$overlaps{$oid} = 0 if not defined($overlaps{$oid});
		$overlaps{$oid}++;
		$faIDs{$1} = $clusterNum;
	    }
	    $clusterCounts[$clusterNum]++;
	}
    }
    close(CIN);
    printf "members=%d\n", $clusterCounts[$clusterNum] if ($verb);
    $clusterNum++;
}

#Concatenated files & print summaries:    
my $totSeqsInClusters=0;
print "cluster\tnumMembers\n" if($verb);
for(my $i=0; $i<$clusterNum; $i++){
    print "$i\t$clusterCounts[$i]\n" if($verb);
    $totSeqsInClusters+=$clusterCounts[$i];
    
    if($clusterCounts[$i] >= $minClusterSize){
	my $outRoot   = "$outDir/$pid-cluster$i";
	system("cat $outRoot\.hmm >> $outDir/$pid\.hmm") and print "WARNING: failed to run [cat $outRoot\.hmm >> $outDir/$pid\.hmm]! \n[$!]\n";
	system("cat $outRoot\.stk >> $outDir/$pid\.stk") and print "WARNING: failed to run [cat $outRoot\.stk >> $outDir/$pid\.stk]! \n[$!]\n";
    }
}

print "CHECK: 
number of input sequences     = [$totInSeqs]
number of clustered sequences = [$totSeqsInClusters] 
NB: some sequences are counted multiple times, local alignments\n\n\n" if($verb);



if($verb){
    system("esl-alistat $outDir/$pid\.stk") and print "WARNING: failed to run [esl-alistat $outDir/$pid\.stk]! \n[$!]\n";
    print "Overlapping families & number of overlaps:\n";
    foreach my $ov (sort keys %overlaps){
	my @ov = split(/:/, $ov);
	printf "$ov\t$overlaps{$ov}\t/%d\t/%d\tJaccard:%0.2f\tOverlapCoeff:%0.2f\n", $clusterCounts[$ov[0]], $clusterCounts[$ov[1]], $overlaps{$ov}/($clusterCounts[$ov[0]]+$clusterCounts[$ov[1]]-$overlaps{$ov}), $overlaps{$ov}/(min($clusterCounts[$ov[0]], $clusterCounts[$ov[1]]));	
    }
}




######################################################################
sub help {
    print STDERR <<EOF;

fasta2families.pl: Given a protein fasta file (and an optional
		   database), return stockholm files and hmms for each group of
                   homologous sequences in the fasta file.
                   I.e. cluster sequences using profile HMMs.
    
Usage:  fasta2families.pl <emblfile>
Options:
    -i|--fa         [fasta file]             Multi-protein fasta file input.
    -n|--nuc        [0|1]                    Input is nucleotide [1] or protein [0]. Default: [0]. 
    -d|--db         [sequence database file] Optional database for searching, if not provided the above fasta file is given.
    -t|--thresh     [float]                  Optional threshold for including sequences in model/output. Default: [-incE $threshold -E $threshold]
    -N|--iterations [integer]                Number of iterations of jackhmmer/nhmmer. Default: [-N $iterations]
    -m|--minClust   [integer]                Minimum cluster size, for saving concatenated alignments, HMMs and tables. Default: [$minClusterSize]
    -o|--outdir     [directory]              Directory to print outputs too, default is current dir. 
    -v|--verbose                             Verbose mode, print lots of stuff. 
    -h|--help                                Show this help.

Dependencies:
    hmmer-3.2
    easel
    
TODO:

    nhmmer support! 
    
EOF
}

#####################################################################
sub max {
  return $_[0] if @_ == 1;
  $_[0] > $_[1] ? $_[0] : $_[1]
}

sub min {
  return $_[0] if @_ == 1;
  $_[0] < $_[1] ? $_[0] : $_[1]
}

#####
#Iteratively run nhmmer, building new alignments & HMMs in each round.
sub jacknhmmer {
    my ($N, $filePrefix, $eThresh, $query, $db) = @_; 
    #jackhmmer -N 5 -o $outRoot\.jackhmmer  -A $outRoot\.stk --domtblout $outRoot\.domtblout -E $threshold --incE $threshold $outDir/$pid\.fa  $dbFile
    my $members=0; 
    for (my $i=0; $i<$N; $i++){
	
	my $nHmmer = "nhmmer -o $filePrefix\.jackhmmer -A $filePrefix\.stk --tblout $filePrefix\.tblout -E $eThresh $query $db";
	$query = "$filePrefix\.hmm";
	system($nHmmer) and die "FATAL: failed to run [$nHmmer]!\n[$!]\n";    
	
	unlink("$filePrefix\.hmm") if (-e "$filePrefix\.hmm");
	my $hmmBuild  = "hmmbuild -o /dev/null  -n \42$filePrefix\42 $filePrefix\.hmm $filePrefix\.stk";  
	system($hmmBuild ) and die "FATAL: failed to run [$hmmBuild]! \n[$!]\n";
	
	#Check number of significant hits, if it's unchanged, then exit loop...
	
    }
    
    return 0;
}







