=comment
Script to extracts genomic sequence from a contig fasta file based on start and end coordinates.
Input: 
- Genome fasta file (contig sequences)
- Table of coordinates to extract
Output:
Fasta file of sequences from the desired genomic coordinates.
=cut

use strict;
use warnings;

my $inputfasta = shift; # IFH
my $inputcoord = shift; # IFH2
my $outfile = shift; # OFH

if(!defined($outfile)) { 
    die "USAGE: perl $0 output file outfile\n";
}

open IFH, "$inputfasta" or die "Can't open $inputfasta: $!\n";
open IFH2, "$inputcoord" or die "Can't open $inputcoord: $!\n";

my ($line,$seq_id,$seq);
my %seqHash;

# save the contig sequences as a hash
while($line = <IFH>) { 
	chomp $line;
    $line =~ s/\n//g;
    $line =~ s/\r//g;
    # description line has to contain the contig name followed by whitespace
 	if($line =~ /^>(.*?)\s/) {

	$seq_id = $1;
	#print $seq_id,"\n";
	$seq = "";   # this gets overwritten each time, every time there's a new seq need to start w an empty string
    } else { 
	$seq .= $line;
	$seqHash{$seq_id} = $seq;
	#print $seqHash{$seq_id},"\n";
	}

}

foreach $seq_id (keys %seqHash) {
	#print "$seq_id\n$seqHash{$seq_id}\n"
}

# read the coordinates and print out the extracted sequences in FASTA format
open IFH2, "$inputcoord" or die "Can't open $inputcoord: $!\n";

open OFH, ">$outfile" or die "Error in opening $outfile for writing\n";

my ($coordseq_id,$sb,$se);

while($line = <IFH2>) { 
    chomp $line;
    
 	($coordseq_id,$sb,$se) = split(/\t/,$line);
	#print "$coordseq_id\t$sb\t$se\n";		
	
 	if(exists $seqHash{$coordseq_id}) {
 	my $sequence = $seqHash{$coordseq_id};
 		if($se>length($sequence) || $sb>length($sequence)) {
 			print "WARNING: coordinates do not match contig length\n";
 		}
 	my $exseq = substr($sequence, $sb - 1, $se-$sb +1);
 	print OFH ">$coordseq_id","_","$sb","_","$se","\n$exseq\n";
     }
     
}	

close IFH;
close IFH2;
close OFH;
