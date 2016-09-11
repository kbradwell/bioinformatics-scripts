=comment
Script extracts fasta sequences based on gene/protein lists.
Input: 
- Genome fasta file (gene sequences)
- List of gene names
Output:
Fasta file of selected genes.
=cut

use strict;
use warnings;

my $inputfasta = shift; # IFH
my $inputlist = shift; # IFH2

if(!defined($inputlist)) { 
    die "USAGE: perl $0 list-file list-file\n";
}

open IFH, "$inputfasta" or die "Can't open $inputfasta: $!\n";
open IFH2, "$inputlist" or die "Can't open $inputlist: $!\n";

my ($line,$defline,$seq,$total,$genename,$other,$counter);
my %seqHash;
my $first_flag = 1;
$total = 0;
while($line = <IFH>) { 
    chomp $line;
    $line =~ s/\n//g;
    $line =~ s/\r//g;
    
    if($line =~ /^>(.*)/) { 
	
	if($first_flag == 1) {
	    $first_flag = 0;
	} else { 
	    #$seq 
	    #print "Defline = $defline\n";
	    $seqHash{$genename} = $seq;
	    # Previous sequence is complete.
	    # Save sequence.
	    # 
	} 
	
	$defline = $1;
        # print $defline;
	if ($defline =~ /\s/) {
	($genename,$other) = split(/\s/, $defline);
        }
	else { 
		# no white space
        $genename = $defline;
        }
        # print "$genename\n";
	$seq = "";
	$total++;
    } else { 
	$seq .= $line;
    }
    
}

close IFH;

# Process Last sequence
$seqHash{$genename} = $seq;

$counter = 0;
foreach $genename (keys %seqHash) {
$counter ++;
}
#print "Total genes in hash:", $counter, "\n";
open IFH2, "$inputlist" or die "Can't open $inputlist: $!\n";

my ($listseq_id);

while($line = <IFH2>) { 
    chomp $line;
    $line =~ s/\n//g;
    $line =~ s/\r//g;
    
 	$listseq_id = $line;		
	
 	if(exists $seqHash{$listseq_id}) {
 	print ">$listseq_id\n$seqHash{$listseq_id}\n";
     }
     
}	

close IFH;
close IFH2;
