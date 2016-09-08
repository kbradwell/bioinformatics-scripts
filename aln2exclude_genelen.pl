=comment
Prints out files that contain sequence lengths differing by >10% from the median length.
Input: 
- List of sequence lengths e.g. output from EMBOSS (infoalign -only -alignlength < alignment filename> -out <alignment filenanme>.infoalign)
Output:
List of file names containing sequences that differ by >10% from the median length.
=cut

use strict;
use warnings;

use List::Util qw(max min sum); 
my $file = shift @ARGV;
open (FH, "< $file") or die "Can't open $file for read: $!";
my @a = <FH>;

# median of the lengths
sub median { $_[0]->[ @{$_[0]} / 2 ] }

my ($median,$min,$max,$align_reduction,$align_expansion);

# max and min lengths
$min     = min(@a);
$max     = max(@a);

$median=  median(\@a);

# print $median;

# percent that the min or max differs from the median
$align_reduction=$min/$median*100;
$align_expansion=$max/$median*100;

if($align_reduction<90 || $align_expansion>110){print $file."\n"}; #get a list of the files to exclude
