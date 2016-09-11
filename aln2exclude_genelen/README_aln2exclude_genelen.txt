PURPOSE
Prints out files that contain sequence lengths differing by >10% from the median length.

USAGE
script <file with lengths of each gene in an alignment>

INPUT
List of sequence lengths e.g. output from infoalign -only -alignlength < alignment filename> -out <alignment filenanme>.infoalign

OUTPUT
List of file names containing sequences that differ by >10% from the median length.