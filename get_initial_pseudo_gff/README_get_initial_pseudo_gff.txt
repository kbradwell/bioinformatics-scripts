PURPOSE
Script takes maf and blastab format files from lastal (LAST, http://last.cbrc.jp/) search of genome assembly contigs and makes a gff of potential pseudogenes (without filtering for
called gene overlaps).

ALGORITHM
1) keep any nr database hits where >40% of the subject length is covered
2) select the hits with the longest genomic sequence, to leave a list of non-overlapping coordinates. If coords are tied (in terms of length), the hit with the best e-value is chosen
3) select only the genomic coodinates that display features of a pseudogene (frameshift or premature stop codon)

USAGE
script <maf infile> <blasttab infile> <GFF outfile>

INPUT
maf and blastab format files from lastal

OUTPUT
GFF file of putative pseudogenes. This can be used to subsequently filter out any coordinates overlapping called genes e.g. with intersect program http://bedtools.readthedocs.io/en/latest/content/tools/intersect.html.
GFF file columns:
Q_ID
"LAST search"
"DNA"
Q_BEG
Q_END
E_VAL
STRAND
"."
Last column, joined by "_": Q_LEN, S_LEN, PERC_S_LEN, IDENT, PSEUDO_FEATURE(F=frameshift,S=premature stop codon), S_ID, S_DESCRIPTION

NOTE
Add blastdbcmd and database (e.g. nr) PATH to the beginning of the script before running.

