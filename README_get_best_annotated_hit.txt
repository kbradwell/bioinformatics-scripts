PURPOSE
Script gets the top 10 gb and ref hits from nr for each gene, performs blastdbcmd for each hit, then returns the top non hypothetical/unknown
protein hit for each gene with e-value <1e-10. 

USAGE
script <BLAST results for genes> <outfile>

INPUT 
BLAST file of nr database hits for genes. nr BLAST results input must contain the usual m8 output but with an extra end
column for database subject length, and be sorted by genename (first column), then e-value.

m8 output = "Q_ID", "S_ID", "IDENT", "ALIGN_LEN", "MISMATCHES", "GAPS", "Q_BEG", "Q_END", "S_BEG", "S_END", "E_VAL", "BIT_SCORE"

OUTPUT
Tab-delimited outout file containing gene name and annotation.
Output columns:
genename "Q_ID"
database hit ID "S_ID"
percent identity "IDENT"
alignment length "ALIGN_LEN"
mismatches "MISMATCHES"
gaps "GAPS"
query start "Q_BEG"
query end "Q_END"
subject start "S_BEG"
subject end "S_END"
e-value "E_VAL"
bit score "BIT_SCORE"
subject length "S_LEN"
percent subject length covered (S_END - S_BEG + 1 / S_LEN * 100) "PERC_S_LEN"
annotation of subject hit "DESCRIPTION"

NOTE
Before running, add the full PATH to the blastdbcmd and the nr database as strings to the beginning of the script.