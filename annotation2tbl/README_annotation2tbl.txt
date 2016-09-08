PURPOSE
Script to enable automation of whole genome annotation. Output is a .tbl file, which is necessary for tbl2asn script input (GenBank data deposition).
The script can be modified to suit particular genome annotation requirements.

USAGE
script <outfile> <mappingFile>
e.g. python annotation2tbl.py miscORG.tbl miscORG.tbl.map
(also have all relevant files containing annotation, with correct file extensions, in pwd).

INPUT
Various files containing different types of annotation information
Example files that all need to be present for this version of the script to work:
tRNA locations:
EXAMPLE_INPUT_annotation2tbl_tRNA.gff3
rRNA locations:
EXAMPLE_INPUT_annotation2tbl_rna.gff
protein-coding gene locations:
EXAMPLE_INPUT_annotation2tbl_genes.gff
Genes fasta:
EXAMPLE_INPUT_annotation2tbl_genes.fasta
EC #s and best nr BLAST hits:
EXAMPLE_INPUT_annotation2tbl_genes_supp_tbl
Intergenic repeats:
EXAMPLE_INPUT_annotation2tbl_ig-reps.outtable
Pseudogenes (all):
EXAMPLE_INPUT_annotation2tbl_pseudogenes.gff
EC number lookup file:
EXAMPLE_INPUT_annotation2tbl_genes.fasta.paths.detail
File with necessary prefixes:
EXAMPLE_INPUT_annotation2tbl_prefixes2use.txt, which contains Locus tag prefix, ProteinID prefix, Transcript ID prefix
e.g.
locus_tag       miscORG
protein_id      gnl|xyzVCU|
transcript_id   gnl|xyzVCU|mrna.

OUTPUT
.tbl file
Lookup table for gene names to locus tags and coordinates of other features to locus tags

