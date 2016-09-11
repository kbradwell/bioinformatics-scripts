PURPOSE
Script parses yn00 results.

USAGE
script <yn00 results>

INPUT
yn00 output file 

OUTPUT
A single "detailed" line with all the NG86 and yn00 results: NG86 omega, dN, dS, and YN00 S, N, t, kappa, omega, dN, dN SE, dS, dS SE
A single "brief" line with just the NG86 dN, dS and omega and yn00 dN, dS and omega.

COPYWRITE

Part of the code was modified from the BioPython Module Bio.Phylo.PAML._parse_yn00:
Copyright (C) 2011 by Brandon Invergo (b.invergo@gmail.com) 
This code is part of the Biopython distribution and governed by its 
license. Please see the LICENSE file that should have been included 
as part of this package. 

Biopython License Agreement

Permission to use, copy, modify, and distribute this software and its
documentation with or without modifications and for any purpose and
without fee is hereby granted, provided that any copyright notices
appear in all copies and that both those copyright notices and this
permission notice appear in supporting documentation, and that the
names of the contributors or copyright holders not be used in
advertising or publicity pertaining to distribution of the software
without specific prior permission.

THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
OR PERFORMANCE OF THIS SOFTWARE.