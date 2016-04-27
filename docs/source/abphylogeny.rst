AbPhylogeny
===========


Overview
--------

AbPhylogeny generates figure-quality phylogenetic trees from antibody sequence data.
Designed with the ability to color individual sequences by attribute, phylogenetic
trees can be drawn that accurately represent data from longitugindal samplings,
different sampling locations (peripheral blood, bone marrow, FNA, etc), or
categorical genetic characteristics like isotype.

AbPhylogeny can take input on any of three levels:

    1. FASTA-formatted sequence files
    2. FASTA-formatted multiple sequence alignment
    3. Newick-formatted tree file

If given sequence files, AbPhylogeny will perform multiple sequence alignment with
MUSCLE, build a tree file from the alignment with FastTree, and draw the tree figure. 
If given an alignment, AbPhylogeny will build the tree file and draw the figure. If
given a tree file, AbPhylogeny will simply draw the figure. In each case, AbPhylogeny
will save all intermediate files to the output directory, so intermediates can be used
to speed up multiple iterations on the same figure. This is especially helpful when 
trying multiple variations (colors, fontsizes, etc) of the same figure.




Examples
--------



