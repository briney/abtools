AbFinder
========


Overview
--------

AbFinder provides methods to mine large datasets of antibody sequences
to rapidly identify sequences with high identity to known antibody
sequences.

Given a MongoDB database and collection, AbFinder computes identity
between one or more 'standard' sequences and all seqeunces in the collection.
Default output is an identity/divergence plot, a hexbin plot of germline
identity (X-axis) and identity to the standard sequence (Y-axis). AbFinder
also updates MongoDB with identity information so that standard identities
can be used in subsequent queries.

Examples
--------

To run, AbFinder needs a MongoDB database and collection, an output directory,
and a FASTA-formatted file of standard sequences::

    $ abfinder -d MyDatabase -c MyCollection -s standards.fasta -o /path/to/output/

Omitting the collection results in AbFinder iteratively processing each collection
in the database.  By default, AbFinder assumes that the standard file contains
amino acid sequences. If you would like to compute nucleotide identity instead,
you can indicate your preference with the ``--nucleotide`` option::

    $ abfinder -d MyDatabase -s standards.fasta -o /path/to/output/ --nucleotide

AbFinder also assumes that the standard file contains heavy chain sequences, and only
heavy chain sequences from MongoDB will be used for comparison. To compare sequences
of a different chain (options are ``'heavy'``, ``'kappa'``, and ``'lambda'``), use
the ``--chain`` option::

    $ abfinder -d MyDatabase -s standards.fasta -o /path/to/output/ --chain kappa

If you do not plan on using the identity scores for any sort of downstream analysis,
you can save some time and skip the MongoDB updates and just make the identity/divergence
figures::

    $ abfinder -d MyDatabase -s standards.fasta -o /path/to/output/ --no-update

There are several other options, mainly related to formatting the identity/divergence
figures. A complete list of all options can be obtained with::

    $ abfinder --help
