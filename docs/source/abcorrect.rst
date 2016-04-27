AbCorrect
=========


Overview
--------

AbCorrect is a full-featured utility for performing error-correction
on antibody repertoire sequencing data. Error correction can be performed
using Unique Antibody IDs (UAIDs; also known as molecular barcodes) or
by identity clustering.

The AbCorrect command-line application is designed to work with antibody
sequence data that has already been annotated with AbStar_. Although other
error correction tools for antibody repertoire sequencing operate on 
raw data (or, in the case
of paired reads, on raw data after read merging), we have found
that annotating the sequences with AbStar before error correction allows us to focus clustering
and consensus/centroid generation on just the VDJ region of the antibody sequence
in the proper orientation. In our hands, this tends to produce more accurate,
reproducible results.

In addition to being provided as a command-line application, the core
functionality can be accessed through the :ref:`AbTools API <API-correct>`,
which allows AbCorrect to be integrated into more sophisticated sequence
processing pipelines. Below are several examples showing how to use AbCorrect 
as a command-line application.


.. _AbStar: https://github.com/briney/abstar/


Examples
--------

The simplest use case for AbCorrect is to perform error correction on a single
JSON file, which is the output from running AbStar on a single FASTA/Q file::

    $ abcorrect -j /path/to/MyData.json -t /path/to/temp/ -o /path/to/output/

This will cluster sequences based on UAIDs (which have been pre-parsed by AbStar)
and generate a 'consensus' sequence for each UAID cluster containing at least one
sequence ('consensus' is in quotes because it isn't truly a consensus for UAID bins
with a single sequence). Output will be a single FASTA file of error-corrected sequences,
located at ``/path/to/output/MyData.fasta``.

To perform the same operation, but only calculate consensus sequences for UAID
clusters with at least 3 sequences::

    $ abcorrect -j /path/to/MyData.json -t /path/to/temp/ -o /path/to/output/ --min 3

If you want to correct errors using UAIDs but you forgot to have AbStar parse them,
you can have AbCorrect parse them by passing the length of the barcode (in 
nucleotides)::

    $ abcorrect -j /path/to/MyData.json -t /path/to/temp/ -o /path/to/output/ --parse-uaids 20 --min 3

This will use the first 20nt of the raw merged read as the UAID. If the UAID is at the end
of the read (for paired reads, this would be the start of R2), use a negative number for
``--parse-uaids``.

To cover the relatively rare case (assuming the UAID length was selected appropriately) where two
sequences were tagged with the same barcode, AbCorrect clusters the sequences within each
UAID bin and builds a consensus/centroid sequence for each subcluster that passes the
``--min`` size threshold. To disable this, you can pass the ``--largest-cluster-only`` option
and AbCorrect will only build a consensus/centroid sequence for the largest cluster within
each UAID bin.

To perform error-correction using identity clustering instead of UAIDs, you can::

    $ abcorrect -j /path/to/MyData.json -t /path/to/temp/ -o /path/to/output/ --no-uaids

This will cluster the sequences at an identity threshold (default is 0.975, or 97.5% identity)
and build a consensus sequemce for each cluster. To cluster with a threshold of 0.96 instead::

    $ abcorrect -j /path/to/MyData.json -t /path/to/temp/ -o /path/to/output/ -I 0.96 --no-uaids

If you have more than one JSON file to be processed, you can pass AbCorrect a directory
that contains one or more JSON files and each JSON file will be iteratively processed::

    $ abcorrect -j /path/to/JSONs/ -t /path/to/temp/ -o /path/to/output/

All of the other options (such as the minimum number of sequences for consensus/centroid 
calculation) remain, although there is currently no way to specifiy different options 
for each JSON file. 


If your AbStar-annotated sequences have already been uploaded to MongoDB, you can still
use AbCorrect to perform error correction. Rather than passing JSON files with ``-j``, you can
pass a MongoDB database name with ``-d`` and a collection name with ``-c``::

    $ abcorrect -d MyDatabase -c MyCollection -t /path/to/temp/ -o /path/to/output/

If you supply just the database name (without a collection), AbCorrect will iteratively process
all collections in the supplied database::

    $ abcorrect -d MyDatabase -t /path/to/temp/ -o /path/to/output/

The above example is querying MyDatabase on your local instance of MongoDB. To do the same 
thing on a remote MongoDB server, you can pass the IP address with ``-i`` (assuming the 
default port of ``27017``::

    $ abcorrect -d MyDatabase -i 123.45.67.89 -t /path/to/temp/ -o /path/to/output/

If your MongoDB server uses a port other than ``27017``, you can provide it using the ``--port`` 
option. And if your remote MongoDB server requires authentication, you can supply the username with
``--user`` and the password with ``--password``. If you don't supply both ``--user`` and 
``--password``, AbCorrect will attempt to connect to the MongoDB database without authentication.

Finally, to make non-redundant set of sequences, AbCorrect provides the ``--nr`` option::

    $ abcorrect -d MyDatabase -t /path/to/temp/ -o /path/to/output/ --nr

This uses ``sort | uniq``, which is much faster than clustering at 100% identity with CD-HIT. 

.. warning::
 
    Using ``--nr`` is is not the same as clustering at 100% identity. Two sequences that are 
    different lengths but are otherwise identical will be collapsed when clustering with CD-HIT 
    but will not be collapsed when using ``sort | uniq``.

