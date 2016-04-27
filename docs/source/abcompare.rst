AbCompare
=========


Overview
--------

AbCompare is used to perform repertoire-level comparison of antibody
sequence data using a variety of similarity and divergence measures. 
Currently, AbCompare compares samples using the frequency of V-gene
use, although other comparison types (such as clonality) are planned.
However, the underlying similarity and divergence functions are accessible
via the AbTools API, so you can compare samples using other characteristics.

Similarity (or divergence) scores are computed by subsampling each dataset 
and computing the score for the subsamples. This process is repeated many
times, and the median score for all of the iterations is returned. In addition
to producing a more accurate representation of the true score, it also makes
is possible to directly compare datasets of different sizes.


Examples
--------

To compute the Marisita-Horn similarity of two collections, both in
the same MongoDB database::

    $ abcompare -d MyDatabase -1 Collection1 -2 Collection2 -o /path/to/output/

If only one collection is provided (via ``-1``), then that collection will
be iteratively compared to all other collections in the database::

    $ abcompare -d MyDatabase -1 Collection1 -o /path/to/output/

If you leave out collections entirely, all collections in the database will be
iteratively compared to all other collections::

    $ abcompare -d MyDatabase -o /path/to/output/

Finally, if you'd like to compare only those collections that share a common
prefix (for example, if your collection names are formatted as ``SubjectName_Timepoint``
and you'd like to compare all the timepoints from a single subject)::

    $ abcompare -d MyDatabase --collection-prefix SubjectName -o /path/to/output/

The default comparison method is Marisita-Horn similiarity, but several other
methods are provided:

    - `Marisita-Horn similarity`_ (``'marisita-horn'``)
    - `Kullback-Leibler divergence`_ (``'kullback-leibler'``)
    - `Jensen-Shannon similarity`_ (``'jensen-shannon'``)
    - `Jaccard similarity`_ (``'jaccard'``)
    - `Bray-Curtis similarity`_ (``'bray-curtis'``)
    - `Renkonen similarity`_ (``'renkonen'``)
    - `Cosine similarity`_ (``'cosine'``)


.. _Marisita-Horn similarity: https://en.wikipedia.org/wiki/Morisita%27s_overlap_index
.. _Kullback-Leibler divergence: https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence
.. _Jensen-Shannon similarity: https://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence
.. _Jaccard similarity: https://en.wikipedia.org/wiki/Jaccard_index
.. _Bray-Curtis similarity: https://en.wikipedia.org/wiki/Bray%E2%80%93Curtis_dissimilarity
.. _Renkonen similarity: https://en.wikipedia.org/wiki/Renkonen_similarity_index
.. _Cosine similarity: https://en.wikipedia.org/wiki/Cosine_similarity

To use an alternate comparison method, pass the method with the ``--method`` option::

    $ abcompare -d MyDatabase -o /path/to/output/ --method jaccard

The number of sequences used in each iteration (``--chunksize``, default is 100,000) and the 
number of iterations (``--iterations``, default is 10,000) can also be changed::

    $ abcompare -d MyDatabase -o /path/to/output --iterations 1000 --chunksize 25000

As with other AbTools applications, there are options for connecting to remote MongoDB
servers (``--ip`` and ``--port``) and MongoDB authentication (``--user`` and ``--password``).
A complete list of AbCompare options can be obtained by::

    $ abcompare --help