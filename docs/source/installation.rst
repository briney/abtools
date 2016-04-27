Install
=======

The easiest way to install AbTools locally (on OSX or Linux) is to use pip::

    $ pip install abtools

If you don't have pip, the Anaconda_ Python distribution contains pip along 
with a ton of useful scientific Python packages and is a great way to get 
started with Python.

AbTools does not run natively on Windows, but Windows users can run AbTools with
Docker_ (AbTools is included in the AbStar Docker image)::

    $ docker pull briney/abstar
    $ docker run -it briney/abstar

Stable_ and development_ versions of AbTools can also be downloaded from Github. 
You can manually install the latest development version of AbTools with::

    $ git clone https://github.com/briney/abtools
    $ cd abtools/
    $ python setup.py install

.. note::

    If installing manually via setup.py and you don't already have scikit-bio installed, 
    you may get an error when setuptools attempts to install scikit-bio. This can be fixed 
    by first installing scikit-bio with pip::

        $ pip install scikit-bio

    and then retrying the manual install of AbTools.


Requirements
------------

* Python 2.7.x (Python 3 compatability is in the works)
* biopython_
* celery_
* ete2_
* matplotlib_
* pandas_
* pymongo_
* `scikit bio`_
* seaborn_


Additional dependencies
-----------------------

Several AbTools components have additional non-Python dependencies:

* ``abtools.alignment`` requires MAFFT_ and MUSCLE_
* ``abtools.correct`` requires CDHIT_ and USEARCH_
* ``abtools.mongodb`` requires MongoDB_
* ``abtools.phylogeny`` requires MUSCLE_ and FastTree_
* ``abtools.s3`` requires s3cmd_

If using Docker, all of the the non-Python dependencies are included.


.. _Docker: https://www.docker.com/
.. _Anaconda: https://www.continuum.io/downloads
.. _stable: https://github.com/briney/abstar/releases
.. _development: https://github.com/briney/abstar
.. _biopython: http://biopython.org/
.. _scikit bio: http://scikit-bio.org/
.. _pandas: http://pandas.pydata.org/
.. _pymongo: https://api.mongodb.org/python/current/
.. _celery: http://www.celeryproject.org/
.. _matplotlib: http://matplotlib.org/
.. _ete2: http://etetoolkit.org/
.. _seaborn: https://stanford.edu/~mwaskom/software/seaborn/
.. _MAFFT: http://mafft.cbrc.jp/alignment/software/
.. _MUSCLE: http://www.drive5.com/muscle/
.. _FastTree: http://meta.microbesonline.org/fasttree/
.. _s3cmd: http://s3tools.org/s3cmd
.. _CDHIT: http://weizhongli-lab.org/cd-hit/
.. _USEARCH: http://www.drive5.com/usearch/
.. _MongoDB: https://www.mongodb.org/
