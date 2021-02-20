=================
Streaming Pileupy
=================


.. image:: https://img.shields.io/pypi/v/streaming_pileupy.svg
        :target: https://pypi.python.org/pypi/streaming_pileupy

.. image:: https://img.shields.io/travis/winni2k/streaming_pileupy.svg
        :target: https://travis-ci.com/winni2k/streaming_pileupy

.. image:: https://readthedocs.org/projects/streaming-pileupy/badge/?version=latest
        :target: https://streaming-pileupy.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status


.. image:: https://pyup.io/repos/github/winni2k/streaming_pileupy/shield.svg
     :target: https://pyup.io/repos/github/winni2k/streaming_pileupy/
     :alt: Updates



Create multi-sample text-pileups of streaming SAM/BAM files.


* Free software: MIT license
* Documentation: https://streaming-pileupy.readthedocs.io.


Features
--------

Streaming Pileupy creates a pileup of a single SAM/BAM file
using the read group SM identifier to split reads by sample:

.. code-block:: bash

    # extract sample names from read group SM tag
    samtools view -H input.bam \
      | grep '^@RG' \
      | perl -pne 's/.*SM:(\S+).*/$1/' \
      | sort | uniq > sample_names.txt

    # create read-group aware pileup
    spileup input.bam sample_names.txt

Base quality filtering
``````````````````````

Bases with less than a certain quality can be filtered with ``-Q``.


Missing features
----------------

* Read beginning and end annotations in pileup output
* Deletion annotations in pileup output
* Filter output bases on BED file


Speed benchmarks
----------------

Speed benchmarks are available at http://warrenwk.com/streaming_pileupy/

The benchmarks are run using Airspeed Velocity:

.. code-black:: bash

    # install asv into conda environment
    conda create -n spileup_asv asv
    conda activate asv

    # run benchmarks
    FIXTURES=$PWD/benchmarks/fixtures asv run

    # publish results to github pages
    asv gh-page


Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
