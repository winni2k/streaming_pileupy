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

.. code-block:: text

    spileup input.bam sample_names.txt


Missing features
----------------

* Base quality filtering
* Read beginning and end annotations in pileup output
* Deletion annotations in pileup output
* Filter output bases on BED file

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
