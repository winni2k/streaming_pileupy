=====
Usage
=====

Streaming Pileupy is meant to be run through its command
line interface ``spileup``. Given a sam file ``in.sam``:

.. code-block:: text

    @HD VN:1.6  SO:unknown
    @SQ SN:chr1 LN:1000000
    @RG ID:0    SM:sample_0
    @RG ID:1    SM:sample_1
    r0  0   chr1    24  0   1M  *   0   0   G   I   RG:Z:0
    r1  0   chr1    24  0   1M  *   0   0   G   I   RG:Z:1

, and a sample name file ``sample_names.txt``:

.. code-block:: text

    sample_0
    sample_1

To run a pileup of the reads in ``in.sam`` for
the samples in ``sample_names.txt``, run the command:

.. code-block:: bash

    spileup input.sam sample_names.txt
