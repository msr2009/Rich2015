.. include:: global.rst

Enrich2 Project To Do List
##########################

Debugging
=========

Write test datasets

* :py:class:`~fqread.FQRead`
* trim_fastq.py
* split_fastq.py
* :py:class:`~seqlib.barcode.BarcodeSeqLib`
* :py:class:`~seqlib.basic.BasicSeqLib`
* :py:class:`~seqlib.overlap.OverlapSeqLib`
* :py:class:`~seqlib.barcodevariant.BarcodeVariantSeqLib`
* :py:class:`~selection.Selection`
* :py:class:`~experiment.Experiment`


Implementation - High Priority
==============================

Make a pass for adding additional status messages

Add basic plotting functions

* :py:class:`~seqlib.seqlib.SeqLib`
	* Diversity heatmap

* :py:class:`~selection.Selection`
	* Sequence-function map

* :py:class:`~experiment.Experiment`
	* Sequence-function map

Implement replicate comparison ``[experiment.py]``

Implement control-based correction ``[experiment.py]``

Implement filtering ``[experiment.py]``


Implementation - Low Priority
=============================

Extend logging to use multiple files (specifically to support sending filtered reads to their own output file)

Define a custom logging message formatting - "root" is unnecessary and confusing ``[enrich.py]``

Polish pass through ``import`` statements

Add multithreading to library read stage (one thread per SeqLib to speed up file I/O step)

More advanced management of peak memory usage
	
	* Maximum size of barcode dictionary before being dumped to disk to a (key-sorted) temp file, then merge the counts from the temp files.
	
	* For barcode_map barcodes, perform the barcode_map-based filtering at each threshold, then resolve merging the unmapped barcodes in the temp files.


Documentation
=============

Document logging behaviour (message types output at ``INFO``, ``DEBUG`` levels) and standard log message format::

	Capitalized message [self.name]

Bare ``.html`` landing page for the documentation

Better description of the output file structure (especially *subdirectory* for :py:meth:`~datacontainer.DataContainer.write_data`)

