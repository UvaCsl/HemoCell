Saving only the CSV output
==========================

This is possible with the extra ``<sim><tcsv>`` parameter which is already added
in the :ref:`cases/pipeflow:Pipe flow` case. This parameter controls a
separate call to ``writeCellInfo_CSV`` that only writes CSV output. Thus, by
increasing ``<sim><tmeas>`` the interval of the HDF5 *and* CSV can be increases,
where the CSV output is then separately set by ``<sim><tcsv>``.

.. note::

  Don't forget to include the right header:

  .. code-block:: c++

     #include "writeCellInfoCSV.h"
