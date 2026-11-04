.. _contributing:

Contributing
============

HemoCell is open-source and community contributions are very welcome, whether
they are bug reports, documentation improvements, new mechanical models, or new
example cases.

Before you contribute
---------------------

* Please make sure that your contribution falls under the HemoCell license (see
  the ``LICENSE`` file in the repository root).
* If you want to resolve a bug, first make sure that it still exists. You can
  build the latest ``master`` branch and verify that the error is reproducible.
* Make sure that the bug you want to report is not already reported in our
  `GitHub issues <https://github.com/UvaCsl/HemoCell/issues>`_ and that no one
  is working on it.
* If you have any questions about the software or if you are facing any issues
  using it, feel free to open a `GitHub discussion
  <https://github.com/UvaCsl/HemoCell/discussions>`_.

Contributing code
-----------------

#. Create a fork of the HemoCell repository.
#. Create a new branch from the ``dev`` branch for the issue you want to
   work on. Please give the branch a name that is relevant to the issue.
#. Modify or add code on your branch, following the
   :ref:`coding conventions <developer/conventions:Coding conventions>`, and
   add tests for new behaviour (:ref:`developer/extending:Adding tests`).
#. Format your changes (``git clang-format``), and run the unit tests
   (``ctest -L unit``) and, for changes to physics, mechanics or
   communication, the validation tests (``ctest -L validation``).
#. Make sure you have not included unrelated changes, and describe
   behaviour changes in ``CHANGELOG.md`` and the user guide.
#. Push your changes and open a pull request from your branch to the ``dev``
   branch. Continuous integration builds HemoCell with GCC and Clang, runs the
   tests on one and two MPI ranks and under sanitizers, and checks the
   formatting of the changed lines.

.. tip::

   If you are adding a new mechanical model or feature, consider also adding a
   short example case and a documentation page so that other users can benefit
   from your work. See :ref:`contributing-docs` below.

.. _contributing-docs:

Contributing to the documentation
----------------------------------

This documentation is built with `Sphinx <https://www.sphinx-doc.org>`_ from the
``doc/user_guide`` directory and with `Doxygen <https://www.doxygen.nl>`_ for
the C++ API reference. Pages are written in either reStructuredText (``.rst``)
or Markdown (``.md`` via MyST).

To build the documentation locally:

.. code-block:: bash

   # from the repository root
   cd doc
   python3 -m venv .venv
   source .venv/bin/activate
   pip install -r requirements.txt

   # build the HTML site (Sphinx invokes Doxygen automatically)
   cd user_guide
   make clean
   make html

The generated site is placed in ``doc/user_guide/_build/html``. Open
``_build/html/index.html`` in a browser to preview your changes.

.. note::

   Doxygen and Graphviz must be installed as system packages (they are not
   Python modules). On Debian/Ubuntu: ``sudo apt-get install doxygen graphviz``.
