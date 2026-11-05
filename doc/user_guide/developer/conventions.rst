Coding conventions
==================

Formatting
----------

C++ code follows ``.clang-format`` in the repository root (clang-format >= 14):
2-space indentation, attached braces, ``Type * name``, unindented namespaces,
120 columns, and includes are never reordered. ``.editorconfig`` carries the
basics (spaces, LF, final newline) to editors without clang-format support.

* All HemoCell sources follow this style. Format what you change before
  committing: ``git clang-format`` (staged changes) or
  ``git clang-format <base commit>``.
* ``scripts/format.sh [--check] [paths]`` formats (or checks) whole files.
* Pull requests are checked by CI on the changed lines only (this keeps the
  check independent of small differences between clang-format versions).
* Commits that only reformat are listed in ``.git-blame-ignore-revs``; enable
  it locally with ``git config blame.ignoreRevsFile .git-blame-ignore-revs``
  (GitHub uses it automatically).
* Vendored and historical code (``external/``, ``legacy/``, the generated
  ``palabos/``) keeps its own formatting. Changes to Palabos belong in
  ``patch/palabos.patch``.

Naming
------

Follow the surrounding code; for new code:

* Types: ``CamelCase`` (``HemoCellParticleField``); functions and variables:
  ``camelCase`` (``applyRepulsionForce``, ``repulsionCutoff``); constants and
  compile-time options: ``UPPER_CASE`` (``FORCE_LIMIT``, ``OUTPUT_FORCE``).
* Everything in ``namespace hemo``, constants included (``hemo::OUTPUT_FORCE``,
  ``hemo::RBC_FROM_SPHERE``); macros only for compile-time options.
* No ``using namespace`` in headers (``.h``/``.hh``): qualify ``std::`` and
  ``plb::`` names there. Source files of the library may use
  ``using namespace std;`` and ``using namespace plb;``; case files declare
  what they import themselves (usually ``hemo``, ``plb`` and ``std``).
* No process-wide state that outlives a simulation: anything static is reset
  in ``resetGlobalState()`` (``core/hemoCell.cpp``), which ``~HemoCell`` calls;
  per-simulation data belongs to an object (as in ``HemoCellStretch``).
* Existing misspellings in the public API (e.g. ``...TimeScaleSeperation``) are
  kept for compatibility; add a correctly spelled alias rather than renaming.

Units
-----

Internally everything is in lattice units (``dx``, ``dt``, ``dm`` from
``hemo::Parameters``, also available as ``param::``). Interfaces mix SI values,
micrometres, iterations and lattice units, so every public function and XML
setting documents its unit, and conversion happens once, where a value enters
(e.g. ``setRepulsion()`` converts ``k`` from pN and the cutoff from
micrometres). Never retune a
material constant, force cap or resolution to hide a regression.

Performance rules
-----------------

HemoCell's cost is dominated by per-particle and per-node loops. In those loops:

* no heap allocation, no ``std::map`` lookups, no XML reads, no logging;
* no validity checks that can be done once: validate inputs at setup (see
  ``HemoCell::sanityCheck()``, run before the first step) and keep the loops
  free of branches that never trigger;
* prefer indices into contiguous storage over pointers, and keep data that is
  used together close together;
* hoist invariants (e.g. ``param::`` values, block sizes) out of loops.

Measure with the built-in profiler (``hemo::global.statistics``) on the same
case, ranks and workload before and after a change, and check that the results
are unchanged (bit-identical where the change should not alter numerics).

Parallel correctness
--------------------

* Collective operations must be called by all ranks in the same order; a
  function that gathers data or writes coordinated output is collective even
  if only rank 0 writes.
* Point-to-point exchanges use a private communicator and receive from known
  source ranks (see :ref:`developer/architecture:Communication`).
* Runs are bit-for-bit reproducible for a given rank count and decomposition.
  Keep it that way: when the order of processing received data affects the
  order of particles or of floating-point sums, process it in a fixed order
  (e.g. by source rank), never in MPI arrival order. Reproducibility is what
  makes exact before/after comparisons possible; checkpoints of identical runs
  are byte-identical (``cmp`` their ``.dat`` files). Records written as raw
  bytes need their padding zeroed for that (see ``HemoCellParticle``).
* Output and statistics must not change the simulation state: read complete
  cells inside a ``HemoCellFields::CompleteCellsScope`` (a copy of the
  particles, discarded afterwards) instead of calling ``syncEnvelopes()`` or
  ``deleteIncompleteCells()``, and leave deleting cells to ``iterate()``.
* Cell IDs are ``plint`` and grow with every periodic crossing; never store
  them in an ``int``.
* A single-rank run does not test communication. Changes to particle transfer,
  ownership, periodicity or decomposition need runs on several ranks, including
  block and periodic crossings, and a restart from a checkpoint.

Errors and logging
------------------

* Use ``hlog`` (stdout and log file, rank 0) and ``hlogfile`` (log file only)
  for messages, ``pcout`` for console-only output; never log per particle.
* Report errors with ``hemo::fatal(...)`` (``config/error.h``), never
  ``exit()``: it prints the message with the MPI rank on the rank that calls it
  (``hlog`` prints nothing on the other ranks) and ends the job with
  ``MPI_Abort``. ``hemo::warning(...)`` reports from any rank. Both stream
  their arguments into the message. Under a ``FatalErrorsThrow`` guard,
  ``fatal()`` throws ``FatalError``, so tests can check error paths (only where
  every rank reaches the error).
* Configuration values are read with ``read<T>()``, which throws
  ``std::invalid_argument`` for missing, empty or invalid values. Optional
  settings catch it and use a default; required ones let it propagate.
* Do not promise validation that the code does not perform.

Documentation
-------------

* Public functions and non-trivial numerical helpers get Doxygen comments
  (``@brief``, ``@param[in/out]``, ``@return``) that state units, coordinate
  frames, array layouts, ownership and whether all ranks must call them.
* Behaviour changes update the user guide (``doc/user_guide``) and
  ``CHANGELOG.md`` in the same change.

Building and testing
--------------------

See ``AGENTS.md`` (build, presets, test commands) and
:ref:`developer/extending:Adding tests`. Before opening a pull request, run
the unit tests (``ctest -L unit``), the validation tests when touching physics,
mechanics or communication (``ctest -L validation``), and the unit tests with
the ``sanitize`` preset.
