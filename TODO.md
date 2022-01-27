# TODO

This document contains remaining TODOs and/or remarks of the existing
implementation. This mostly refers to implementations that do resolve issues,
but might not have been tested extensively in other configurations, or other
systems.

## Pre-inlet communication

The messages send between the pre-inlet and main domain have been significantly
reduced by performing bounding-box checks between the involved ranks. This makes
sure that only atomic blocks that are overlapping perform communication. This
has only been tested with pre-inlets aligned along the positive `x` direction,
i.e. `Direction::Xpos`. The other configurations have *not* been tested yet.

## Solidification communication

The original solidification of platelets near binding sites was fully
sequential. Although the approach works suitable for single core simulations,
multi-core simulations ran into issues (see issue #86), where the binding sites
were not properly propagated into neighbouring atomic blocks.

A patch is added that inserts an explicit synchronisation moment between
labelling and removal of tagged particles (e.g. cells that were close to a
binding site), such that information near the boundaries of atomic blocks is
shared with their neighbours. Although this resolves the original issues for a
simple case, such as described in `cases/solidify_example`, more extensive tests
are required for large(r) scale problems. In preliminary tests, issues were
still observed if regions were marked as binding site that were already
overlapping with cells, creating broken or torn cells in the simulation domain.
