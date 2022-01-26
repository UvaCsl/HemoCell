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
