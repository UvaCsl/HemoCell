# Capillary examples

This example contains two different examples highlighting capillary flows. The
examples consider a wedge-like geometry and a two-capillary bifurcating
configuration.

## Usage

- For the wedge example, set `const auto scenario` to `Scenario::Wedge`, compile
  the problem, and run with the corresponding configuration file:
  `wedge_config.xml`.
- For the two capillary example, set `const auto scenario` to
  `Scenario::Bifurcation`, compile the problem, and run with the corresponding
  configuration file: `bifurcation.xml`.

## Wedge example

This simulation is based on the following article: [Lykov et al. 2017 Probing eukaryotic cell mechanics via mesoscopic simulations.](http://doi.org/10.1371/journal.pcbi.1005726).
In specific, it recreates the geometry as presented in [Figure 3](http://doi.org/10.1371/journal.pcbi.1005726.g003).

- The gaps between the obstacles are 10, 12, and 15um, scaled with the size of
  the WBC diameter (now WBC diameter 8um with a gap of 6um).
- The gap should be larger than the wbc rigid core diameter!!
- The length is 25 um, the larger height is 10 um and the smaller one is 3.4 um.
  The spacing between the rows of obstacles is 60 um Pressure drop: 0.67Pa/μm ->
  in cpp dpdz = 6.7e5

```
       ^    |¯   -   _                                         |¯   -   _
       |    |             ¯    _                               |
       |    |                        |   ^                     |
     10 um  |                        |   | 3.4 um              |
       |    |                        |   |                     |
       |    |             _    ¯         v                     |
       v    |_   -  ¯                                          |_   -  ¯
            ^
            |
          10~15 um                    <---------60um---------->
            |
            v
            |¯   -   _                                         |¯   -   _
            |             ¯    _                               |
            |                        |                         |
            |                        |                         |
            |                        |                         |
            |             _    ¯                               |
            |_   -  ¯                                          |_   -  ¯
            <------------25um------->
```

## Two capillary example

The domain has a single inlet region that splits into two identical capillaries
that combine near the end of the channel. The dividing solid region in the
center has been rounded using ellipsoids (not shown in the ASCII sketch).

```
┌───────────────────────────────────────────────────────────────────────────────┐
│┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼│
│┼┼┼┼┼┼┼┼┼┼────────────────────────────────────────────────────────────┼┼┼┼┼┼┼┼┼│
│┼┼┼┼┼┼┼┼┼│                                                            │┼┼┼┼┼┼┼┼│
│┼┼┼┼┼┼┼┼┼│                                                            │┼┼┼┼┼┼┼┼│
├─────────┘                                                            └────────┤
│                                                                               │
│                                                                               │
│             ┌────────────────────────────────────────────────────┐            │
│             │┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼│            │
│             │┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼│            │
│             │┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼│            │
│             │┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼│            │
│             │┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼│            │
│             └────────────────────────────────────────────────────┘            │
│                                                                               │
│                                                                               │
├────────┐                                                             ┌────────┤
│┼┼┼┼┼┼┼┼│                                                             │┼┼┼┼┼┼┼┼│
│┼┼┼┼┼┼┼┼┼─────────────────────────────────────────────────────────────┼┼┼┼┼┼┼┼┼│
│┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼┼│
└───────────────────────────────────────────────────────────────────────────────┘
```
