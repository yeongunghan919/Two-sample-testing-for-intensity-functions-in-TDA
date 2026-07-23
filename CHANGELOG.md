# Changelog

### Fixed

- **Persistence image baseline weighting**
  - **Date discovered:** July 23, 2026
  - In Torus Simulation, the persistence image baseline was implemented using a constant weighting function instead of the root weighting function specified in the manuscript.
  - The implementation was corrected to use the intended root weighting function.
  - The corrected PI power values differ only slightly from the originally reported values.
  - The qualitative comparisons, rankings among methods, and conclusions of the manuscript remain unchanged.

