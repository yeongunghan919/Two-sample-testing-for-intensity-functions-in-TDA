# Changelog

## Unreleased

### Fixed

- **Persistence image baseline weighting**
  - **Date discovered:** July 23, 2026
  - The persistence image baseline was implemented using a constant weighting function instead of the root weighting function specified in the manuscript.
  - The implementation was corrected to use the intended root weighting function.
  - All simulation experiments involving the persistence image baseline should be rerun.
  - The corrected PI power values differ only slightly from the originally reported values.
  - The qualitative comparisons, rankings among methods, and conclusions of the manuscript remain unchanged.

### Revision checklist

- [ ] Rerun all simulations involving the persistence image baseline.
- [ ] Verify both type-I error and power results.
- [ ] Replace all affected figures.
- [ ] Replace all affected tables.
- [ ] Update numerical PI results mentioned in the main text or appendix.
- [ ] Confirm that the manuscript description matches the corrected implementation.
- [ ] Update the public reproducibility code.
- [ ] Mention the correction transparently in the revision response letter.

### Suggested response-letter note

> During the review process, we identified an inconsistency in the implementation of the persistence image baseline. Specifically, the simulations used a constant weighting function instead of the root weighting function specified in the manuscript. We corrected the implementation, reran all experiments involving persistence images, and updated the corresponding figures and numerical results. The resulting changes are minor and do not affect the conclusions of the paper.
