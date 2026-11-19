# edgeicg

Spillover-aware interaction-changed genes at three spatial scales: touching cells (contact), short-range
neighbours (paracrine) and the regional neighbourhood (niche).

For every receiver type, sender type and gene, `edgeicg()` tests whether the receiver's expression changes
with its exposure to the sender. Transcripts leaked from touching neighbours are modelled from the data and
removed per cell before testing.

At niche scale, co-located senders are collinear, so the analysis runs in two stages:
1. A marginal screen controls the FDR. It is calibrated by a toroidal-shift null.
2. Exact fine-mapping over 1- and 2-sender models attributes each hit to a sender, with a credible set.

## Scale
- Counts are read only in gene blocks, through a provider: `counts_provider()` for an in-memory matrix, or `counts_provider_store()` for a GiottoDisk store.
- Memory is bounded by `block_genes`.
- Blocks run sequentially or on `mirai` daemons (`workers = `), with identical results.

## Use
```r
library(edgeicg)
sim <- simulate_edgeicg()
fit <- edgeicg(sim$counts, sim$coords, sim$types, lib = sim$lib, min_shift = 150)
fit
hits <- edgeicg_classify(fit, edgeicg_results(fit))
```

For real data:
- Pass touching-cell edges from the cell outlines with `contact_from_polygons()` as `spill_edges`.
- Pass the library size over all genes as `lib`.

## Status
0.1.0 is a port of the 2026 ICG2 prototype. With `drop_degenerate = FALSE`, it reproduces the prototype exactly on the kidney and breast crops.
