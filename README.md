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
Every statistic of the model is a sum over cells, so each pass is a sweep of tasks that return partial sums.

- **Counts** are read through a provider:
  - `counts_provider()` for an in-memory matrix;
  - `counts_provider_store()` (or a GiottoDisk expression store passed directly) for data on disk.
- **In memory** (the default), a task is a gene block over all cells. Memory is bounded by `block_genes`.
- **Above `window_cells` cells** (the default for a store: one window per worker), cells are visited in spatially ordered windows and the sums are added up.
  - Memory then no longer grows with the section.
  - A store is baked once in window order with `GiottoDisk::storeWrite()` and read in contiguous cell ranges, once per pass.
- **Tasks** run sequentially or on `mirai` daemons (`workers = `). Results are identical either way, and identical with and without windows (up to rounding).
- **Rust kernels:** with the [edgeicgrs](https://github.com/drieslab/edgeicgrs) package installed, `backend = "auto"` computes the matrix products in Rust (faer). That is 20–25× faster than R's reference BLAS for these shapes, on every platform.

On the whole Atera kidney section (399k cells, 3,000 genes):
- 5 min in memory with 8 daemons;
- 6 min from a GiottoDisk store with 4 daemons (41 GB peak);
- 10 min with gene blocks from the store (31 GB).

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

## Effect sizes and confounding
- **Effect size per test:** `lfc`, the log2 fold change between receivers at the 10th and 90th percentile of their exposure. By default, results need |lfc| ≥ 0.25 (`edgeicg_results(min_lfc = 0.25)`).
- **Why a floor:** at whole-section power, significance alone also picks up small associations caused by zonation or by sub-states within a cell type. On the Atera kidney section the floor cut null-gene calls from 27% to 6% without losing planted effects. It is free: no refit.
- **Spatial covariates (optional):** `covariates = spatial_basis(coords)`, usually with `covariate_scope = "niche"`. They adjust for smooth spatial trends and fix the niche-scale confounding (pseudo-null 4.7% → 0.6%) at about 1.4× the run time.

## Giotto
`Giotto::findInteractionChangedFeats(g, cluster_column = ..., method = "edgeicg")` runs this engine on a giotto object. The expression can be in memory or GiottoDisk-backed.

EdgeICG succeeds the original interaction-changed feature (ICF) test of Giotto, which compared receivers with and without a neighbour of the sender type. It is the default method of `findInteractionChangedFeats()`; the original test stays available as `method = "classic"`.

## Status
0.1.0. Results are identical for both backends (R and Rust), sequential or parallel, and for one window or many (up to rounding).
