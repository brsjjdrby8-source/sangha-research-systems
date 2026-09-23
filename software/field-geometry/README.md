# SANGHA Field Geometry v0.1.0

A runnable reference implementation for object-field geometry, distribution moments,
and pointwise confidence intervals. SANGHA Research Systems.

Typed observations → exhaustive geometry → explicit estimand → cluster uncertainty
→ reviewable artifacts → provenance receipt.

This is a new reference kernel, not a promoted copy of Cherie's QuPath kernel.
Analytic and independent directional-reference tests establish its current gate.
No legacy QuPath parity, biological validation, 3D reconstruction, or calibrated
finite-sample coverage is claimed.

## Start

Python 3.11+ and GNU Make are required. From the extracted project directory:

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -e .
make doctor
make test
make demo
```

Open `runs/demo/report.html`. The demo contains synthetic polygons only.
This source directory generates the demo afresh; the downloadable release ZIP also includes completed outputs.
The delivered archive includes a completed example in `runs/demo`, so to run a
fresh example use `make demo OUT=runs/my-demo` instead. Existing output directories
are refused rather than overwritten.

For TIFF labels: `python -m pip install tifffile`. NumPy `.npy` labels require no
additional package. Unsegmented intensity images are not accepted as objects.

## Make grammar

```text
command := make verb assignment*
verb := help | doctor | validate | measure | estimate | run | verify
      | replay | import-labels | test | demo | package
assignment := DATA=path | CONFIG=path | OUT=path | MEASURE_OUT=path
            | REPLAY_OUT=path | LABELS=path | IMPORTED=path | PYTHON=executable
```

| Verb | Input | Output / contract |
|---|---|---|
| `doctor` | Python runtime | Versions and source hash |
| `validate` | DATA + CONFIG | Strict schema keys, calibration, identity, boundary checks |
| `measure` | DATA + CONFIG | Exhaustive geometry and angular samples in MEASURE_OUT |
| `estimate` / `run` | DATA + CONFIG | Full results, CSVs, HTML, gates and receipt in OUT |
| `verify` | OUT | Byte checks against its unsigned receipt |
| `replay` | Saved OUT inputs | Fresh REPLAY_OUT and scientific artifact parity |
| `import-labels` | LABELS manifest | Canonical IMPORTED dataset and import receipt |
| `test` | Source + fixtures | Analytic/reference/inference/integrity gates |
| `demo` | Bundled synthetic dataset | Complete analysis and verification |
| `package` | Source + passing tests | Versioned deterministic-entry ZIP and SHA-256 |

```sh
make validate DATA=my.dataset.json CONFIG=examples/default.config.json
make measure DATA=my.dataset.json MEASURE_OUT=runs/case-measure
make run DATA=my.dataset.json OUT=runs/case
make verify OUT=runs/case
make replay OUT=runs/case
make import-labels LABELS=examples/labels.manifest.json IMPORTED=imported.dataset.json
make package
```

`measure` and `estimate` share the kernel. Estimation recomputes from the input,
rather than trusting an intermediate CSV. `replay` is a deterministic rerun in the
current environment, not a portable environment reconstruction. Receipt hashes
can differ if source/environment or the original JSON formatting differ; scientific
output bytes are compared explicitly. No cleanup target deletes user runs.

## Observation contract

`examples/demo.dataset.json` is the full input example. Each field has:

- `specimen_id`, `field_id`: nonempty strings, jointly unique.
- `area`: sampled reference-frame area in physical `unit` squared.
- `pixel_size_x`, `pixel_size_y`: physical unit per pixel, positive and finite.
- `sampling`: `complete_objects_by_reference_point`.

Each object has those field IDs, a unique `object_id`, pixel-coordinate `points`
(an Nx2 boundary array), and `complete: true`. Coordinates need not be ordered.
Objects must be selected by a reference point in a sampling frame, with enough
surrounding image to recover the complete boundary. User declarations cannot be
verified from polygon coordinates alone. Convenience-selected fields do not become
representative merely by applying a bootstrap.

For labeled images, `examples/labels.manifest.json` supplies calibration and a
half-open frame `[x0,y0,x1,y1]` in pixel-edge coordinates. Each nonzero integer label
is treated as one object (the adapter does not relabel disconnected components).
Selection is by pixel-cell centroid. A selected label touching the acquired image
edge is refused: acquire a larger guard region. Labels are interpreted as unions
of pixel cells; their convex hull is retained. Background is zero. Physical field
area is frame pixel area times x/y calibration. All fields, including empty fields,
are retained. Fully empty datasets cannot have an object-size distribution.

Do not include the same object in overlapping sampled frames. This version cannot
detect biological duplication across images. Import receipts record source image
hashes; analysis receipts retain the converted dataset. Keep the import receipt
with the original images to preserve the entire provenance chain.

## Geometry contract

EXACT means exhaustive supplied polygon-hull geometry evaluated in floating point:

- Maximum Feret: maximum pairwise hull-vertex distance.
- Minimum Feret: minimum support width over hull-edge normals.
- Mean Feret: convex-hull perimeter / pi, analytically averaged over directions.
- Feret ratio: max / min.
- Hull area and hull perimeter: explicitly hull properties, not concave object
  area/perimeter. Feret widths depend only on the convex hull.

The finite angular grid is uniform on `[0, pi)`; its distribution moments are
**quadrature estimates**, not exact integrals. Increase `angles` to check numerical
stability. Raster discretization and segmentation are separate error sources.
The maximum-width implementation has O(h^2) pairwise storage for h hull vertices;
this reference version is intended for bounded object outlines, not huge meshes.

## Statistical contract

All moments are empirical plug-in estimates. Variance uses the probability-weighted
second central moment, not a Bessel-corrected sample variance. Skewness and kurtosis
are correspondingly uncorrected. Configurable raw moments 1..K and central moments
2..K support K=4..12; the default is 4. Order 0 is 1 and central order 1 is 0, so they
are omitted. Higher moments can overflow or be dominated by rare large objects.

For each scalar geometry feature: mean, variance, SD, CV, skewness, excess kurtosis,
raw/central moments, 5/25/50/75/95th percentiles, and IQR. Quantiles use the inverse
weighted empirical CDF, without interpolation. Dimensionless CV/shape summaries are
undefined at zero variance or mean where applicable. Undefined values are null
with explicit status, not silently replaced by zero.

The angular profile is a separate joint object × direction distribution, not the
distribution of mean object width. Directions remain attached to objects, with
within-object and between-object directional variance reported separately.

Weighting options:

- `specimen_equal` (default): choose a specimen uniformly, then an observed object
  uniformly within that specimen. An entirely empty specimen makes this object
  estimand undefined and is refused. Density/burden are means of specimen ratios.
- `object_pooled`: pool observed objects with equal weights. Density/burden use
  pooled count or power-sum divided by pooled area. Empty specimens remain in the
  density denominator. Unequal acquisition effort changes this target.

Density is objects / physical area. Burden is the sum of powers of **object mean
Feret** / physical area; it is not the power moment of the angular distribution.
Neither weighting corrects section-selection bias to a 3D particle population.

The default CI resamples independent specimens with replacement, retaining every
field, object and direction together. Repeated specimen draws retain multiplicity.
It does not automatically resample nested fields/objects, which can double-count
variation for some designs. Fixed specimens, single-specimen field surveys,
stratification, paired conditions, and spatial-block designs need different
resampling contracts and are not implemented here.

Each finite statistic gets a pointwise percentile interval at `confidence` (default
0.95). These are approximate bootstrap CIs, not simultaneous coverage of the whole
profile. With one specimen the point estimate remains available, but CIs are null.
If any replicate is undefined, the corresponding interval is withheld rather than
silently conditioning on successful replicates. Fewer than 10 specimens, degenerate
intervals, and sparse bootstrap tails are flagged. This is not a coverage guarantee.

`fields.csv` is descriptive per-field data; it does not pretend object replication
provides a biological CI for a single field. Population CIs are in `statistics.csv`.
Sampling intervals exclude segmentation, pixel calibration and angular quadrature
uncertainty. For higher moments, adequate population moments and regularity are
assumptions, not verified properties of the sample.

## Outputs and units

- `objects.csv`: calibrated exhaustive object measurements.
- `directional_widths.npy`: object × angle widths, matching objects.csv row order.
- `fields.csv`: counts, areas, mean-Feret moments, and burden (including empty fields).
- `statistics.csv`, `results.json`: population estimates, pointwise CIs and statuses.
- `report.html`: local readable result surface, with downloadable tables.
- `dataset.json`, `config.json`: frozen inputs for replay.
- `gates.json`: admission and inference status.
- `receipt.json`: source, dependency versions, input and output SHA-256 values.

Receipts are unsigned integrity records, not security attestations. `verify` checks
listed file bytes; it does not prove scientific validity or trustworthiness of a
modified receipt.

A length feature's raw/central kth moment has unit^k; its variance has unit^2;
its mean/SD/quantiles/IQR retain unit. Hull area starts at unit^2. Feret ratio,
CV, skewness and excess kurtosis are dimensionless. Density is unit^-2;
mean-Feret power-k burden is unit^(k-2). CSV includes the unit expression.

## Validation and extensions

Tests cover rectangles, circle-polygon convergence, anisotropic calibration,
rotation/translation invariance, independent directional reference parity, known
moments, variance decomposition, weighting/multiplicity, input-order invariance,
single-specimen refusal of CIs, invalid inputs, label pixel-cell geometry, border
refusal, reproducible outputs, and receipt tamper detection.

Planned adapters: QuPath/Fiji tables and contours, spatial-block or design-specific
resampling, simultaneous intervals, comparative contrasts, uncertainty propagation,
model-specific 3D inverse estimates, and adaptive field acquisition. These must
preserve the reference kernel and pass numerical parity before promotion. None is
claimed as implemented in v0.1. The runtime contract is in `operator-contract.v1.json`. Exact tested dependency
versions are in `requirements-tested.txt`; broader supported ranges are in
`pyproject.toml`. The Make/CLI contract is ready for a later GUI or
Nextflow adapter without duplicating the numerical operators.

References:
- Cauchy projection formula: https://arxiv.org/abs/1604.05815
- Wicksell model assumptions: https://pmc.ncbi.nlm.nih.gov/articles/PMC4890128/
- SciPy ConvexHull: https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.ConvexHull.html
- Percentile bootstrap semantics: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.bootstrap.html
