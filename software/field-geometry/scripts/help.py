print('''SANGHA FIELD GEOMETRY — Make grammar

make <verb> [DATA=dataset.json] [CONFIG=config.json] [OUT=runs/case]

  doctor        Verify runtime imports and report versions/source digest
  validate      Validate input, calibration, object identity and geometry
  measure       Export exhaustive geometry to MEASURE_OUT (default OUT-measure)
  estimate      Compute moment profile and specimen-cluster intervals
  run           Alias for estimate
  verify        Check output hashes against the unsigned receipt
  replay        Recompute saved inputs in REPLAY_OUT and compare scientific outputs
  import-labels Convert LABELS manifest to IMPORTED dataset
  test          Run analytic, reference parity, inference and receipt tests
  demo          Run bundled synthetic example and verify outputs
  package       Require tests, then build versioned ZIP + SHA-256

Examples:
  make demo
  make run DATA=my.dataset.json CONFIG=my.config.json OUT=runs/case01
  make replay OUT=runs/case01
  make measure DATA=my.dataset.json MEASURE_OUT=runs/measurement01

Outputs must be new directories. No implicit overwrite or destructive clean target.
EXACT means exhaustive supplied polygon geometry, not exact biological boundaries
or exact angular integration. Angular moments are explicitly quadrature-based.''')
