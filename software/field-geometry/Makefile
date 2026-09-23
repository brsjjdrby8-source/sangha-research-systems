# SANGHA operator grammar. Override values at invocation, never edit recipes.
SHELL := /bin/bash
.SHELLFLAGS := -eu -o pipefail -c
.DEFAULT_GOAL := help
PYTHON ?= python3
DATA ?= examples/demo.dataset.json
CONFIG ?= examples/default.config.json
OUT ?= runs/demo
MEASURE_OUT ?= $(OUT)-measure
REPLAY_OUT ?= $(OUT)-replay
LABELS ?= examples/labels.manifest.json
IMPORTED ?= runs/imported.dataset.json
export PYTHONPATH := $(CURDIR)/src$(if $(PYTHONPATH),:$(PYTHONPATH))
CLI = $(PYTHON) -m sangha_geometry.cli

.PHONY: help doctor validate measure estimate run verify replay test demo import-labels package
help:
	@$(PYTHON) scripts/help.py

doctor:
	$(CLI) doctor

validate:
	$(CLI) validate "$(DATA)" --config "$(CONFIG)"

measure: validate
	$(CLI) measure "$(DATA)" --config "$(CONFIG)" --out "$(MEASURE_OUT)"

# Estimate always recomputes the exhaustive kernel from the hashed source.
# It never trusts an unverified intermediate measurement table.
estimate: validate
	$(CLI) run "$(DATA)" --config "$(CONFIG)" --out "$(OUT)"

run: estimate

verify:
	$(CLI) verify "$(OUT)"

replay: verify
	$(CLI) run "$(OUT)/dataset.json" --config "$(OUT)/config.json" --out "$(REPLAY_OUT)"
	$(PYTHON) scripts/compare_runs.py "$(OUT)" "$(REPLAY_OUT)"

test:
	$(PYTHON) -m unittest discover -s tests -v

demo:
	$(MAKE) run DATA=examples/demo.dataset.json CONFIG=examples/default.config.json OUT="$(OUT)"
	$(MAKE) verify OUT="$(OUT)"

import-labels:
	$(CLI) from-labels "$(LABELS)" --out "$(IMPORTED)"

package: test
	$(PYTHON) scripts/package.py
