# One-command entry points. Real-data stages require the open SEA-AD
# MERFISH download (~500 MB) into data/seaad (not committed).

.PHONY: install test download-data analyze freeze-dataset real-pipeline clean

install:
	pip install -e ".[dev,spatial]"

test:
	python -m pytest -q

# Download + verify the public SEA-AD MTG MERFISH h5ad (no DUA required).
download-data:
	python scripts/download_seaad.py --cache-dir data/seaad

# Neighbor enrichment + niche detection on the staged real cells;
# writes small CSV/JSON outputs to reports/.
analyze:
	python scripts/analyze_seaad.py --cache-dir data/seaad

# Freeze donor covariates, QC report, and the hashed Paper 1 manifest.
freeze-dataset:
	python scripts/freeze_paper1_dataset.py --cache-dir data/seaad

# One command: download (if needed), run the real pipeline, and freeze metadata.
real-pipeline: download-data analyze freeze-dataset

clean:
	rm -rf data/seaad
