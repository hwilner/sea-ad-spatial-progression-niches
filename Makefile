# One-command entry points. Real-data stages require the open SEA-AD
# MERFISH download (~500 MB) into data/seaad (not committed).

.PHONY: install test download-data analyze real-pipeline clean

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

# One command: download (if needed) then run the real pipeline.
real-pipeline: download-data analyze

clean:
	rm -rf data/seaad
