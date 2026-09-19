"""The data-loading interface must fail loudly without SEA-AD access."""

import pytest

from seaad_niches import io


def test_load_merfish_requires_access(tmp_path):
    with pytest.raises(RuntimeError, match="SEA-AD data access"):
        io.load_merfish_cells(cache_dir=tmp_path)


def test_load_snrnaseq_requires_access(tmp_path):
    with pytest.raises(RuntimeError, match="SEA-AD data access"):
        io.load_snrnaseq_reference(cache_dir=tmp_path)


def test_error_message_points_to_docs():
    assert "docs/DATA_ACCESS.md" in io.ACCESS_ERROR
