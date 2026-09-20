# SEA-AD data sources (exact versions and access paths)

All sources below are **open access** (AWS Open Data Registry /
CELLxGENE). No controlled-access material is used by this repository.

## MERFISH (spatial transcriptomics, middle temporal gyrus)

- **Dataset:** SEA-AD MTG MERFISH, combined all-donors processed h5ad.
- **File:** `middle-temporal-gyrus/all_donors-h5ad/SEAAD_MTG_MERFISH.2024-12-11.h5ad`
  (release date 2024-12-11; 497,737,587 bytes).
- **Access path:** public HTTPS, no credentials:
  `https://sea-ad-spatial-transcriptomics.s3.amazonaws.com/middle-temporal-gyrus/all_donors-h5ad/SEAAD_MTG_MERFISH.2024-12-11.h5ad`
- **Registry:** AWS Open Data Registry, "Seattle Alzheimer's Disease Brain
  Cell Atlas (SEA-AD)": https://registry.opendata.aws/allen-sea-ad-atlas
- **Citation / DOI:** Gabitto, Travaglini et al. (2024), "Integrated
  multimodal cell atlas of Alzheimer's disease", *Nature Neuroscience* 27,
  2366–2383. doi:10.1038/s41593-024-01774-5
- **Integrity:** byte size checked on download; SHA-256
  `3e3dac22446a8ce66afd07c209cd74df260b2f3d7c4085ffe4732390d2054a24`
  (verified against a reference download, pinned in
  `src/seaad_niches/download.py`) and recorded in the cache
  `MANIFEST.json` written by `scripts/download_seaad.py`.
- **Contents used:** per-cell table — spatial coordinates (x, y per
  section; stored in `obsm['X_spatial_raw']`), donor ID, section ID,
  mapped subclass/supertype annotations.

## snRNA-seq reference (MTG)

- **Dataset:** SEA-AD MTG snRNA-seq / snMultiome (processed, open).
- **Access paths:** AWS Open Data Registry (processed single-cell data
  link) and the CELLxGENE SEA-AD collection:
  https://cellxgene.cziscience.com/collections/1ca90a2d-2943-483d-b678-b809bf464c30
- **Citation / DOI:** same primary publication as above
  (doi:10.1038/s41593-024-01774-5).
- **Status:** not currently required by the pipeline — the MERFISH cells
  already carry subclass/supertype annotations mapped to the SEA-AD
  taxonomy. `io.load_snrnaseq_reference` remains a documented stub.

## Terms of use

Open data under the Allen Institute Terms of Use; cite both the primary
publication and the specific dataset. Raw sequencing data (AD Knowledge
Portal, controlled access) are **out of scope** for this repository.
