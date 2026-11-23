# Paper 1 analytical dataset QC report

Criteria version: `paper1-freeze-v1`
Covariate table version: `paper1-donor-covariates-v1`
QC report version: `paper1-qc-v1`

## Inclusion and technical QC thresholds

- Donors require complete, donor-constant sex, age-at-death, PMI, and Continuous Pseudo-progression Score metadata.
- Donors require at least 1000 analyzed cells.
- Sections require at least 1000 cells.
- Cells must have finite x/y coordinates, a non-missing subclass, Genes detected >= 1.0, Number of spots >= 1.0, and Supertype confidence >= 0.0.
- SEA-AD ages recorded as `90+` are encoded as age 90 with `age_censored_90plus = True` in the donor covariate table.

## Outcome

- Donors included: 27 / 27
- Donors excluded: ['none']
- Sections passing QC: 69 / 69
- Sections failing QC: ['none']
- Cells passing cell-level QC: 1887728 / 1887729
- QC-passing cells in passing sections: 1887728 / 1887729

No downstream niche-progression model is fit in this QC report.
