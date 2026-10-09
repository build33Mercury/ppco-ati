# Visium HD cross-domain prospective protocol

Status: protocol only; no data analyzed here.

Question: can ATI certificates and augmentation decisions correctly predict which spatial count aggregations are reconstructable from another aggregation, using retained 2 um counts as truth?

Latent object: gene-specific nonnegative counts x_g on a specified 2 um grid. Source A is a predeclared binning/aggregation matrix. Target B is an independently fixed alternate binning or region aggregation. Observation y_A,g=A x_g and y_B,g=B x_g. This is a first-moment linear experiment, not automatically a covariance theorem.

Pre-specify before analysis: exact dataset/release/checksums; masks; genes; source/target bin geometry; held-out units; endpoints; predictive baseline if any; normalization.

Primary outcomes: certificate correctness against direct 2 um truth; exact-map error when ATI0; non-identifiability witness validation; minimal augmentation and reconstruction after augmentation; actual storage bytes and metadata for raw, source-only, and sufficient-sketch representations.

Invalidating conditions: nonlinear normalization, adaptive segmentation, inconsistent masks, or processing that invalidates the linear operator model.
