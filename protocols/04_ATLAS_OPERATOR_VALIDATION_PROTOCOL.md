# Outcome-blind atlas operator registry protocol

Do not use HCP E2, CoRR, ABIDE LOCK outcomes, or target performance to choose atlas inclusion.

For every atlas: source/release/version; license; volumetric/surface space; template/resolution/affine; labels and semantics; cortical/subcortical scope; holes/overlap; interpolation policy; latent grid; parcel weights/normalization; operator SHA-256; coordinate compatibility; exclusion reason.

Validity tests: mass/weight conservation where expected, no unexplained empty rows, reproducible construction from specified source assets, grid/affine equality or specified transform, tolerance/resampling sensitivity, and independent reconstruction of sampled operators.

Eligible pair outputs: ranks, exact/non-identifiability status, row-space residuals, principal angles, missing dimensions, witnesses, certificate JSON, and stability range. Quotient equal-row-space representations before graph reduction. Do not claim an arbitrary finite atlas catalog is a lattice.
