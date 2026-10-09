"""ATI reference implementation v0.4.0-rc7."""
from .core import (
    validate_operator,array_sha256,numerical_rank,rowspace_basis,rowspace_projector,principal_angles,
    interoperability_deficit,find_kernel_witness,exact_recovery_analysis,classify_transfer,tolerance_sweep,
    exact_recovery_map,trace_bounded_projection_error_bound,minimal_augmentation,design_minimal_augmentation,
    universal_representation,design_universal_representation,build_interoperability_graph,gram_operator_rank_info,gram_containment_analysis,
    generate_certificate,
)
from .verifier import verify_certificate
from .exact_verifier import exact_rational_classify
from .identified_sets import scalar_variance_identified_interval,identified_interval_record,verify_interval_endpoint
from .stability import perturbation_stability_sweep
from .numerics import NUMERICAL_POLICY_ID
__version__="0.4.0rc7"
