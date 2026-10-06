import json
from pathlib import Path
import numpy as np
import jsonschema
from ati.identified_sets import scalar_variance_identified_interval, identified_interval_record


def test_bounded_interval_public_fields_validate_schema():
    root=Path(__file__).resolve().parents[1]
    schema=json.loads((root/'schema'/'ati_scalar_interval_v0_2.schema.json').read_text())
    A=np.array([[1.,0.]])
    C=np.array([[1.]])
    b=np.array([1.,1.])
    r=identified_interval_record(scalar_variance_identified_interval(A,C,b,tau=2.0,return_endpoint_witnesses=False))
    jsonschema.validate(r,schema)
