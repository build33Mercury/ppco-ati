import json
from pathlib import Path
import numpy as np
import jsonschema
from ati import generate_certificate


def test_certificate_validates_against_schema():
    root=Path(__file__).resolve().parents[1]
    schema=json.loads((root/'schema'/'ati_certificate_v0_3.schema.json').read_text())
    A=np.eye(3); B=np.array([[1.,0.,0.]])
    cert=generate_certificate(A,B)
    jsonschema.validate(cert,schema)
