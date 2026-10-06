import json
import numpy as np
from ati import generate_certificate, verify_certificate

A = np.array([[1.0, 0.0, 0.0]])
B = np.array([[0.0, 1.0, 0.0]])
cert = generate_certificate(A, B, source_id="toy_source", target_id="toy_target")
check = verify_certificate(cert, A, B)
print(json.dumps({"ati_class": cert["claim"]["ati_class"], "verified": check["verified"], "mode": check["mode"]}, indent=2))
