#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
WHEEL="$(find "$ROOT/dist" -maxdepth 1 -type f -name 'ati_reference-0.4.0rc7-*.whl' | head -n1)"
if [[ -z "$WHEEL" ]]; then
  echo "Missing prebuilt wheel in $ROOT/dist" >&2
  exit 2
fi
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
python -m venv --system-site-packages "$TMP/venv"
"$TMP/venv/bin/python" -m pip install "$WHEEL" --no-deps -q
cd "$TMP"
"$TMP/venv/bin/python" - <<'PY'
import ati, sys
from pathlib import Path
p=Path(ati.__file__).resolve()
assert str(p).startswith(str(Path(sys.prefix).resolve())), (p, sys.prefix)
print("INSTALLED_ATI_ORIGIN="+str(p))
PY
"$TMP/venv/bin/python" "$ROOT/examples/example_certificate.py"
"$TMP/venv/bin/python" - <<'PY'
import json, subprocess, sys, numpy as np
from pathlib import Path
A=np.eye(3); B=np.array([[1.,0.,0.]])
np.save('A_source.npy',A); np.save('B_target_operator.npy',B)
subprocess.check_call([str(Path(sys.executable).parent/'ati-certify'),'A_source.npy','B_target_operator.npy','--out','cert.json'])
x=json.loads(Path('cert.json').read_text())
assert x['verification']['independent_verifier_result']['verified']
C=np.array([[1.]])
b=np.array([1.,1.])
A2=np.array([[1.,0.]])
np.save('A_interval.npy',A2); np.save('C_interval.npy',C); np.save('target_vector.npy',b)
subprocess.check_call([str(Path(sys.executable).parent/'ati-identify-scalar'),'A_interval.npy','C_interval.npy','target_vector.npy','--tau','2','--out','interval.json'])
y=json.loads(Path('interval.json').read_text())
assert y['status']=='BOUNDED_TRACE_CONSTRAINED'
assert y['endpoint_verification']['lower_endpoint']['verified']
assert y['endpoint_verification']['upper_endpoint']['verified']
print('CLEANROOM_CLI_PASS')
PY
