from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
from .core import generate_certificate
from .verifier import verify_certificate


def main(argv=None):
    p = argparse.ArgumentParser(description="Generate and independently verify a fail-closed ATI numerical certificate from NumPy operator files.")
    p.add_argument("source", help="source operator .npy")
    p.add_argument("target", help="target operator .npy")
    p.add_argument("--source-id", default="source")
    p.add_argument("--target-id", default="target")
    p.add_argument("--rtol", type=float, default=1e-10)
    p.add_argument("--atol", type=float, default=1e-12)
    p.add_argument("--out", required=True, help="output JSON certificate")
    args = p.parse_args(argv)
    A = np.load(args.source, allow_pickle=False)
    B = np.load(args.target, allow_pickle=False)
    cert = generate_certificate(A, B, args.source_id, args.target_id, args.rtol, args.atol)
    check = verify_certificate(cert, A, B)
    cert["verification"]["independent_verifier_result"] = check
    if not check["verified"]:
        raise SystemExit(f"certificate failed independent verification: {check}")
    out = Path(args.out)
    out.write_text(json.dumps(cert, indent=2, sort_keys=True), encoding="utf-8")
    print(out)


if __name__ == "__main__":
    main()
