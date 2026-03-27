#!/usr/bin/env python3
"""Compare w2w output files (amn, mmn, eig) between test run and baseline.

Usage: python3 compare_outputs.py <test_dir> <baseline_dir> <prefix>
  e.g. python3 compare_outputs.py . baseline_bugfix GaAs-WANN
       python3 compare_outputs.py . baseline_bugfix WANN
"""
import sys, os

def compare_file(new_path, base_path, label):
    if not os.path.exists(new_path):
        print(f"  {label}: MISSING (new)")
        return False
    if not os.path.exists(base_path):
        print(f"  {label}: MISSING (baseline)")
        return False

    with open(new_path) as f:
        new_lines = f.readlines()
    with open(base_path) as f:
        base_lines = f.readlines()

    if len(new_lines) != len(base_lines):
        print(f"  {label}: LINE COUNT MISMATCH: {len(new_lines)} vs {len(base_lines)}")
        return False

    max_err = 0.0
    max_rel_err = 0.0
    ndiff = 0
    ntotal = 0
    max_err_info = ""
    first_diffs = []

    for i, (a, b) in enumerate(zip(new_lines, base_lines)):
        ta = a.split()
        tb = b.split()
        for j, (va, vb) in enumerate(zip(ta, tb)):
            try:
                fa, fb = float(va), float(vb)
                ntotal += 1
                err = abs(fa - fb)
                denom = max(abs(fa), abs(fb), 1e-30)
                rel = err / denom
                if err > 0:
                    ndiff += 1
                if err > max_err:
                    max_err = err
                    max_err_info = f"line {i+1}, col {j+1}: new={va} base={vb}"
                if rel > max_rel_err:
                    max_rel_err = rel
                if err > 1e-15 and len(first_diffs) < 5:
                    first_diffs.append(f"    line {i+1}: new={va} base={vb} (err={err:.2e})")
            except ValueError:
                pass

    bit_identical = (max_err == 0.0)
    status = "BIT-IDENTICAL" if bit_identical else "DIFFERS"
    print(f"  {label}: {status}")
    print(f"    lines={len(new_lines)}, values_compared={ntotal}")
    print(f"    max_abs_err={max_err:.2e}, max_rel_err={max_rel_err:.2e}, ndiff={ndiff}")
    if max_err > 0:
        print(f"    worst at: {max_err_info}")
    if first_diffs:
        print(f"    first differences:")
        for d in first_diffs:
            print(d)
    return bit_identical


def main():
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)

    test_dir = sys.argv[1]
    base_dir = sys.argv[2]
    prefix = sys.argv[3]

    print(f"Comparing {prefix} outputs:")
    print(f"  test_dir:  {test_dir}")
    print(f"  baseline:  {base_dir}")
    print()

    all_ok = True
    for ext in ['amn', 'mmn', 'eig']:
        new_path = os.path.join(test_dir, f"{prefix}.{ext}")
        base_path = os.path.join(base_dir, f"{prefix}.{ext}")
        ok = compare_file(new_path, base_path, ext.upper())
        all_ok = all_ok and ok
        print()

    if all_ok:
        print("RESULT: ALL FILES BIT-IDENTICAL")
    else:
        print("RESULT: DIFFERENCES DETECTED")

    return 0 if all_ok else 1


if __name__ == "__main__":
    main()
