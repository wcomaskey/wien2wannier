#!/usr/bin/env python3
"""
wien2wannier Performance Benchmark Suite

This script benchmarks the w2w program before and after algorithmic optimizations.
It creates synthetic test cases of varying sizes and measures execution time.

Usage:
    python benchmark_w2w.py --build       # Build original and optimized versions
    python benchmark_w2w.py --run         # Run benchmarks
    python benchmark_w2w.py --compare     # Compare results and generate plots
    python benchmark_w2w.py --all         # Do all of the above
    python benchmark_w2w.py --test-lo     # Test LO coefficient module changes only

Requirements:
    - gfortran or ifort compiler
    - LAPACK/BLAS libraries
    - Python 3.6+ with matplotlib, numpy (for plotting)

Installation of dependencies (Ubuntu/Debian):
    sudo apt-get install gfortran liblapack-dev libblas-dev make

Installation of dependencies (macOS with Homebrew):
    brew install gcc lapack

Installation of dependencies (Windows with MSYS2):
    pacman -S mingw-w64-x86_64-gcc-fortran mingw-w64-x86_64-lapack mingw-w64-x86_64-openblas
"""

import os
import sys
import subprocess
import time
import argparse
import shutil
import re
import json
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple, Optional

# Configuration
def resolve_path(p: Path) -> Path:
    """Resolve path, handling WSL path conversion if needed."""
    s = str(p)
    # Check if running in WSL with a Windows-style path
    if sys.platform.startswith('linux') and len(s) > 2 and s[1] == ':':
        # Convert C:\path\to\file to /mnt/c/path/to/file
        drive = s[0].lower()
        rest = s[2:].replace('\\', '/')
        return Path(f'/mnt/{drive}{rest}')
    return p

SCRIPT_DIR = resolve_path(Path(__file__).parent.resolve())
PROJECT_ROOT = SCRIPT_DIR.parent
SRC_W2W = PROJECT_ROOT / "SRC_w2w"
RESULTS_DIR = SCRIPT_DIR / "benchmark_results"

# Versions to compare
VERSIONS = {
    "original": {
        "tag": "original",
        "description": "Original implementation (before optimizations)",
    },
    "optimized": {
        "tag": "optimized",
        "description": "Optimized implementation (after all fixes)",
    },
}

# Test configurations (varying problem sizes)
TEST_CONFIGS = [
    {"name": "tiny",   "nk": 8,   "nb": 10, "lmax": 2, "nneq": 1},
    {"name": "small",  "nk": 27,  "nb": 20, "lmax": 3, "nneq": 2},
    {"name": "medium", "nk": 64,  "nb": 40, "lmax": 4, "nneq": 3},
    {"name": "large",  "nk": 125, "nb": 60, "lmax": 5, "nneq": 4},
]


class BenchmarkRunner:
    """Manages benchmark execution and result collection."""

    def __init__(self, verbose: bool = True):
        self.verbose = verbose
        self.results: Dict = {}
        RESULTS_DIR.mkdir(exist_ok=True)

    def log(self, msg: str):
        if self.verbose:
            print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}")

    def run_command(self, cmd: List[str], cwd: Path = None, timeout: int = 600) -> Tuple[int, str, str]:
        """Run a command and return (returncode, stdout, stderr)."""
        self.log(f"Running: {' '.join(cmd)}")
        try:
            result = subprocess.run(
                cmd,
                cwd=cwd or PROJECT_ROOT,
                capture_output=True,
                text=True,
                timeout=timeout
            )
            return result.returncode, result.stdout, result.stderr
        except subprocess.TimeoutExpired:
            return -1, "", "Command timed out"
        except Exception as e:
            return -1, "", str(e)

    def build_version(self, version: str, extra_flags: str = "") -> bool:
        """Build a specific version of w2w."""
        self.log(f"Building {version} version...")

        build_dir = RESULTS_DIR / f"build_{version}"
        build_dir.mkdir(exist_ok=True)

        # Copy source files
        for f in SRC_W2W.glob("*.[fF]"):
            shutil.copy(f, build_dir)
        shutil.copy(SRC_W2W / "Makefile", build_dir)

        # Copy make.sys if it exists
        make_sys = PROJECT_ROOT / "make.sys"
        if not make_sys.exists():
            make_sys = PROJECT_ROOT / "make.sys.example"
        if make_sys.exists():
            shutil.copy(make_sys, build_dir.parent / "make.sys")

        # Clean and build
        self.run_command(["make", "clean"], cwd=build_dir)
        ret, stdout, stderr = self.run_command(
            ["make", f"FOPT=-O2 {extra_flags}".strip()],
            cwd=build_dir,
            timeout=300
        )

        if ret != 0:
            self.log(f"Build failed for {version}:")
            self.log(stderr)
            return False

        # Check executable exists
        exe = build_dir / "w2w"
        if not exe.exists():
            self.log(f"Executable not found: {exe}")
            return False

        self.log(f"Build successful: {exe}")
        return True

    def create_test_case(self, config: Dict, test_dir: Path) -> bool:
        """Create synthetic test input files for a given configuration."""
        test_dir.mkdir(parents=True, exist_ok=True)

        nk = config["nk"]
        nb = config["nb"]
        lmax = config["lmax"]
        nneq = config["nneq"]

        # Create minimal struct file
        struct_content = f"""test structure
P   LATTICE,NONEQUIV.ATOMS:{nneq:3d}
MODE OF CALC=RELA
  5.000000  5.000000  5.000000 90.000000 90.000000 90.000000
ATOM  -1: X=0.00000000 Y=0.00000000 Z=0.00000000
          MULT= 1          ISPLIT= 2
Si         NPT= 781  R0=0.00010000 RMT= 2.00000000   Z: 14.0
LOCAL ROT MATRIX:    1.0000000 0.0000000 0.0000000
                     0.0000000 1.0000000 0.0000000
                     0.0000000 0.0000000 1.0000000
"""
        for i in range(1, nneq):
            struct_content += f"""ATOM  -{i+1}: X={0.25*i:.8f} Y={0.25*i:.8f} Z={0.25*i:.8f}
          MULT= 1          ISPLIT= 2
Si         NPT= 781  R0=0.00010000 RMT= 2.00000000   Z: 14.0
LOCAL ROT MATRIX:    1.0000000 0.0000000 0.0000000
                     0.0000000 1.0000000 0.0000000
                     0.0000000 0.0000000 1.0000000
"""
        struct_content += f"""   {nneq} NUMBER OF SYMMETRY OPERATIONS
 1 0 0 0.00000000
 0 1 0 0.00000000
 0 0 1 0.00000000
       1
"""
        (test_dir / "test.struct").write_text(struct_content)

        # Create inwf file
        inwf_content = f"""Mmn Amn
{nb} {nb+5}
1
0 0.0 0.0
"""
        (test_dir / "test.inwf").write_text(inwf_content)

        # Create nnkp file (k-point neighbors)
        nnkp_content = f"""begin kpoints
{nk}
"""
        for i in range(nk):
            kx = (i % 3) / 3.0
            ky = ((i // 3) % 3) / 3.0
            kz = (i // 9) / 3.0
            nnkp_content += f"  {kx:.10f}  {ky:.10f}  {kz:.10f}\n"

        nnkp_content += """end kpoints

begin nnkpts
8
"""
        # Add k-point neighbor pairs
        pair_idx = 0
        for k1 in range(1, nk + 1):
            for b in range(8):  # 8 neighbors
                k2 = ((k1 + b - 1) % nk) + 1
                bqx = b % 2
                bqy = (b // 2) % 2
                bqz = b // 4
                nnkp_content += f"  {k1}  {k2}  {bqx}  {bqy}  {bqz}\n"
                pair_idx += 1
        nnkp_content += "end nnkpts\n"
        (test_dir / "test.nnkp").write_text(nnkp_content)

        # Create energy file
        energy_content = ""
        for i in range(nneq):
            energy_content += f"ATOM {i+1}\nDATA\n"
        for k in range(nk):
            kx = (k % 3) / 3.0
            ky = ((k // 3) % 3) / 3.0
            kz = (k // 9) / 3.0
            energy_content += f"   {kx:.5f}   {ky:.5f}   {kz:.5f}         "
            energy_content += f"        1       1    1.000000                  {100:6d}{nb:6d}\n"
            for b in range(1, nb + 1):
                energy_content += f"  {b}  {-0.5 + 0.1*b:12.7f}\n"
        (test_dir / "test.energy").write_text(energy_content)

        # Create fermi file
        (test_dir / "test.fermi").write_text("0.0\n")

        # Create vsp file (spherical potential)
        vsp_content = f"HEADER LINE               {lmax}\nDATA\n\n"
        for i in range(nneq):
            vsp_content += f"ATOM {i+1}\n"
            vsp_content += "   " + "   ".join([f"{-10.0 + 0.01*j:19.12E}" for j in range(4)]) + "\n" * 200
            vsp_content += "\n\n\n\n"
        (test_dir / "test.vsp").write_text(vsp_content)

        # Create def file
        def_content = f"""  5,'test.inwf',      'old',    'formatted',0
  6,'test.outputwf',  'unknown','formatted',0
  7,'test.amn',       'unknown','formatted',0
  8,'test.mmn',       'unknown','formatted',0
 10,'test.vector',    'unknown','unformatted',9000
 11,'test.nnkp',      'old',    'formatted',0
 12,'test.eig',       'unknown','formatted',0
 18,'test.vsp',       'old',    'formatted',0
 20,'test.struct',    'old',    'formatted',0
 50,'test.energy',    'old',    'formatted',0
 51,'test.fermi',     'old',    'formatted',0
"""
        (test_dir / "test.def").write_text(def_content)

        self.log(f"Created test case: {config['name']} (nk={nk}, nb={nb})")
        return True

    def parse_timing(self, output: str) -> Dict[str, float]:
        """Parse timing information from w2w output."""
        timings = {}

        # Look for ptime output format: "Times for XXX (sec): Y.YYY wall; Z.ZZZ CPU"
        pattern = r"Times for\s+(\w+)\s+\(sec\):\s+([\d.]+)\s+wall;\s+([\d.]+)\s+CPU"
        for match in re.finditer(pattern, output):
            name = match.group(1)
            wall_time = float(match.group(2))
            cpu_time = float(match.group(3))
            timings[name] = {"wall": wall_time, "cpu": cpu_time}

        # Look for internal l2Mmn timing
        pattern = r"CPU time used for atom\s+(\d+)\s+=\s+([\d.E+-]+)\s+([\d.E+-]+)\s+([\d.E+-]+)\s+([\d.E+-]+)"
        atom_times = []
        for match in re.finditer(pattern, output):
            atom_times.append({
                "atom": int(match.group(1)),
                "total": float(match.group(2)),
                "almgen": float(match.group(3)),
                "radint_ylm": float(match.group(4)),
                "overlap": float(match.group(5)),
            })
        if atom_times:
            timings["atom_breakdown"] = atom_times

        return timings

    def run_benchmark(self, version: str, config: Dict, num_runs: int = 3) -> Optional[Dict]:
        """Run benchmark for a specific version and configuration."""
        exe = RESULTS_DIR / f"build_{version}" / "w2w"
        if not exe.exists():
            self.log(f"Executable not found: {exe}")
            return None

        test_dir = RESULTS_DIR / f"test_{config['name']}"
        if not test_dir.exists():
            self.create_test_case(config, test_dir)

        results = []
        for run in range(num_runs):
            self.log(f"  Run {run+1}/{num_runs}...")

            # Clean previous output
            for f in test_dir.glob("test.outputwf"):
                f.unlink()

            start_time = time.time()
            ret, stdout, stderr = self.run_command(
                [str(exe), "test.def"],
                cwd=test_dir,
                timeout=600
            )
            elapsed = time.time() - start_time

            if ret != 0:
                self.log(f"  Run failed: {stderr[:200]}")
                continue

            # Read output file for detailed timing
            output_file = test_dir / "test.outputwf"
            output = ""
            if output_file.exists():
                output = output_file.read_text()

            timings = self.parse_timing(output)
            timings["total_elapsed"] = elapsed
            results.append(timings)

        if not results:
            return None

        # Average results
        avg_result = {
            "version": version,
            "config": config,
            "num_runs": len(results),
            "total_elapsed_avg": sum(r.get("total_elapsed", 0) for r in results) / len(results),
            "total_elapsed_min": min(r.get("total_elapsed", 0) for r in results),
            "total_elapsed_max": max(r.get("total_elapsed", 0) for r in results),
            "runs": results,
        }

        return avg_result

    def run_all_benchmarks(self) -> Dict:
        """Run benchmarks for all versions and configurations."""
        all_results = {}

        for version in VERSIONS:
            self.log(f"\n=== Benchmarking {version} ===")
            all_results[version] = {}

            for config in TEST_CONFIGS:
                self.log(f"\nConfig: {config['name']}")
                result = self.run_benchmark(version, config)
                if result:
                    all_results[version][config['name']] = result
                    self.log(f"  Average time: {result['total_elapsed_avg']:.3f}s")

        # Save results
        results_file = RESULTS_DIR / f"benchmark_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        with open(results_file, 'w') as f:
            json.dump(all_results, f, indent=2, default=str)
        self.log(f"\nResults saved to: {results_file}")

        return all_results

    def compare_results(self, results: Dict = None) -> str:
        """Generate comparison report from benchmark results."""
        if results is None:
            # Load most recent results
            result_files = sorted(RESULTS_DIR.glob("benchmark_*.json"))
            if not result_files:
                return "No benchmark results found."
            with open(result_files[-1]) as f:
                results = json.load(f)

        report = []
        report.append("=" * 70)
        report.append("WIEN2WANNIER PERFORMANCE COMPARISON")
        report.append("=" * 70)
        report.append("")

        # Table header
        report.append(f"{'Config':<12} {'Original':>12} {'Optimized':>12} {'Speedup':>10}")
        report.append("-" * 50)

        for config in TEST_CONFIGS:
            name = config['name']
            orig = results.get('original', {}).get(name, {}).get('total_elapsed_avg', 0)
            opt = results.get('optimized', {}).get(name, {}).get('total_elapsed_avg', 0)

            if orig > 0 and opt > 0:
                speedup = orig / opt
                report.append(f"{name:<12} {orig:>10.3f}s {opt:>10.3f}s {speedup:>9.2f}x")
            elif orig > 0:
                report.append(f"{name:<12} {orig:>10.3f}s {'N/A':>12} {'N/A':>10}")
            elif opt > 0:
                report.append(f"{name:<12} {'N/A':>12} {opt:>10.3f}s {'N/A':>10}")

        report.append("")
        report.append("=" * 70)

        return "\n".join(report)

    def generate_plots(self, results: Dict = None):
        """Generate comparison plots (requires matplotlib)."""
        try:
            import matplotlib.pyplot as plt
            import numpy as np
        except ImportError:
            self.log("matplotlib/numpy not available, skipping plots")
            return

        if results is None:
            result_files = sorted(RESULTS_DIR.glob("benchmark_*.json"))
            if not result_files:
                return
            with open(result_files[-1]) as f:
                results = json.load(f)

        configs = [c['name'] for c in TEST_CONFIGS]
        nk_values = [c['nk'] for c in TEST_CONFIGS]

        orig_times = []
        opt_times = []

        for config in TEST_CONFIGS:
            name = config['name']
            orig_times.append(results.get('original', {}).get(name, {}).get('total_elapsed_avg', 0))
            opt_times.append(results.get('optimized', {}).get(name, {}).get('total_elapsed_avg', 0))

        # Bar chart comparison
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        x = np.arange(len(configs))
        width = 0.35

        ax1.bar(x - width/2, orig_times, width, label='Original', color='steelblue')
        ax1.bar(x + width/2, opt_times, width, label='Optimized', color='coral')
        ax1.set_xlabel('Test Configuration')
        ax1.set_ylabel('Time (seconds)')
        ax1.set_title('Execution Time Comparison')
        ax1.set_xticks(x)
        ax1.set_xticklabels(configs)
        ax1.legend()
        ax1.grid(axis='y', alpha=0.3)

        # Scaling plot
        ax2.loglog(nk_values, orig_times, 'o-', label='Original', color='steelblue')
        ax2.loglog(nk_values, opt_times, 's-', label='Optimized', color='coral')
        ax2.set_xlabel('Number of k-points')
        ax2.set_ylabel('Time (seconds)')
        ax2.set_title('Scaling with Problem Size')
        ax2.legend()
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        plot_file = RESULTS_DIR / f"benchmark_plot_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
        plt.savefig(plot_file, dpi=150)
        self.log(f"Plot saved to: {plot_file}")
        plt.close()


def check_prerequisites() -> Tuple[bool, List[str]]:
    """Check if all build prerequisites are available."""
    issues = []
    is_windows = sys.platform.startswith('win')

    def find_executable(name: str) -> bool:
        """Cross-platform executable finder."""
        return shutil.which(name) is not None

    # Check for Fortran compiler
    fc_found = False
    for fc in ['gfortran', 'ifort']:
        if find_executable(fc):
            fc_found = True
            break
    if not fc_found:
        issues.append("No Fortran compiler found (gfortran or ifort)")

    # Check for make
    make_names = ['make', 'mingw32-make'] if is_windows else ['make']
    make_found = any(find_executable(m) for m in make_names)
    if not make_found:
        issues.append("'make' not found")

    # Check for BLAS/LAPACK (best effort, platform-dependent)
    blas_found = False

    if is_windows:
        # On Windows, check common MSYS2/MinGW paths
        mingw_paths = [
            'C:/msys64/mingw64/lib',
            'C:/msys64/ucrt64/lib',
            'C:/mingw64/lib',
        ]
        for lib_path in mingw_paths:
            if os.path.exists(lib_path):
                files = os.listdir(lib_path) if os.path.isdir(lib_path) else []
                if any('blas' in f.lower() or 'lapack' in f.lower() for f in files):
                    blas_found = True
                    break
        # Also check if openblas is in PATH
        if not blas_found and find_executable('libopenblas.dll'):
            blas_found = True
    else:
        # Unix-like systems
        for lib_path in ['/usr/lib', '/usr/local/lib', '/lib', '/usr/lib/x86_64-linux-gnu']:
            if os.path.exists(lib_path):
                try:
                    files = os.listdir(lib_path)
                    if any('blas' in f.lower() or 'lapack' in f.lower() for f in files):
                        blas_found = True
                        break
                except PermissionError:
                    continue

        # Also check via ldconfig if available (Linux only)
        if not blas_found and find_executable('ldconfig'):
            result = subprocess.run(
                'ldconfig -p 2>/dev/null | grep -E "(blas|lapack)"',
                shell=True, capture_output=True, text=True
            )
            if result.returncode == 0 and result.stdout.strip():
                blas_found = True

    if not blas_found:
        if is_windows:
            issues.append("BLAS/LAPACK libraries not found (install via MSYS2: pacman -S mingw-w64-x86_64-openblas)")
        else:
            issues.append("BLAS/LAPACK libraries not found (install liblapack-dev libblas-dev)")

    return len(issues) == 0, issues


def test_lo_coefficients():
    """Test the LO coefficient fix by checking module compilation."""
    print("\n" + "=" * 60)
    print("TEST: Local Orbital Coefficient Module")
    print("=" * 60)

    # Check that the module has been updated
    modules_file = SRC_W2W / "modules.f"
    if not modules_file.exists():
        print("ERROR: modules.f not found")
        return False

    content = modules_file.read_text()

    # Check for required changes
    checks = [
        ("allocatable pi12lo", "pi12lo(:,:)" in content or "allocatable.*pi12lo" in content),
        ("allocatable pe12lo", "pe12lo(:,:)" in content or "allocatable.*pe12lo" in content),
        ("allocatable pilolo", "pilolo(:,:,:)" in content or "allocatable.*pilolo" in content),
        ("init_loabc subroutine", "subroutine init_loabc" in content),
        ("cleanup_loabc subroutine", "subroutine cleanup_loabc" in content),
    ]

    all_passed = True
    for name, passed in checks:
        status = "PASS" if passed else "FAIL"
        print(f"  [{status}] {name}")
        if not passed:
            all_passed = False

    # Try to compile
    print("\n  Attempting compilation...")
    build_dir = RESULTS_DIR / "build_test_lo"
    build_dir.mkdir(parents=True, exist_ok=True)

    # Copy files
    for f in SRC_W2W.glob("*.[fF]"):
        shutil.copy(f, build_dir)
    shutil.copy(SRC_W2W / "Makefile", build_dir)

    make_sys = PROJECT_ROOT / "make.sys"
    if not make_sys.exists():
        make_sys = PROJECT_ROOT / "make.sys.example"
    if make_sys.exists():
        shutil.copy(make_sys, build_dir.parent / "make.sys")

    # Try compilation
    result = subprocess.run(
        ["make", "clean"],
        cwd=build_dir,
        capture_output=True,
        text=True
    )

    result = subprocess.run(
        ["make"],
        cwd=build_dir,
        capture_output=True,
        text=True,
        timeout=300
    )

    if result.returncode == 0:
        print("  [PASS] Compilation successful")
    else:
        print("  [FAIL] Compilation failed:")
        print(result.stderr[:500])
        all_passed = False

    return all_passed


def main():
    parser = argparse.ArgumentParser(description="wien2wannier benchmark suite")
    parser.add_argument("--build", action="store_true", help="Build all versions")
    parser.add_argument("--run", action="store_true", help="Run benchmarks")
    parser.add_argument("--compare", action="store_true", help="Compare results")
    parser.add_argument("--plot", action="store_true", help="Generate plots")
    parser.add_argument("--test-lo", action="store_true", help="Test LO coefficient fix")
    parser.add_argument("--check", action="store_true", help="Check build prerequisites")
    parser.add_argument("--all", action="store_true", help="Run all steps")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")

    args = parser.parse_args()

    if not any([args.build, args.run, args.compare, args.plot, args.test_lo, args.check, args.all]):
        parser.print_help()
        return

    runner = BenchmarkRunner(verbose=args.verbose or args.all)

    # Check prerequisites if requested or before building
    if args.check or args.build or args.all:
        print("\n" + "=" * 60)
        print("CHECKING BUILD PREREQUISITES")
        print("=" * 60)
        ok, issues = check_prerequisites()
        if ok:
            print("  [PASS] All prerequisites found")
        else:
            print("  [WARN] Missing prerequisites:")
            for issue in issues:
                print(f"    - {issue}")
            print("\n  To install on Ubuntu/Debian/WSL:")
            print("    sudo apt-get install gfortran liblapack-dev libblas-dev make")
            print("\n  To install on macOS:")
            print("    brew install gcc lapack")
            print("\n  To install on Windows (MSYS2):")
            print("    pacman -S mingw-w64-x86_64-gcc-fortran mingw-w64-x86_64-openblas make")
            if not args.check:
                print("\n  Continuing anyway (build may fail)...")

        if args.check:
            return

    if args.test_lo or args.all:
        test_lo_coefficients()

    if args.build or args.all:
        runner.log("\n=== Building versions ===")
        for version in VERSIONS:
            runner.build_version(version)

    if args.run or args.all:
        results = runner.run_all_benchmarks()
    else:
        results = None

    if args.compare or args.all:
        print(runner.compare_results(results))

    if args.plot or args.all:
        runner.generate_plots(results)


if __name__ == "__main__":
    main()
