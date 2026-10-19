#!/usr/bin/env python3
"""Behavioral tests for the public wien2wannier orchestration driver."""

from __future__ import annotations

import json
import importlib.machinery
import importlib.util
import os
import re
from pathlib import Path
import stat
import subprocess
import sys
import tempfile
import textwrap
import unittest


ROOT = Path(__file__).resolve().parents[1]
CLI = ROOT / "SRC" / "wien2wannier"


class Wien2WannierCliTests(unittest.TestCase):
    def setUp(self) -> None:
        self._temp = tempfile.TemporaryDirectory()
        self.root = Path(self._temp.name)
        self.case_dir = self.root / "Si"
        self.case_dir.mkdir()
        self.out = self.root / "result"
        self.bin_dir = self.root / "bin"
        self.bin_dir.mkdir()
        self.log = self.root / "calls.jsonl"
        self._write_case("Si")
        self._write_fake_backends()

    def tearDown(self) -> None:
        self._temp.cleanup()

    def _write_case(self, case: str, *, nnkp: bool = True) -> None:
        for suffix, data in {
            "struct": (
                "test structure\n"
                "P   LATTICE,NONEQUIV.ATOMS:  0\n"
                "MODE OF CALC=RELA\n"
                "  1.000000  1.000000  1.000000 90.000000 90.000000 90.000000\n"
            ),
            "vector": "vector-data\n",
            "energy": "energy-data\n",
            "vsp": "potential-data\n",
            "fermi": "0.125\n",
        }.items():
            (self.case_dir / f"{case}.{suffix}").write_text(data)
        (self.case_dir / f"{case}.klist").write_text(
            "         1         0         0         0         2  1.0"
            " -7.0  1.5         0 k, div: ( 2 1 1)\n"
            "         2         1         0         0         2  1.0\n"
            "END\n"
        )
        (self.case_dir / f"{case}.win").write_text(
            "num_bands = 8\nnum_wann = 4\nwrite_hr = true\n"
            "dis_num_iter = 200\n"
        )
        if nnkp:
            (self.case_dir / f"{case}.nnkp").write_text(
                "begin recip_lattice\n1 0 0\n0 1 0\n0 0 1\nend recip_lattice\n"
                "begin kpoints\n2\n0.0 0.0 0.0\n0.5 0.0 0.0\n"
                "end kpoints\nbegin nnkpts\n1\n1 2 0 0 0\n2 1 1 0 0\n"
                "end nnkpts\n"
            )

    def _executable(self, name: str, body: str) -> Path:
        path = self.bin_dir / name
        path.write_text("#!/usr/bin/env python3\n" + textwrap.dedent(body))
        path.chmod(path.stat().st_mode | stat.S_IXUSR)
        return path

    def _write_fake_backends(self) -> None:
        common_log = """
            import json, os, pathlib, sys
            log = pathlib.Path(os.environ['W2W_TEST_LOG'])
            with log.open('a') as stream:
                stream.write(json.dumps({
                    'exe': pathlib.Path(sys.argv[0]).name,
                    'argv': sys.argv[1:],
                    'cwd': os.getcwd(),
                    'omp': os.environ.get('OMP_NUM_THREADS'),
                    'openblas': os.environ.get('OPENBLAS_NUM_THREADS'),
                    'edges': os.environ.get('W2W_EDGES'),
                    'seedk': os.environ.get('W2W_SEEDK'),
                    'seedloss': os.environ.get('W2W_SEEDLOSS'),
                    'struct_env': os.environ.get('W2W_STRUCT'),
                }) + '\\n')
        """
        self._executable(
            "w2w",
            common_log
            + """
            deffile = pathlib.Path(sys.argv[-1])
            lines = deffile.read_text().splitlines()
            inwf = next(line.split("'", 2)[1] for line in lines
                        if line.lstrip().startswith("5,"))
            seed = pathlib.Path(inwf).stem
            mode = pathlib.Path(inwf).read_text().splitlines()[0]
            pathlib.Path(seed + '.outputwf').write_text('mock output\\n')
            pathlib.Path(seed + '.amn').write_text(
                'mock raw\\n           3           2           2\\n')
            pathlib.Path(seed + '.mmn').write_text('mock raw mmn\\n')
            pathlib.Path(seed + '.eig').write_text('mock raw eig\\n')
            pathlib.Path(seed + '.error').write_text('')
            if mode == 'HYBR':
                pathlib.Path(seed + '-hyb.amn').write_text(
                    'mock handoff\\n           2           2           2\\n')
                pathlib.Path(seed + '-hyb.mmn').write_text('mock mmn\\n')
                pathlib.Path(seed + '-hyb.eig').write_text('mock eig\\n')
                pathlib.Path(seed + '-hyb.anchors').write_text('mock anchors\\n')
            """,
        )
        self._executable("w2wc", (self.bin_dir / "w2w").read_text().split("\n", 1)[1])
        self._executable(
            "w2w_hybrid",
            common_log
            + """
            # supported named interface: --seed --nw --out --pool-factor --solver
            args = sys.argv[1:]
            if '--continuity-chk' in args:
                print('  subspace continuity: mock 0.9 for ' + args[-1])
                sys.exit(0)
            out = args[args.index('--out') + 1]
            pathlib.Path(out + '.amn').write_text(
                'mock bridge handoff\\n           3           2           3\\n')
            pathlib.Path(out + '.mmn').write_text('mock mmn\\n')
            pathlib.Path(out + '.eig').write_text('mock eig\\n')
            if '--window' in args:
                pathlib.Path(out + '-win.amn').write_text('mock window\\n           3           2           2\\n')
                pathlib.Path(out + '-win.mmn').write_text('mock mmn\\n')
                pathlib.Path(out + '-win.eig').write_text('mock eig\\n')
                pathlib.Path(out + '-win.win').write_text('num_bands = 3\\nnum_wann = 2\\n')
            """,
        )
        self._executable(
            "write_win_backend",
            common_log
            + """
            print('num_bands = 1')
            print('num_wann = 1')
            print('write_hr = true')
            print('begin unit_cell_cart')
            print('bohr')
            print('1 0 0\\n0 1 0\\n0 0 1')
            print('end unit_cell_cart')
            print('begin kpoints')
            print('0 0 0\\n0.5 0 0')
            print('end kpoints')
            """,
        )
        self._executable(
            "wannier90.x",
            common_log
            + """
            args = sys.argv[1:]
            if args and args[0] == '-pp':
                seed = args[1]
                pathlib.Path(seed + '.nnkp').write_text(
                    'begin recip_lattice\\n1 0 0\\n0 1 0\\n0 0 1\\nend recip_lattice\\n'
                    'begin kpoints\\n2\\n0 0 0\\n0.5 0 0\\nend kpoints\\n'
                    'begin nnkpts\\n1\\n1 2 0 0 0\\n2 1 1 0 0\\nend nnkpts\\n')
            else:
                pathlib.Path(args[-1] + '.wout').write_text('mock localized\\n')
            """,
        )

    def run_cli(self, *args: str, env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
        run_env = os.environ.copy()
        run_env["W2W_TEST_LOG"] = str(self.log)
        run_env["PATH"] = str(self.bin_dir) + os.pathsep + run_env.get("PATH", "")
        if env:
            run_env.update(env)
        return subprocess.run(
            [str(CLI), *args],
            cwd=self.root,
            env=run_env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )

    def calls(self) -> list[dict[str, object]]:
        if not self.log.exists():
            return []
        return [json.loads(line) for line in self.log.read_text().splitlines()]

    def common_args(self, *extra: str) -> tuple[str, ...]:
        return (
            "--input", str(self.case_dir),
            "--seed", "model",
            "--output-dir", str(self.out),
            "--backend-dir", str(self.bin_dir),
            *extra,
        )

    def test_help_and_version_need_no_scientific_environment(self) -> None:
        help_run = self.run_cli("--help", env={"PATH": "/usr/bin:/bin"})
        self.assertEqual(help_run.returncode, 0, help_run.stderr)
        for heading in (
            "Input and output", "Automatic selection", "Run stages",
            "Resources", "Diagnostics", "Advanced backend options",
        ):
            self.assertIn(heading, help_run.stdout)
        for option in (
            "--input", "--seed", "--output-dir", "--method", "--num-wann",
            "--k-grid", "--stage", "--wannier90", "--threads", "--parallel",
            "--mpi-ranks", "--mpi-bind", "--check", "--dry-run", "--verbose", "--quiet",
            "--version",
        ):
            self.assertIn(option, help_run.stdout)
        version_run = self.run_cli("--version", env={"PATH": "/usr/bin:/bin"})
        self.assertEqual(version_run.returncode, 0, version_run.stderr)
        self.assertRegex(version_run.stdout, r"wien2wannier \d+\.\d+")
        self.assertRegex(version_run.stdout, r"source revision: (?!unavailable)")
        self.assertIn("supported:", version_run.stdout)
        self.assertIn("detected backends:", version_run.stdout)

    def test_check_validates_mesh_without_writing_or_executing(self) -> None:
        result = self.run_cli(*self.common_args("--check", "--k-grid", "2", "1", "1"))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("input kind: wien", result.stdout)
        self.assertIn("mesh: 2 x 1 x 1 (2 k-points)", result.stdout)
        self.assertIn("automatic num_wann", result.stdout)
        self.assertFalse(self.out.exists())
        self.assertEqual(self.calls(), [])

    def test_conflicting_num_wann_aliases_are_rejected(self) -> None:
        result = self.run_cli(*self.common_args("--check", "--nw", "4", "--num-wann", "5"))
        self.assertEqual(result.returncode, 2)
        self.assertIn("conflicting values", result.stderr)
        self.assertFalse(self.out.exists())

    def test_wrong_requested_mesh_fails_before_writes(self) -> None:
        result = self.run_cli(*self.common_args("--check", "--k-grid", "3", "1", "1"))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("does not match existing", result.stderr)
        self.assertFalse(self.out.exists())
        self.assertEqual(self.calls(), [])

    def test_same_count_permuted_nnkp_is_rejected(self) -> None:
        (self.case_dir / "Si.nnkp").write_text(
            "begin recip_lattice\n1 0 0\n0 1 0\n0 0 1\nend recip_lattice\n"
            "begin kpoints\n2\n0.5 0.0 0.0\n0.0 0.0 0.0\n"
            "end kpoints\nbegin nnkpts\n1\n1 2 1 0 0\n2 1 0 0 0\n"
            "end nnkpts\n"
        )
        result = self.run_cli(*self.common_args("--check"))
        self.assertEqual(result.returncode, 2)
        self.assertIn("coordinates/order", result.stderr)
        self.assertFalse(self.out.exists())

    def test_nnkp_of_a_mesh_with_nonterminating_fractions_is_accepted(self) -> None:
        # Three k-points along x, printed with eight decimals as Wannier90 does.
        # The neighbour vectors of the three k-points then differ by 6e-8 in
        # these units; an exact comparison of values rounded to 5e-8 refused it.
        loader = importlib.machinery.SourceFileLoader("w2w_cli_nnkp_module", str(CLI))
        spec = importlib.util.spec_from_loader(loader.name, loader)
        assert spec is not None
        driver = importlib.util.module_from_spec(spec)
        sys.modules[loader.name] = driver
        try:
            loader.exec_module(driver)
            nnkp = self.case_dir / "third.nnkp"
            nnkp.write_text(
                "begin recip_lattice\n6.2831853 0 0\n0 6.2831853 0\n0 0 6.2831853\nend recip_lattice\n"
                "begin kpoints\n3\n0.00000000 0.0 0.0\n0.33333333 0.0 0.0\n0.66666667 0.0 0.0\n"
                "end kpoints\nbegin nnkpts\n2\n"
                "1 2 0 0 0\n1 3 -1 0 0\n2 3 0 0 0\n2 1 0 0 0\n3 1 1 0 0\n3 2 0 0 0\n"
                "end nnkpts\n"
            )
            points, nntot = driver.inspect_nnkp(nnkp)
            self.assertEqual((len(points), nntot), (3, 2))
            # a neighbour list whose vectors really differ is still refused
            nnkp.write_text(nnkp.read_text().replace("3 1 1 0 0", "3 1 2 0 0"))
            with self.assertRaises(driver.CliError):
                driver.inspect_nnkp(nnkp)
        finally:
            sys.modules.pop(loader.name, None)

    def test_duplicate_klist_is_rejected_even_with_a_valid_div_header(self) -> None:
        lines = (self.case_dir / "Si.klist").read_text().splitlines()
        lines[1] = "         2         0         0         0         2  1.0"
        (self.case_dir / "Si.klist").write_text("\n".join(lines) + "\n")
        result = self.run_cli(*self.common_args("--check"))
        self.assertEqual(result.returncode, 2)
        self.assertIn("duplicate", result.stderr)
        self.assertFalse(self.out.exists())

    def test_nonfinite_and_incomplete_nnkp_records_are_rejected(self) -> None:
        reciprocal = "begin recip_lattice\n1 0 0\n0 1 0\n0 0 1\nend recip_lattice\n"
        malformed = (
            (
                reciprocal + "begin kpoints\n2\nnan 0 0\n0.5 0 0\nend kpoints\n"
                "begin nnkpts\n1\n1 2 0 0 0\n2 1 0 0 0\nend nnkpts\n",
                "finite",
            ),
            (
                reciprocal + "begin kpoints\n2\n0 0 0\n0.5 0 0\nend kpoints\n"
                "begin nnkpts\n1\n1 2 0 0 0\nend nnkpts\n",
                "neighbour",
            ),
        )
        for data, message in malformed:
            with self.subTest(message=message):
                (self.case_dir / "Si.nnkp").write_text(data)
                result = self.run_cli(*self.common_args("--check"))
                self.assertEqual(result.returncode, 2)
                self.assertIn(message, result.stderr)
                self.assertFalse(self.out.exists())

    def test_nnkp_neighbour_record_and_slot_order_are_validated(self) -> None:
        reciprocal = "begin recip_lattice\n1 0 0\n0 1 0\n0 0 1\nend recip_lattice\n"
        malformed = (
            (
                reciprocal + "begin kpoints\n2\n0 0 0\n0.5 0 0\nend kpoints\n"
                "begin nnkpts\n1\n2 1 0 0 0\n1 2 0 0 0\nend nnkpts\n",
                "record order",
            ),
            (
                reciprocal + "begin kpoints\n2\n0 0 0\n0.5 0 0\nend kpoints\n"
                "begin nnkpts\n2\n1 2 0 0 0\n1 1 0 1 0\n"
                "2 2 0 1 0\n2 1 1 0 0\nend nnkpts\n",
                "slot",
            ),
        )
        for data, message in malformed:
            with self.subTest(message=message):
                (self.case_dir / "Si.nnkp").write_text(data)
                result = self.run_cli(*self.common_args("--check"))
                self.assertEqual(result.returncode, 2)
                self.assertIn(message, result.stderr)
                self.assertFalse(self.out.exists())

    def test_native_coordinate_conversion_covers_centering_and_ortho_threshold(self) -> None:
        loader = importlib.machinery.SourceFileLoader("w2w_cli_test_module", str(CLI))
        spec = importlib.util.spec_from_loader(loader.name, loader)
        assert spec is not None
        module = importlib.util.module_from_spec(spec)
        sys.modules[loader.name] = module
        loader.exec_module(module)
        raw = [(module.Fraction(1, 8), module.Fraction(2, 8), module.Fraction(3, 8))]
        expected = {
            "F": (0.3125, 0.25, 0.1875),
            "B": (0.25, 0.125, 0.0),
            "CXY": (-0.0625, 0.1875, 0.375),
            "CYZ": (0.125, -0.0625, 0.3125),
            "CXZ": (-0.125, 0.25, 0.25),
        }
        for lattice, want in expected.items():
            struct = self.root / f"{lattice}.struct"
            struct.write_text(
                "fixture\n" + f"{lattice:<4} LATTICE,NONEQUIV.ATOMS:  0\n"
                + "MODE OF CALC=RELA\n"
                + "  1.000000  1.000000  1.000000 90.000000 90.000000 90.001000\n"
            )
            got = module.native_klist_coordinates(struct, raw)[0]
            for actual, reference in zip(got, want):
                self.assertAlmostEqual(actual, reference, places=12)
        pstruct = self.root / "P-near.struct"
        pstruct.write_text(
            "fixture\nP   LATTICE,NONEQUIV.ATOMS:  0\nMODE OF CALC=RELA\n"
            "  1.000000  1.000000  1.000000 90.000000 90.000000 90.001000\n"
        )
        got = module.native_klist_coordinates(pstruct, raw)[0]
        self.assertAlmostEqual(got[0], 0.12499563665783164, places=12)

    def test_soc_dry_run_resolves_extended_def_contract_without_writes(self) -> None:
        for suffix in ("vectorsoup", "vectorsodn", "energysoup", "vspup", "vspdn"):
            (self.case_dir / f"Si.{suffix}").write_text(suffix + "\n")
        (self.case_dir / "Si.inso").write_text("WFFIL\n4 0 0\n-10 1.5\n0 0 1\n")
        result = self.run_cli(
            *self.common_args("--spin", "soc", "--dry-run", "--k-grid", "2", "1", "1")
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("backend: w2wc", result.stdout)
        self.assertRegex(result.stdout, r"unit\s+9\s+.*vectorsodn")
        self.assertRegex(result.stdout, r"unit\s+19\s+.*vspdn")
        self.assertRegex(result.stdout, r"unit\s+22\s+.*inso")
        self.assertFalse(self.out.exists())
        self.assertEqual(self.calls(), [])

    def test_all_runs_integrated_auto_hybrid_and_publishes_complete_handoff(self) -> None:
        result = self.run_cli(
            *self.common_args(
                "--k-grid", "2", "1", "1", "--parallel", "threads", "--threads", "2"
            )
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        inwf = (self.out / "model.inwf").read_text()
        self.assertTrue(inwf.startswith("HYBR\n0 -1\n4\n"))
        self.assertIn("0.0 0 T T 0.5", inwf)
        self.assertIn("\n0\n", inwf)
        self.assertIn("froz_emax = 1.0", inwf)
        self.assertIn("solver = cg", inwf)
        self.assertIn("colored = false", inwf)
        deffile = (self.out / "model.def").read_text()
        self.assertIn(str(self.case_dir / "Si.vector"), deffile)
        self.assertIn("11,'model.nnkp'", deffile)
        calls = self.calls()
        self.assertEqual([call["exe"] for call in calls], ["w2w"])
        # An OpenMP-built BLAS reads only OMP_NUM_THREADS; w2w has no OpenMP regions.
        self.assertEqual(calls[0]["omp"], "2")
        self.assertEqual(calls[0]["openblas"], "2")
        handoff_win = (self.out / "model-hyb.win").read_text()
        self.assertIn("num_bands = 2", handoff_win)
        self.assertIn("num_wann = 2", handoff_win)
        self.assertNotIn("dis_num_iter", handoff_win)
        self.assertTrue((self.out / "model-hyb.amn").is_file())
        self.assertFalse((self.out / "model-hyb.wout").exists())

    def test_froz_emax_auto_reaches_the_hybr_keyword_with_explicit_num_wann(self) -> None:
        # the published release protocol: the ceiling walked down from E_F+8 eV (HYBR froz_emax = auto)
        result = self.run_cli(
            *self.common_args(
                "--k-grid", "2", "1", "1", "--num-wann", "2", "--froz-emax", "auto",
                "--parallel", "threads", "--threads", "2",
            )
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        inwf = (self.out / "model.inwf").read_text()
        self.assertIn("\n2\n", inwf)
        self.assertIn("froz_emax = auto", inwf)

    def test_froz_emax_words_need_num_wann(self) -> None:
        for word in ("auto", "none"):
            with self.subTest(word=word):
                result = self.run_cli(*self.common_args("--check", "--k-grid", "2", "1", "1", "--froz-emax", word))
                self.assertEqual(result.returncode, 2)
                self.assertIn("--num-wann", result.stderr)
                self.assertFalse(self.out.exists())
        result = self.run_cli(*self.common_args("--check", "--froz-emax", "high"))
        self.assertEqual(result.returncode, 2)
        self.assertIn("froz-emax", result.stderr)

    def test_explicit_wannier90_adds_localization_after_handoff(self) -> None:
        w90 = self.bin_dir / "wannier90.x"
        result = self.run_cli(*self.common_args("--wannier90", str(w90)))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([call["exe"] for call in self.calls()], ["w2w", "wannier90.x"])
        self.assertEqual(self.calls()[-1]["argv"], ["model-hyb"])
        self.assertTrue((self.out / "model-hyb.wout").is_file())

    def test_prepare_generates_neighbors_without_running_w2w(self) -> None:
        (self.case_dir / "Si.nnkp").unlink()
        (self.case_dir / "Si.win").unlink()
        result = self.run_cli(
            *self.common_args(
                "--stage", "prepare", "--wannier90-preprocess", str(self.bin_dir / "wannier90.x")
            )
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(
            [call["exe"] for call in self.calls()],
            ["write_win_backend", "wannier90.x"],
        )
        self.assertEqual(self.calls()[-1]["argv"], ["-pp", "model"])
        self.assertTrue((self.out / "model.nnkp").is_file())
        self.assertFalse((self.out / "model-hyb.amn").exists())

    def test_prepare_builds_missing_win_while_reusing_existing_neighbors(self) -> None:
        (self.case_dir / "Si.win").unlink()
        result = self.run_cli(*self.common_args("--stage", "prepare"))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([call["exe"] for call in self.calls()], ["write_win_backend"])
        self.assertTrue((self.out / "model.win").is_file())
        self.assertTrue((self.out / "model.nnkp").is_file())

    def test_backend_failure_stops_before_optional_localization(self) -> None:
        self._executable(
            "w2w",
            """
            import json, os, pathlib, sys
            with pathlib.Path(os.environ['W2W_TEST_LOG']).open('a') as stream:
                stream.write(json.dumps({'exe': 'w2w', 'argv': sys.argv[1:]}) + '\\n')
            raise SystemExit(7)
            """,
        )
        result = self.run_cli(
            *self.common_args("--wannier90", str(self.bin_dir / "wannier90.x"))
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("exit code 7", result.stderr)
        self.assertEqual([call["exe"] for call in self.calls()], ["w2w"])
        self.assertFalse((self.out / "model-hyb.wout").exists())

    def test_quiet_keeps_successful_backend_warnings_visible(self) -> None:
        self._executable(
            "w2w",
            """
            import pathlib, sys
            seed = pathlib.Path(sys.argv[-1]).stem
            print('routine iteration progress')
            print('WARNING: backend reached a guarded fallback')
            print('native library diagnostic', file=sys.stderr)
            pathlib.Path(seed + '.outputwf').write_text(
                'ordinary scientific progress\\nWARNING: integrated capacity fallback\\n'
                '  -- NOT converged at the requested sweep cap\\n')
            pathlib.Path(seed + '.error').write_text('')
            pathlib.Path(seed + '.amn').write_text('raw\\n3 2 2\\n')
            pathlib.Path(seed + '.mmn').write_text('raw\\n')
            pathlib.Path(seed + '.eig').write_text('raw\\n')
            pathlib.Path(seed + '-hyb.amn').write_text('hyb\\n2 2 2\\n')
            pathlib.Path(seed + '-hyb.mmn').write_text('hyb\\n')
            pathlib.Path(seed + '-hyb.eig').write_text('hyb\\n')
            pathlib.Path(seed + '-hyb.anchors').write_text('hyb\\n')
            """,
        )
        result = self.run_cli(*self.common_args("--quiet"))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("routine iteration progress", result.stdout + result.stderr)
        self.assertIn("WARNING: backend reached a guarded fallback", result.stderr)
        self.assertIn("native library diagnostic", result.stderr)
        self.assertIn("WARNING: integrated capacity fallback", result.stderr)
        self.assertIn("NOT converged", result.stderr)
        self.assertNotIn("ordinary scientific progress", result.stderr)

    def test_disentanglement_status_heads_the_handoff_win(self) -> None:
        self._executable(
            "w2w",
            """
            import pathlib, sys
            seed = pathlib.Path(sys.argv[-1]).stem
            pathlib.Path(seed + '.outputwf').write_text('progress\\n')
            pathlib.Path(seed + '.error').write_text('')
            pathlib.Path(seed + '.amn').write_text('raw\\n3 2 2\\n')
            pathlib.Path(seed + '.mmn').write_text('raw\\n')
            pathlib.Path(seed + '.eig').write_text('raw\\n')
            pathlib.Path(seed + '-hyb.amn').write_text('hyb\\n2 2 2\\n')
            pathlib.Path(seed + '-hyb.mmn').write_text('hyb\\n')
            pathlib.Path(seed + '-hyb.eig').write_text('hyb\\n')
            pathlib.Path(seed + '-hyb.anchors').write_text('hyb\\n')
            pathlib.Path(seed + '-hyb.dis').write_text(
                '! wien2wannier disentanglement: stopped at the sweep cap (2000 sweeps, '
                'Omega_I 49.8 Ang^2): POOR CONVERGENCE: still descending\\n'
                '! wien2wannier subspace continuity:  7.200E-04 (smallest singular value '
                'of the hand-off overlaps), at k 257 -> 193; WARNING: below 1e-2\\n')
            """,
        )
        result = self.run_cli(*self.common_args("--quiet"))
        self.assertEqual(result.returncode, 0, result.stderr)
        handoff_win = (self.out / "model-hyb.win").read_text().splitlines()
        self.assertTrue(handoff_win[0].startswith("! wien2wannier disentanglement: stopped at the sweep cap"))
        self.assertTrue(handoff_win[1].startswith("! wien2wannier subspace continuity:"))
        self.assertIn("num_wann = 2", handoff_win)
        self.assertIn("POOR CONVERGENCE", result.stderr)
        self.assertIn("WARNING: below 1e-2", result.stderr)

    def test_overlaps_emits_validated_raw_products_without_hybrid_or_localize(self) -> None:
        result = self.run_cli(*self.common_args("--stage", "overlaps"))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((self.out / "model.inwf").read_text().splitlines()[0], "AUTO")
        for suffix in ("amn", "mmn", "eig"):
            self.assertTrue((self.out / f"model.{suffix}").is_file())
        self.assertFalse((self.out / "model-hyb.amn").exists())
        self.assertEqual([call["exe"] for call in self.calls()], ["w2w"])

    def test_overlaps_rejects_a_zero_exit_backend_without_raw_products(self) -> None:
        self._executable("w2w", "raise SystemExit(0)\n")
        result = self.run_cli(*self.common_args("--stage", "overlaps"))
        self.assertEqual(result.returncode, 2)
        self.assertIn("raw overlap", result.stderr)

    def test_nonempty_native_error_file_blocks_handoff_publication(self) -> None:
        self._executable(
            "w2w",
            """
            import pathlib, sys
            seed = pathlib.Path(sys.argv[-1]).stem
            pathlib.Path(seed + '.error').write_text('native STOP diagnostic\\n')
            pathlib.Path(seed + '.amn').write_text('raw\\n3 2 2\\n')
            pathlib.Path(seed + '.mmn').write_text('raw\\n')
            pathlib.Path(seed + '.eig').write_text('raw\\n')
            pathlib.Path(seed + '-hyb.amn').write_text('hyb\\n2 2 2\\n')
            pathlib.Path(seed + '-hyb.mmn').write_text('hyb\\n')
            pathlib.Path(seed + '-hyb.eig').write_text('hyb\\n')
            """,
        )
        result = self.run_cli(*self.common_args())
        self.assertEqual(result.returncode, 2)
        self.assertIn("native error file", result.stderr)
        self.assertFalse((self.out / "model-hyb.win").exists())

    def test_stale_handoff_is_never_accepted_as_a_new_backend_result(self) -> None:
        self.out.mkdir()
        for suffix in ("amn", "mmn", "eig", "win"):
            (self.out / f"model-hyb.{suffix}").write_text("stale\n")
        self._executable("w2w", "raise SystemExit(0)\n")
        result = self.run_cli(*self.common_args())
        self.assertEqual(result.returncode, 2)
        self.assertIn("pre-existing scientific output", result.stderr)
        self.assertEqual(self.calls(), [])

    def test_force_removes_stale_handoff_before_running_backend(self) -> None:
        self.out.mkdir()
        for suffix in ("amn", "mmn", "eig", "win"):
            (self.out / f"model-hyb.{suffix}").write_text("stale\n")
        self._executable("w2w", "raise SystemExit(0)\n")
        result = self.run_cli(*self.common_args("--force"))
        self.assertEqual(result.returncode, 2)
        self.assertIn("missing or empty", result.stderr)
        for suffix in ("amn", "mmn", "eig", "win"):
            self.assertFalse((self.out / f"model-hyb.{suffix}").exists())

    def test_ambiguous_collinear_input_requires_an_explicit_channel(self) -> None:
        for suffix in ("vectorup", "vectordn"):
            (self.case_dir / f"Si.{suffix}").write_text(suffix + "\n")
        result = self.run_cli(*self.common_args("--check"))
        self.assertEqual(result.returncode, 2)
        self.assertIn("select --spin up or --spin dn", result.stderr)
        self.assertFalse(self.out.exists())
        self.assertEqual(self.calls(), [])

    def test_empty_soc_and_complex_placeholders_do_not_change_auto_detection(self) -> None:
        (self.case_dir / "Si.inso").write_text("")
        (self.case_dir / "Si.in1c").write_text("")
        result = self.run_cli(*self.common_args("--check"))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("spin: none", result.stdout)
        self.assertIn("vector format: real", result.stdout)

    def test_mpi_dry_run_uses_mpi_binary_and_launcher_without_computing(self) -> None:
        self._executable("w2w_mpi", "raise SystemExit('must not run during dry-run')\n")
        self._executable("mpiexec", "raise SystemExit('must not run during dry-run')\n")
        result = self.run_cli(
            *self.common_args(
                "--dry-run", "--parallel", "mpi", "--mpi-ranks", "2",
                "--mpi-launcher", str(self.bin_dir / "mpiexec"),
            )
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertRegex(result.stdout, r"mpiexec\s+-n\s+2\s+.*w2w_mpi\s+model\.def")
        self.assertFalse(self.out.exists())
        self.assertEqual(self.calls(), [])

    def test_dry_run_reports_resolved_scientific_and_thread_controls(self) -> None:
        result = self.run_cli(
            *self.common_args(
                "--dry-run", "--parallel", "threads", "--threads", "2",
                "--solver", "cg", "--hybrid-option", "gauge=ptg",
                "--hybrid-option", "edges=strict",
            ),
            env={"W2W_SEEDLOSS": "abel"},
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("solver=cg", result.stdout)
        self.assertIn("gauge=ptg", result.stdout)
        self.assertIn("edges=strict", result.stdout)
        self.assertIn("seedloss=abel", result.stdout)
        self.assertIn("OMP=2, BLAS=2", result.stdout)
        self.assertFalse(self.out.exists())

    def _load_module(self):
        loader = importlib.machinery.SourceFileLoader("w2w_cli_alias_module", str(CLI))
        spec = importlib.util.spec_from_loader(loader.name, loader)
        assert spec is not None
        module = importlib.util.module_from_spec(spec)
        sys.modules[loader.name] = module
        loader.exec_module(module)
        return module

    def test_complex_mpi_backend_accepts_the_wien2k_name(self) -> None:
        # A WIEN2k installation names the complex MPI build w2wc_mpi.
        self._executable("w2wc_mpi", "raise SystemExit('must not run during dry-run')\n")
        self._executable("mpiexec", "raise SystemExit('must not run during dry-run')\n")
        wienroot = self.root / "wienroot"
        wienroot.mkdir()
        result = self.run_cli(
            *self.common_args(
                "--dry-run", "--parallel", "mpi", "--mpi-ranks", "2",
                "--vector-format", "complex",
                "--mpi-launcher", str(self.bin_dir / "mpiexec"),
            ),
            env={"WIENROOT": str(wienroot)},
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertRegex(result.stdout, r"mpiexec\s+-n\s+2\s+.*w2wc_mpi\s+model\.def")
        self.assertEqual(self.calls(), [])

    def test_backend_alias_search_keeps_directory_priority(self) -> None:
        module = self._load_module()
        first = self.root / "first"
        second = self.root / "second"
        for directory in (first, second):
            directory.mkdir()
        def put(directory: Path, name: str) -> None:
            path = directory / name
            path.write_text("#!/bin/sh\n")
            path.chmod(0o755)
        empty = self.root / "empty"
        empty.mkdir()
        old_path = os.environ.get("PATH", "")
        os.environ["PATH"] = str(empty)  # keep a real WIEN2k on PATH out of the search
        try:
            put(first, "w2wc_mpi")
            put(second, "w2w_mpic")
            found = module.locate_executable("w2w_mpic", [first, second])
            self.assertEqual(found, ((first / "w2wc_mpi").resolve(), "w2wc_mpi"))
            put(first, "w2w_mpic")
            found = module.locate_executable("w2w_mpic", [first, second])
            self.assertEqual(found, ((first / "w2w_mpic").resolve(), "w2w_mpic"))
            # the alias is one-way: other names take no alternative
            self.assertIsNone(module.locate_executable("w2wc_mpi", [second]))
            with self.assertRaisesRegex(module.CliError, "w2w_mpic or w2wc_mpi"):
                module.require_executable("w2w_mpic", "scientific backend", [empty])
        finally:
            os.environ["PATH"] = old_path

    def test_mpi_ranks_receive_the_thread_budget_in_omp_and_blas(self) -> None:
        source = (self.bin_dir / "w2w").read_text().split("\n", 1)[1]
        self._executable("w2w_mpi", source)
        self._open_mpi_standin()
        result = self.run_cli(
            *self.common_args(
                "--k-grid", "2", "1", "1", "--parallel", "mpi", "--mpi-ranks", "2",
                "--threads", "2", "--mpi-launcher", str(self.bin_dir / "mpiexec"),
                "--allow-oversubscribe",  # 2 x 2 threads must not depend on the host size
            )
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        calls = self.calls()
        self.assertEqual(calls[0]["exe"], "w2w_mpi")
        self.assertEqual(calls[0]["omp"], "2")
        self.assertEqual(calls[0]["openblas"], "2")
        # Open MPI would bind each of the two ranks to one core; the launcher unbinds them.
        launched = (self.root / "launcher_argv.txt").read_text().split()
        self.assertEqual(launched[:4], ["-n", "2", "--bind-to", "none"])

    def _open_mpi_standin(self, version: str = "mpiexec (Open MPI) 5.0.9") -> None:
        # Answers --version like Open MPI, records its argv, and runs the backend once
        # with the inherited environment.
        self._executable(
            "mpiexec",
            "import os, pathlib, sys\n"
            "args = sys.argv[1:]\n"
            "if args == ['--version']:\n"
            f"    print({version!r})\n"
            "    sys.exit(0)\n"
            "log = pathlib.Path(os.environ['W2W_TEST_LOG']).with_name('launcher_argv.txt')\n"
            "log.write_text(' '.join(args))\n"
            "assert args[0] == '-n'\n"
            "rest = args[2:]\n"
            "if rest[:2] == ['--bind-to', 'none']:\n"
            "    rest = rest[2:]\n"
            "os.execv(rest[0], rest)\n",
        )

    def _mpi_dry_run(self, *extra: str) -> subprocess.CompletedProcess[str]:
        self._executable("w2w_mpi", "raise SystemExit('must not run during dry-run')\n")
        return self.run_cli(
            *self.common_args(
                "--dry-run", "--parallel", "mpi", "--mpi-ranks", "2",
                "--mpi-launcher", str(self.bin_dir / "mpiexec"), "--allow-oversubscribe",
                *extra,
            )
        )

    def test_mpi_bind_auto_unbinds_open_mpi_ranks_with_threads(self) -> None:
        self._open_mpi_standin()
        result = self._mpi_dry_run("--threads", "2")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertRegex(result.stdout, r"mpiexec\s+-n\s+2\s+--bind-to\s+none\s+\S*w2w_mpi\s+model\.def")
        self.assertIn("MPI binding: --bind-to none (Open MPI 5.0.9, 2 threads per rank)", result.stdout)
        self.assertFalse(self.out.exists())

    def test_mpi_bind_auto_leaves_other_launchers_alone(self) -> None:
        self._open_mpi_standin(version="HYDRA build details:")
        result = self._mpi_dry_run("--threads", "2")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertRegex(result.stdout, r"mpiexec\s+-n\s+2\s+\S*w2w_mpi\s+model\.def")
        self.assertIn("MPI binding: left to the launcher (not identified as Open MPI", result.stdout)

    def test_mpi_bind_skips_the_probe_when_it_cannot_matter(self) -> None:
        # A launcher that must not be executed at all: one thread per rank, or an explicit choice.
        self._executable(
            "mpiexec",
            "import os, pathlib\n"
            "pathlib.Path(os.environ['W2W_TEST_LOG']).with_name('probed').write_text('x')\n"
            "raise SystemExit('must not run')\n",
        )
        for extra, expected in (
            ((), "left to the launcher (one thread per rank)"),
            (("--threads", "2", "--mpi-bind", "launcher"), "left to the launcher (--mpi-bind launcher)"),
            (("--threads", "2", "--mpi-bind", "none"), "--bind-to none (--mpi-bind none)"),
        ):
            with self.subTest(extra=extra):
                result = self._mpi_dry_run(*extra)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn(f"MPI binding: {expected}", result.stdout)
                has_option = re.search(r"mpiexec\s+-n\s+2\s+--bind-to\s+none\s", result.stdout)
                self.assertEqual(bool(has_option), extra[-1:] == ("none",))
                self.assertFalse((self.root / "probed").exists())

    def test_slurm_per_task_allocation_accepts_one_thread_on_each_rank(self) -> None:
        self._executable("w2w_mpi", "raise SystemExit('must not run during dry-run')\n")
        self._executable("mpiexec", "raise SystemExit('must not run during dry-run')\n")
        result = self.run_cli(
            *self.common_args(
                "--dry-run", "--parallel", "mpi", "--mpi-ranks", "2", "--threads", "1",
                "--mpi-launcher", str(self.bin_dir / "mpiexec"),
            ),
            env={"SLURM_CPUS_PER_TASK": "1", "SLURM_NTASKS": "2", "SLURM_CPUS_ON_NODE": "2"},
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("2 rank(s), 1 thread(s)/rank", result.stdout)
        self.assertFalse(self.out.exists())

    def test_enclosing_mpi_world_is_rejected_before_writes(self) -> None:
        for variable in (
            "OMPI_COMM_WORLD_SIZE", "PMI_SIZE", "PMIX_SIZE", "MV2_COMM_WORLD_SIZE"
        ):
            with self.subTest(variable=variable):
                result = self.run_cli(*self.common_args("--check"), env={variable: "2"})
                self.assertEqual(result.returncode, 2)
                self.assertIn("MPI world", result.stderr)
                self.assertFalse(self.out.exists())

    def test_check_rejects_invalid_scientific_controls(self) -> None:
        bad_cases = (
            (("--band-max", "-5"), "band-max"),
            (("--hybrid-option", "mystery=1"), "HYBR keyword"),
            (("--projection-threshold", "2"), "projection-threshold"),
            (("--pool-factor", "nan"), "pool-factor"),
            (("--froz-emax", "inf"), "froz-emax"),
        )
        for extra, message in bad_cases:
            with self.subTest(extra=extra):
                result = self.run_cli(*self.common_args("--check", *extra))
                self.assertEqual(result.returncode, 2)
                self.assertIn(message, result.stderr)
                self.assertFalse(self.out.exists())

    def test_window_route_reports_continuity_after_localization(self) -> None:
        result = self.run_cli(*self.common_args(
            "--method", "window", "--num-wann", "2",
            "--wannier90", str(self.bin_dir / "wannier90.x")))
        self.assertEqual(result.returncode, 0, result.stderr)
        calls = [c for c in self.calls() if c["exe"] in ("wannier90.x", "w2w_hybrid")]
        localize = [i for i, c in enumerate(calls)
                    if c["exe"] == "wannier90.x" and c["argv"] == ["model-hyb-win"]]
        check = [i for i, c in enumerate(calls)
                 if c["exe"] == "w2w_hybrid" and c["argv"] == ["--continuity-chk", "model-hyb-win"]]
        self.assertTrue(localize, calls)
        self.assertTrue(check, calls)
        self.assertGreater(check[0], localize[0])
        self.assertIn("subspace continuity", result.stdout)

    def test_check_localize_requires_a_complete_existing_handoff(self) -> None:
        result = self.run_cli(
            *self.common_args(
                "--stage", "localize", "--wannier90", str(self.bin_dir / "wannier90.x"),
                "--check",
            )
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("existing handoff", result.stderr)
        self.assertFalse(self.out.exists())

    def test_dotted_seed_preserves_the_full_basename(self) -> None:
        result = self.run_cli(
            "--input", str(self.case_dir), "--seed", "model.v2", "--output-dir",
            str(self.out), "--backend-dir", str(self.bin_dir),
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue((self.out / "model.v2-hyb.amn").is_file())
        self.assertTrue((self.out / "model.v2-hyb.win").is_file())

    def test_fortran_def_escapes_apostrophes_in_input_paths(self) -> None:
        quoted = self.root / "O'Brien"
        self.case_dir.rename(quoted)
        self.case_dir = quoted
        result = self.run_cli(*self.common_args("--case", "Si"))
        self.assertEqual(result.returncode, 0, result.stderr)
        deffile = (self.out / "model.def").read_text()
        self.assertIn("O''Brien/Si.vector", deffile)

    def test_bridge_forwards_pool_factor(self) -> None:
        bridge = self.root / "bridge-pool"
        bridge.mkdir()
        for suffix in ("amn", "mmn", "eig", "nnkp"):
            (bridge / f"raw.{suffix}").write_text(suffix + "\n")
        (bridge / "raw.win").write_text("num_bands=5\nnum_wann=4\n")
        result = self.run_cli(
            "--input", str(bridge), "--input-kind", "bridge", "--bridge-seed", "raw",
            "--seed", "model", "--output-dir", str(self.out), "--backend-dir",
            str(self.bin_dir), "--num-wann", "3", "--pool-factor", "1.75",
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--pool-factor", self.calls()[0]["argv"])
        index = self.calls()[0]["argv"].index("--pool-factor")
        self.assertEqual(self.calls()[0]["argv"][index + 1], "1.75")

    def test_bridge_thread_mode_avoids_nested_blas_threads(self) -> None:
        bridge = self.root / "bridge-threads"
        bridge.mkdir()
        for suffix in ("amn", "mmn", "eig", "nnkp"):
            (bridge / f"raw.{suffix}").write_text(suffix + "\n")
        (bridge / "raw.win").write_text("num_bands=5\nnum_wann=4\n")
        source = (self.bin_dir / "w2w_hybrid").read_text().split("\n", 1)[1]
        self._executable("w2w_hybrid_acc", source)
        result = self.run_cli(
            "--input", str(bridge), "--input-kind", "bridge", "--bridge-seed", "raw",
            "--seed", "model", "--output-dir", str(self.out), "--backend-dir",
            str(self.bin_dir), "--num-wann", "3", "--parallel", "threads", "--threads", "2",
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.calls()[0]["exe"], "w2w_hybrid_acc")
        self.assertEqual(self.calls()[0]["omp"], "2")
        self.assertEqual(self.calls()[0]["openblas"], "1")

    def test_bridge_preserves_binary_sidecars_and_accepts_bin_only_payloads(self) -> None:
        bridge = self.root / "bridge-sidecars"
        bridge.mkdir()
        for suffix in ("eig", "nnkp", "win"):
            (bridge / f"raw.{suffix}").write_text(suffix + "\n")
        for suffix in ("amn.bin", "mmn.bin"):
            (bridge / f"raw.{suffix}").write_bytes(b"native-sidecar\n")
        result = self.run_cli(
            "--input", str(bridge), "--input-kind", "bridge", "--bridge-seed", "raw",
            "--seed", "model", "--output-dir", str(self.out), "--backend-dir",
            str(self.bin_dir), "--num-wann", "3", "--stage", "prepare",
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        for suffix in ("amn.bin", "mmn.bin"):
            destination = self.out / f"raw.{suffix}"
            self.assertTrue(destination.is_symlink())
            self.assertEqual(destination.read_bytes(), b"native-sidecar\n")

    def test_bridge_rejects_a_stale_local_sidecar_missing_from_source(self) -> None:
        bridge = self.root / "bridge-stale-sidecar"
        bridge.mkdir()
        for suffix in ("amn", "mmn", "eig", "nnkp", "win"):
            (bridge / f"raw.{suffix}").write_text(suffix + "\n")
        self.out.mkdir()
        (self.out / "raw.amn.bin").write_bytes(b"stale\n")
        result = self.run_cli(
            "--input", str(bridge), "--input-kind", "bridge", "--bridge-seed", "raw",
            "--seed", "model", "--output-dir", str(self.out), "--backend-dir",
            str(self.bin_dir), "--num-wann", "3", "--stage", "prepare",
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("stale bridge cache", result.stderr)
        self.assertEqual((self.out / "raw.amn.bin").read_bytes(), b"stale\n")

    def test_cli_options_override_inherited_backend_environment(self) -> None:
        result = self.run_cli(
            *self.common_args("--hybrid-option", "edges=strict", "--hybrid-option", "seedk=opt16"),
            env={"W2W_EDGES": "warn", "W2W_SEEDK": "gamma", "W2W_SEEDLOSS": "legacy"},
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        call_env = (self.out / "model.inwf").read_text()
        self.assertIn("edges = strict", call_env)
        self.assertIn("seedk = opt16", call_env)
        self.assertIsNone(self.calls()[0]["edges"])
        self.assertIsNone(self.calls()[0]["seedk"])
        self.assertEqual(self.calls()[0]["seedloss"], "legacy")

    def test_relative_backend_input_environment_is_resolved_before_chdir(self) -> None:
        structure = self.root / "environment.struct"
        structure.write_text("environment structure\n")
        result = self.run_cli(
            *self.common_args(), env={"W2W_STRUCT": structure.name}
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.calls()[0]["struct_env"], str(structure.resolve()))

    def test_bridge_input_output_overlap_is_rejected_even_with_force(self) -> None:
        bridge = self.root / "bridge-source"
        bridge.mkdir()
        original: dict[str, bytes] = {}
        for suffix in ("amn", "mmn", "eig", "nnkp", "win"):
            path = bridge / f"raw.{suffix}"
            path.write_text(f"original-{suffix}\n")
            original[suffix] = path.read_bytes()
        alias = self.root / "bridge-alias"
        alias.symlink_to(bridge, target_is_directory=True)
        for input_path, output_path in ((bridge, bridge), (alias, bridge)):
            with self.subTest(input=input_path):
                result = self.run_cli(
                    "--input", str(input_path), "--input-kind", "bridge", "--bridge-seed", "raw",
                    "--seed", "model", "--output-dir", str(output_path), "--backend-dir",
                    str(self.bin_dir), "--num-wann", "3", "--stage", "prepare", "--force",
                )
                self.assertEqual(result.returncode, 2)
                self.assertIn("physically separate", result.stderr)
                for suffix, data in original.items():
                    self.assertEqual((bridge / f"raw.{suffix}").read_bytes(), data)

    def test_force_never_unlinks_a_hardlinked_bridge_source(self) -> None:
        bridge = self.root / "bridge-hardlink"
        bridge.mkdir()
        for suffix in ("amn", "mmn", "eig", "nnkp", "win"):
            (bridge / f"raw.{suffix}").write_text(f"source-{suffix}\n")
        self.out.mkdir()
        os.link(bridge / "raw.amn", self.out / "raw.amn")
        result = self.run_cli(
            "--input", str(bridge), "--input-kind", "bridge", "--bridge-seed", "raw",
            "--seed", "model", "--output-dir", str(self.out), "--backend-dir",
            str(self.bin_dir), "--num-wann", "3", "--stage", "prepare", "--force",
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((bridge / "raw.amn").read_text(), "source-amn\n")
        self.assertTrue(os.path.samefile(bridge / "raw.amn", self.out / "raw.amn"))

    def test_bridge_mode_runs_supported_standalone_route(self) -> None:
        bridge = self.root / "bridge"
        bridge.mkdir()
        for suffix in ("amn", "mmn", "eig", "nnkp"):
            (bridge / f"raw.{suffix}").write_text(suffix + "\n")
        (bridge / "raw.win").write_text("num_bands=5\nnum_wann=4\ndis_num_iter=9\n")
        result = self.run_cli(
            "--input", str(bridge), "--input-kind", "bridge", "--bridge-seed", "raw",
            "--seed", "model", "--output-dir", str(self.out), "--backend-dir",
            str(self.bin_dir), "--num-wann", "3", "--bridge-arg=--gauge",
            "--bridge-arg=auto8",
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([call["exe"] for call in self.calls()], ["w2w_hybrid"])
        self.assertEqual(
            self.calls()[0]["argv"],
            ["--seed", "raw", "--nw", "3", "--out", "model-hyb",
             "--froz-emax", "1.0", "--pool-factor", "1.5", "--solver", "cg",
             "--gauge", "auto8"],
        )
        self.assertTrue((self.out / "raw.amn").is_symlink())
        self.assertIn("num_wann = 3", (self.out / "model-hyb.win").read_text())

    def test_bridge_adapts_to_legacy_positional_backend(self) -> None:
        bridge = self.root / "legacy-bridge"
        bridge.mkdir()
        for suffix in ("amn", "mmn", "eig", "nnkp"):
            (bridge / f"raw.{suffix}").write_text(suffix + "\n")
        (bridge / "raw.win").write_text("num_bands=5\nnum_wann=4\n")
        self._executable(
            "w2w_hybrid",
            """
            # usage: w2w_hybrid <seed>; pool factor =; arg 20 (solver)
            import json, os, pathlib, sys
            with pathlib.Path(os.environ['W2W_TEST_LOG']).open('a') as stream:
                stream.write(json.dumps({'exe': pathlib.Path(sys.argv[0]).name,
                                         'argv': sys.argv[1:]}) + '\\n')
            out = sys.argv[3]
            pathlib.Path(out + '.amn').write_text(
                'legacy bridge handoff\\n           3           2           3\\n')
            pathlib.Path(out + '.mmn').write_text('legacy mmn\\n')
            pathlib.Path(out + '.eig').write_text('legacy eig\\n')
            """,
        )
        result = self.run_cli(
            "--input", str(bridge), "--input-kind", "bridge", "--bridge-seed", "raw",
            "--seed", "model", "--output-dir", str(self.out), "--backend-dir",
            str(self.bin_dir), "--num-wann", "3", "--pool-factor", "1.75",
            "--solver", "cg", "--colored-gs",
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(
            self.calls()[0]["argv"],
            ["raw", "3", "model-hyb", "occ99", "auto", "1.0", "0.0",
             "2000", "1e-10", "10", "1", "0", "1.0", "scdm", "0", "0",
             "none", "none", "1.75", "cg"],
        )
        self.assertTrue((self.out / "model-hyb.amn").is_file())

        rejected = self.run_cli(
            "--input", str(bridge), "--input-kind", "bridge", "--bridge-seed", "raw",
            "--seed", "second", "--output-dir", str(self.root / "second"), "--backend-dir",
            str(self.bin_dir), "--num-wann", "3", "--bridge-arg=--gauge",
        )
        self.assertEqual(rejected.returncode, 2)
        self.assertIn("positional bridge backend", rejected.stderr)

    def test_bridge_rejects_an_unknown_backend_interface_before_writes(self) -> None:
        bridge = self.root / "unknown-bridge"
        bridge.mkdir()
        for suffix in ("amn", "mmn", "eig", "nnkp", "win"):
            (bridge / f"raw.{suffix}").write_text(suffix + "\n")
        self._executable("w2w_hybrid", "raise SystemExit(0)\n")
        result = self.run_cli(
            "--input", str(bridge), "--input-kind", "bridge", "--bridge-seed", "raw",
            "--seed", "model", "--output-dir", str(self.out), "--backend-dir",
            str(self.bin_dir), "--num-wann", "3", "--check",
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("unrecognized bridge backend interface", result.stderr)
        self.assertFalse(self.out.exists())
        self.assertEqual(self.calls(), [])

    def test_bridge_mpi_is_rejected_before_output_creation(self) -> None:
        bridge = self.root / "bridge"
        bridge.mkdir()
        for suffix in ("amn", "mmn", "eig", "nnkp", "win"):
            (bridge / f"raw.{suffix}").write_text(suffix + "\n")
        result = self.run_cli(
            "--input", str(bridge), "--input-kind", "bridge", "--bridge-seed", "raw",
            "--seed", "model", "--output-dir", str(self.out), "--backend-dir",
            str(self.bin_dir), "--num-wann", "3", "--parallel", "mpi", "--mpi-ranks", "2",
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("bridge backend does not support MPI", result.stderr)
        self.assertFalse(self.out.exists())
        self.assertEqual(self.calls(), [])

    def test_launcher_only_install_exposes_public_command(self) -> None:
        install_dir = self.root / "install path"
        result = subprocess.run(
            ["make", "-s", "install-cli", f"target-dir={install_dir}"],
            cwd=ROOT,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        installed = install_dir / "wien2wannier"
        self.assertTrue(installed.is_file())
        self.assertTrue(os.access(installed, os.X_OK))
        help_run = subprocess.run(
            [str(installed), "--help"], text=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, check=False,
        )
        self.assertEqual(help_run.returncode, 0, help_run.stderr)


if __name__ == "__main__":
    unittest.main()
