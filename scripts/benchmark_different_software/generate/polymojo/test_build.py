"""Setup/cache tests; no Mojo installation or package downloads required."""

from concurrent.futures import ThreadPoolExecutor
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import build


class BuildTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory(prefix="poly mojo ")
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.prefix = self.root / "environment"
        self.packages = self.prefix / "lib" / "mojo"
        self.packages.mkdir(parents=True)
        for name in build.PACKAGES:
            (self.packages / name).write_text(name)
        compiler = self.prefix / "bin" / "mojo"
        compiler.parent.mkdir()
        compiler.write_text(
            f"#!{sys.executable}\n"
            "import os, pathlib, sys\n"
            "if sys.argv[1] == '--version':\n"
            "    print(os.environ.get('TEST_VERSION', 'Mojo 1.2.0 (test)'))\n"
            "else:\n"
            "    with open(os.environ['TEST_BUILD_LOG'], 'a') as log: log.write('build\\n')\n"
            "    output = pathlib.Path(sys.argv[sys.argv.index('-o') + 1])\n"
            "    output.write_text('test executable')\n"
            "    output.chmod(0o755)\n"
            "    if os.environ.get('TEST_FAIL'): sys.exit(3)\n"
        )
        compiler.chmod(0o755)
        self.log = self.root / "builds"
        self.source = self.root / "driver.mojo"
        self.source.write_text("test source")
        self.env = {
            "MOJO_COMMAND": str(compiler),
            "TEST_BUILD_LOG": str(self.log),
            "PATH": os.defpath,
        }

    def prepare(self):
        return build.build_driver(
            source=self.source, cache_dir=self.root / "cache", environ=self.env
        )

    def count(self):
        return len(self.log.read_text().splitlines()) if self.log.exists() else 0

    def test_installed_prefix_and_repeat(self):
        first = self.prepare()
        self.assertTrue(os.access(first, os.X_OK))
        self.assertEqual(self.prepare(), first)
        self.assertEqual(self.count(), 1)

    def test_relative_cli_path_reuses_default_cache(self):
        checkout = self.root / "checkout"
        folder = checkout / "generate" / "polymojo"
        environment = folder / "environment"
        environment.mkdir(parents=True)
        helper = folder / "build.py"
        helper.write_bytes(Path(build.__file__).read_bytes())
        (folder / "driver.mojo").write_bytes(self.source.read_bytes())
        first = subprocess.check_output(
            [sys.executable, str(helper)], cwd=checkout, env=self.env, text=True
        ).strip()
        second = subprocess.check_output(
            [sys.executable, "../build.py"], cwd=environment, env=self.env, text=True
        ).strip()
        self.assertEqual(first, second)
        self.assertEqual(Path(first).parents[2], checkout / "results")
        self.assertEqual(self.count(), 1)

    def test_concurrent_setup_compiles_once(self):
        with ThreadPoolExecutor(max_workers=2) as pool:
            binaries = list(pool.map(lambda _: self.prepare(), range(2)))
        self.assertEqual(binaries[0], binaries[1])
        self.assertEqual(self.count(), 1)

    def test_library_and_source_changes_invalidate_cache(self):
        first = self.prepare()
        (self.packages / "polynomials.mojoc").write_text("updated library")
        second = self.prepare()
        self.source.write_text("updated bridge")
        third = self.prepare()
        self.assertEqual(len({first, second, third}), 3)
        self.assertEqual(self.count(), 3)

    def test_explicit_packages_do_not_fall_back(self):
        self.env["POLY_MOJO_PATH"] = str(self.root / "missing")
        with self.assertRaisesRegex(build.SetupError, "Precompiled Poly Mojo packages"):
            self.prepare()
        self.assertEqual(self.count(), 0)

    def test_wrong_compiler_version(self):
        self.env["TEST_VERSION"] = "Mojo 1.1.0 (test)"
        with self.assertRaisesRegex(build.SetupError, "requires Mojo 1.2.0"):
            self.prepare()
        self.assertEqual(self.count(), 0)

    def test_missing_compiler(self):
        self.env["MOJO_COMMAND"] = str(self.root / "missing-compiler")
        with self.assertRaisesRegex(build.SetupError, "Mojo not found"):
            self.prepare()

    def test_failed_compilation_is_not_cached(self):
        self.env["TEST_FAIL"] = "1"
        with self.assertRaises(subprocess.CalledProcessError):
            self.prepare()
        self.assertEqual(list((self.root / "cache").glob("*/driver")), [])
        del self.env["TEST_FAIL"]
        self.assertTrue(self.prepare().is_file())
        self.assertEqual(self.count(), 2)


if __name__ == "__main__":
    unittest.main()
