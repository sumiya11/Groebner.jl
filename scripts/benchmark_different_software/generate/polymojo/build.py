"""Compile the benchmark bridge once against already installed Mojo packages."""

import fcntl
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import tempfile

SOURCE = Path(__file__).resolve().with_name("driver.mojo")
PACKAGES = ("apn_mojo.mojoc", "polynomials.mojoc")


class SetupError(RuntimeError):
    """An actionable error in the optional benchmark environment."""


def package_directory(compiler, environ):
    explicit = environ.get("POLY_MOJO_PATH")
    if explicit:
        candidates = [Path(explicit).expanduser().resolve()]
    else:
        candidates = []
        if environ.get("CONDA_PREFIX"):
            candidates.append(Path(environ["CONDA_PREFIX"]) / "lib" / "mojo")
        candidates.append(compiler.parent.parent / "lib" / "mojo")
    for folder in candidates:
        if all((folder / name).is_file() for name in PACKAGES):
            return folder.resolve()
    locations = ", ".join(map(str, candidates))
    raise SetupError(
        "Precompiled Poly Mojo packages not found in " + locations
        + "; activate the package environment or set POLY_MOJO_PATH."
    )


def build_driver(*, source=SOURCE, cache_dir=None, environ=None):
    environ = dict(os.environ if environ is None else environ)
    command = environ.get("MOJO_COMMAND", "mojo")
    executable = shutil.which(command, path=environ.get("PATH", os.defpath))
    if executable is None:
        raise SetupError("Mojo not found; activate its environment or set MOJO_COMMAND.")
    compiler = Path(executable).resolve()
    version = subprocess.check_output([str(compiler), "--version"], env=environ, text=True).strip()
    if not re.search(r"\bMojo 1\.2\.0(?:\s|$)", version):
        raise SetupError(f"This bridge requires Mojo 1.2.0; found {version}.")
    packages = package_directory(compiler, environ)
    source = Path(source).resolve()
    identity = {
        "compiler": str(compiler),
        "version": version,
        "modular_home": environ.get("MODULAR_HOME", ""),
        "platform": platform.system(),
        "machine": platform.machine(),
        "source": hashlib.sha256(source.read_bytes()).hexdigest(),
        "packages": {
            name: hashlib.sha256((packages / name).read_bytes()).hexdigest()
            for name in PACKAGES
        },
    }
    key = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()[:24]
    if cache_dir is None:
        cache_dir = SOURCE.parents[2] / "results" / "polymojo-build"
    folder = Path(cache_dir).resolve() / key
    folder.mkdir(parents=True, exist_ok=True)
    binary = folder / "driver"
    metadata = folder / "build.json"
    # Separate setup invocations may race; only one compiles each cache entry.
    with (folder / "build.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if binary.is_file() and os.access(binary, os.X_OK) and metadata.is_file():
            if json.loads(metadata.read_text()) == identity:
                return binary
        with tempfile.TemporaryDirectory(dir=folder) as temporary:
            output = Path(temporary) / "driver"
            subprocess.run(
                [str(compiler), "build", "-I", str(packages), str(source), "-o", str(output)],
                env=environ,
                stdout=sys.stderr,
                check=True,
            )
            if not output.is_file() or not os.access(output, os.X_OK):
                raise SetupError("Mojo did not produce an executable benchmark bridge.")
            record = Path(temporary) / "build.json"
            record.write_text(json.dumps(identity, indent=2) + "\n")
            output.replace(binary)
            record.replace(metadata)
    return binary


def main():
    try:
        print(build_driver())
    except (SetupError, OSError, subprocess.CalledProcessError) as error:
        print(f"Poly Mojo setup: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
