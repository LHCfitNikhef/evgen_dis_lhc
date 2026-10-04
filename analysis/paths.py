"""Path configuration for the benchmark's Python analysis.

Values come from ``config.sh`` -- the single source -- read through
``tools/print_config.sh``, so this module never re-declares a default that
could drift out of step with the shell scripts and the Makefiles.

    from paths import SHERPA_RUNS, POWHEG_RES

Override the same way as everywhere else, e.g. ``PHYSICS24=/project/faser`` or
``BENCH_SHERPA_RUNS=/scratch/runs`` in the environment.
"""

import functools
import os
import subprocess

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@functools.lru_cache(maxsize=1)
def _config():
    """Read config.sh once, as KEY=VALUE lines."""
    env = dict(os.environ, BENCH_REPO=REPO)
    out = subprocess.run(
        [os.path.join(REPO, "tools", "print_config.sh")],
        check=True, capture_output=True, text=True, env=env,
    ).stdout
    return dict(line.split("=", 1) for line in out.splitlines() if "=" in line)


def get(name):
    """Return one configured path, or raise if config.sh does not define it."""
    try:
        return _config()[name]
    except KeyError:
        raise KeyError(
            f"{name} is not defined in config.sh (known: "
            f"{', '.join(sorted(_config()))})"
        ) from None


def __getattr__(name):          # module-level attribute access, PEP 562
    if name.startswith("_"):
        raise AttributeError(name)
    try:
        return get(name)
    except KeyError as exc:
        raise AttributeError(str(exc)) from None
