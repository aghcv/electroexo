#!/usr/bin/env bash
# A single executable for ELECTROEXO_PYTHON on Wahab's containerized Python.
set -eo pipefail

repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
environment="${ELECTROEXO_ENV:-$repo_dir/.venv}"
python_module="${ELECTROEXO_WAHAB_MODULE:-python3/2026.0-py312}"

# A container venv's Python symlink may not resolve on the host. Check the
# activation files here, then let crun execute Python inside the container.
if [[ ! -f "$environment/pyvenv.cfg" || ! -f "$environment/bin/activate" ]]; then
    printf 'ElectroExo environment is missing or incomplete: %s\n' "$environment" >&2
    printf 'Create the Python 3.12 environment using docs/wahab_integration.md.\n' >&2
    exit 2
fi
environment="$(cd -- "$environment" && pwd)"

# Site initialization scripts do not support nounset reliably.
if ! command -v module >/dev/null 2>&1; then
    if [[ -r /etc/profile.d/modules.sh ]]; then
        source /etc/profile.d/modules.sh
    else
        printf 'Wahab module initialization is unavailable; run on a Wahab node.\n' >&2
        exit 2
    fi
fi
module load container_env "$python_module"
set -u

if ! command -v crun >/dev/null 2>&1; then
    printf 'crun is unavailable after loading %s.\n' "$python_module" >&2
    exit 2
fi

exec crun -p "$environment" python "$@"
