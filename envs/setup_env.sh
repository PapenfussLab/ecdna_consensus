#!/usr/bin/env bash
set -euo pipefail

PREFIX="${1:?usage: bash envs/setup_env.sh /path/to/micromamba_env/ecdna_consensus}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ -d "$PREFIX/conda-meta" ]; then
  micromamba env update -y -p "$PREFIX" -f "$HERE/environment_consensus.yml"
else
  micromamba env create -y -p "$PREFIX" -f "$HERE/environment_consensus.yml"
fi

micromamba run -p "$PREFIX" Rscript "$HERE/install_mskilab.R"
