#!/usr/bin/env bash
#SBATCH --job-name=halsynteny-palaeo
#SBATCH --time=24:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --output=halsynteny-%j.out
set -euo pipefail
module load cactus/3.3.0 2>/dev/null || true
export PATH="${HAL_ENV:-$PWD/.hal-env/bin}:$PATH"
exec bash experiments/halsynteny-baseline/commands.sh
