#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
../monaco |& tee one_phase_run.log
