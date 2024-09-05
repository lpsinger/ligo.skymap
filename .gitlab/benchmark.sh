#!/bin/bash

set -e

# From https://docs.gitlab.com/ee/ci/jobs/index.html#use-a-script-to-improve-display-of-collapsible-sections
function section_start () {
  local section_title="${1}"
  local section_description="${2:-$section_title}"

  echo -e "section_start:`date +%s`:${section_title}[collapsed=true]\r\e[0K${section_description}"
}

# From https://docs.gitlab.com/ee/ci/jobs/index.html#use-a-script-to-improve-display-of-collapsible-sections
function section_end () {
  local section_title="${1}"

  echo -e "section_end:`date +%s`:${section_title}\r\e[0K"
}

section_start "clone" "Cloning benchmark repo"
git clone https://git.ligo.org/leo-singer/ligo-skymap-benchmark.git
section_end "clone"

section_start "venv" "Creating Python virtual environment"
python3.12 -m venv env
source env/bin/activate
pip install pandas
section_end "venv"

if [[ -n "${1}" ]]; then
section_start "install_binary" "Installing ligo.skymap from binary package"
pip install "${1}"
section_end "install_binary"
section_start "benchmark_binary" "Running benchmark on binary package"
pushd ligo-skymap-benchmark
./benchmark.py
./plot.py
mv benchmark.svg ../
popd
section_end "benchmark_binary"
fi

section_start "install_source" "Installing ligo.skymap from source"
source /opt/intel/oneapi/compiler/latest/env/vars.sh --include-intel-llvm
source /opt/intel/oneapi/advisor/latest/advisor-vars.sh
CC=icc CXX=icpc AR=xiar LDSHARED="icc -shared" \
    CFLAGS="-O3 -g -axCORE-AVX512,CORE-AVX2,AVX,SSE4.2" \
    LIGO_SKYMAP_USE_ITTNOTIFY=1 pip install -ve .
section_end "install_source"

section_start "advisor" "Profiling with Intel Advisor"
pushd ligo-skymap-benchmark
advisor --collect=roofline --start-paused --module-filter-mode=include \
    --module-filter=core.abi3.so,libgsl.so.25,libm.so.6,libc.so.6 \
    --project-dir=advisor -- ./benchmark.py
advisor --report=roofline --project-dir=advisor --report-output=../roofline.html
popd
section_end "advisor"

section_start "perf" "Profiling with perf"
git clone https://github.com/brendangregg/FlameGraph
pushd ligo-skymap-benchmark
perf record -g -- ./benchmark.py
perf script | ../FlameGraph/stackcollapse-perf.pl > out.perf-folded
../FlameGraph/flamegraph.pl out.perf-folded > ../perf.svg
popd
section_end "perf"
