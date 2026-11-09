#!/usr/bin/env bash
# Build and run CommonSense inside the docker/Dockerfile environment.
#
#   docker/run.sh image                  build the Docker image (once)
#   docker/run.sh build [Debug|Release]  configure + compile into the build volume
#   docker/run.sh sims                   compile sc_study/sc_sim and sc_study/bidir_sim (Linux binaries)
#   docker/run.sh two_party CONFIG.json  run a binary (two_party, two_party_iblt, two_party_i64)
#   docker/run.sh shell                  interactive shell in the container
#   docker/run.sh exec CMD...            run any command in the container (cwd = repo root)
#
# The repo is mounted at /work, so configs, results and source edits are shared
# with the host. The CMake build tree lives in a Docker volume at /build: the
# third-party build steps need a case-sensitive filesystem, which the macOS
# default is not. Binaries are therefore at /build/bin inside the container.
set -euo pipefail

IMAGE=commonsense-dev
VOLUME=commonsense-build
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

run() {
  local tty=()
  [ -t 0 ] && [ -t 1 ] && tty=(-it)
  docker run --rm ${tty[@]+"${tty[@]}"} -v "$REPO":/work -v "$VOLUME":/build "$IMAGE" "$@"
}

cmd="${1:-}"
shift || true
case "$cmd" in
  image)
    docker build -t "$IMAGE" "$REPO/docker" ;;
  build)
    run bash -c "cmake -S . -B /build -DCMAKE_BUILD_TYPE=${1:-Release} && cmake --build /build -j\$(nproc)" ;;
  sims)
    run bash -c "mkdir -p /build/sc_study && \
      g++ -O3 -std=c++17 -o /build/sc_study/sc_sim sc_study/sc_sim.cpp -lpthread && \
      g++ -O3 -std=c++17 -o /build/sc_study/bidir_sim sc_study/bidir_sim.cpp -lpthread && \
      echo 'built /build/sc_study/{sc_sim,bidir_sim}'" ;;
  two_party|two_party_iblt|two_party_i64)
    run "/build/bin/$cmd" "$@" ;;
  shell)
    run bash ;;
  exec)
    run "$@" ;;
  *)
    sed -n '2,13p' "$0" | sed 's/^# \{0,1\}//'
    exit 1 ;;
esac
