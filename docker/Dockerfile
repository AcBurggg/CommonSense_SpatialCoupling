# Reproducible Linux build environment for CommonSense (two_party, two_party_iblt,
# two_party_i64) and the sc_study simulators.
FROM gcc:14.2.0

RUN apt-get update \
 && apt-get install -y --no-install-recommends cmake python3 python3-numpy python3-scipy \
 && rm -rf /var/lib/apt/lists/*

# The repo is bind-mounted here; the CMake build tree (including the third-party
# libraries it downloads at configure time) lives in a named volume at /build.
WORKDIR /work
ENV CXX=g++ CC=gcc
