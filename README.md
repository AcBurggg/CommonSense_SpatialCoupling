# CommonSense: Efficient Set Intersection (SetX) Protocol Based on Compressed Sensing

## Prerequestites

You need g++-14 installed in your system. 

For other dependencies, the build script will automatically install them in correct places. 

## Build & Run

To build the artifacts for CDMA, run the following command: 
```bash
CXX=g++-14 cmake -S . -B build
cmake --build build -j
```

This will produce 3 executable files under `build/bin` directory. 

To run an simple example, use the following instructions: 
```
build/bin/two_party      example12.json
build/bin/two_party_iblt example12.json
build/bin/two_party_i64  example12.json
```

### Docker (macOS, or any machine without g++-14)

`docker/` contains a pinned GCC 14 environment. Only Docker is required:
```bash
docker/run.sh image                                   # build the image (once)
docker/run.sh build                                   # configure + compile (Release)
docker/run.sh two_party config/sc_unidirectional.json # run a config
docker/run.sh sims                                    # build the sc_study simulators
docker/run.sh shell                                   # interactive shell
```
The repo is mounted at `/work`, so configs and results are shared with the host;
relative `result filename`s are resolved from the repo root, and their directory must exist.
The build tree is kept in the Docker volume `commonsense-build` (binaries at `/build/bin`),
because the third-party builds need a case-sensitive filesystem.

## Reproducing Experiments

### Synthetic Data (Small Experiments)

We use a python script to automatically generate multiple configuration files from a configuration template (using libONIAK). The tool needs the absolute path of current directory to run. 

For example, if you want to run the configurations generated from `threshold.json`, you need to first use `pwd` to obtain the current absolute path, and then paste it to `root path` field. 

And also, **make sure that result folder exists**, e.g. in the following configuration template, we need to ensure the folder `result/one_side_quantize` exists under your previous configured root path. 

```json
{
  "script": [
    "ROOT/build/bin/two_party",
    "JSON"
  ],
  "root path": "/media/mydrive/CDMA", # <--- modify this
  "max process": 30,
  "seed": [
    "ITER_RANDOM",
    1000
  ],
  "universe": 2000000,
  "A size": 1000000,
  "B size": [
    "ITER_ALL",
    1010000,
    1020000,
    1050000,
    1100000,
    1200000
  ],
  "A minus B size": 0,
  "k": 8,
  "df": [
    "ITER_WITH",
    "B size", {
      "1010000": [7.75, 7.8, 7.82, 7.84, 7.86, 7.88, 7.9, 7.95, 8.0],
      "1020000": [6.42, 6.44, 6.46, 6.48, 6.5, 6.52, 6.54, 6.56, 6.58],
      "1050000": [5.02, 5.04, 5.06, 5.08, 5.1, 5.12, 5.14, 5.16, 5.18],
      "1100000": [4.12, 4.14, 4.16, 4.18, 4.2, 4.22, 4.24, 4.26, 4.28],
      "1200000": [3.32, 3.34, 3.36, 3.38, 3.4, 3.42, 3.44, 3.46, 3.48]
    }
  ],
  "d": [
    "EVAL",
    "int((${B size} - ${A size}) * ${df})"
  ],
  "max num peels": [
    "EVAL",
    "(${B size} - ${A size}) * 2.5"
  ],
  "ta": 100,
  "max rounds": 1000,
  "max comm rounds": 1,
  "counting": true,
  "resolving round": -1,
  "result filename": "ROOT/result/one_side_quantize/bs{B size}_d{d}_sd{seed}.json"
}
```

Then, run: 

```bash
python libONIAK/auto_run_experiment.py config/threshold.json
```

### Ethereum Data

#### Download Ethereum Data

We provide pre-processed ethereum data on google drive. 

```bash
https://drive.google.com/drive/folders/1615yrvJrhqMQTi6FGU0nP67yOWQlll96?usp=sharing
```

#### Run Ethereum Experiments

```bash
P='YOUR PATH TO DOWNLOADED DATA' bash script/run_geth_experiment.sh
```

#### Cleaning Ethereum Data

In case you want to compare other ethereum states, we provide a brief guide to download and cleanup these data. 

1. Download data from [public node](https://www.publicnode.com/snapshots). 
    For example, at the current of writing this guide, the newest snapshots can be downloaded using the following commands. 
    ```bash
    wget https://snapshots.publicnode.com/ethereum-geth-base-0-20353804.tar.lz4
    wget https://snapshots.publicnode.com/ethereum-geth-part-20353805-23093071.tar.lz4
    ```
2. Decompress these files, using `tar` with `lz4` algorithm. 
    ```bash
    tar -I lz4 -xf ethereum-geth-base-0-20353804.tar.lz4
    tar -I lz4 -xf ethereum-geth-part-20353805-23093071.tar.lz4
    ```
3. Run snapshot dump to export state infomation. By the time of writing, we use geth version `1.15.7-stable-827d3fcc`. 
    ```bash
    geth snapshot dump 2> /dev/null > state-23093071.txt
    ```
4. To prepare for diff, run external sorting. We provide a tool for that [(omni)](https://github.com/Y-jiji/omni-rs.git), but other methods may serve the same purpose. 
    ```bash
    omni external-sort --input state-23093071-sha256.txt --output state-23093071-sha256-sorted.txt
    ```
5. Compare it to another sorted state file to acquire the intersection and differences between two different states. We provide a tool for that [(omni)](https://github.com/Y-jiji/omni-rs.git), but other methods may serve the same purpose. 
    ```bash
    omni diff-sorted-string --input-a state-22020359-sha256-sorted.txt --input-b state-23093071-sha256-sorted.txt --output-a-minus-b state-22020359-minus-23093071.txt --output-b-minus-a state-23093071-minus-22020359.txt --output-intersect state-22020359-intersect-23093071.txt
    ```

## Acknowledgement

- Hash Function: https://github.com/wangyi-fudan/wyhash
- JSON Processor: https://github.com/nlohmann/json
- RANS Compression: https://github.com/loxxous/rans-tricks (we self-hosted our copy on https://github.com/Y-jiji/rans-fix)
- Updatable Priority Queue: https://github.com/Ten0/updatable_priority_queue
- BCH codec: https://github.com/mborgerding/bch_codec (we self-hosted our copy on https://github.com/Y-jiji/bch_codec)