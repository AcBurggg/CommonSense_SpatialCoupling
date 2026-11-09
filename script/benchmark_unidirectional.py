#!/usr/bin/env python3
"""Paired, reproducible comparison of the unidirectional CommonSense decoders."""

import argparse
import csv
import json
import math
import re
import statistics
import subprocess
from pathlib import Path


VARIANTS = {
    "original_mp": {
        "decoder": "mp",
        "ratios": (7.2, 7.4, 7.6, 7.8, 8.0, 8.2),
    },
    "uniform_bp": {
        "decoder": "bp",
        "ratios": (5.0, 5.2, 5.4, 5.6, 5.8, 6.0),
    },
    "spatial_bp": {
        "decoder": "bp",
        "ratios": (4.0, 4.2, 4.4, 4.6, 4.8, 5.0),
        "sc window fraction": 0.05,
        "sc tle": True,
    },
}

RUN_FIELDS = (
    "variant",
    "seed",
    "counter_ratio",
    "counters",
    "primary_decoder_success",
    "end_to_end_success",
    "communication_bits",
    "communication_over_lower_bound",
    "protocol_seconds",
    "communication_rounds",
    "peak_memory_kib",
)

SUMMARY_FIELDS = (
    "variant",
    "counter_ratio",
    "counters",
    "trials",
    "primary_decoder_success_rate",
    "end_to_end_success_rate",
    "mean_communication_bits",
    "mean_bits_on_primary_success",
    "mean_communication_over_lower_bound",
    "median_protocol_seconds",
    "median_communication_rounds",
    "median_peak_memory_kib",
    "minimum_ratio_at_95pct_primary_success",
)


def parse_peak_memory_kib(value):
    match = re.search(r"VmHWM:\s*(\d+)\s*kB", value or "")
    return int(match.group(1)) if match else ""


def lower_bound_bits(n, d):
    return (
        math.lgamma(n + 1)
        - math.lgamma(d + 1)
        - math.lgamma(n - d + 1)
    ) / math.log(2)


def write_csv(path, rows, fields):
    with path.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trials", type=int, default=5, help="paired dataset seeds per ratio (default: 5)")
    parser.add_argument("--seed-start", type=int, default=1, help="first dataset seed (default: 1)")
    parser.add_argument("--a-size", type=int, default=1_000_000, help="size of A (default: 1000000)")
    parser.add_argument("--delta-size", type=int, default=10_000, help="size of B\\A (default: 10000)")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("/work/result/unidirectional_benchmark"),
        help="directory for configs, raw JSON, logs, and CSVs",
    )
    args = parser.parse_args()
    if args.trials < 1 or args.a_size < 1 or args.delta_size < 1:
        parser.error("--trials, --a-size, and --delta-size must be positive")
    minimum_spatial_counters = round(
        args.delta_size * min(VARIANTS["spatial_bp"]["ratios"])
    )
    if round(0.05 * minimum_spatial_counters) < 4 * 8:
        parser.error(
            "delta size is too small for the 5% spatial window at the minimum "
            "sweep ratio; increase --delta-size or widen the spatial ratio sweep"
        )

    output_dir = args.output_dir.resolve()
    configs_dir = output_dir / "configs"
    results_dir = output_dir / "runs"
    logs_dir = output_dir / "logs"
    for directory in (configs_dir, results_dir, logs_dir):
        directory.mkdir(parents=True, exist_ok=True)

    b_size = args.a_size + args.delta_size
    bound = lower_bound_bits(b_size, args.delta_size)
    records = []
    binary = Path("/build/bin/two_party")

    for variant, settings in VARIANTS.items():
        for ratio in settings["ratios"]:
            counters = round(args.delta_size * ratio)
            for trial in range(args.trials):
                seed = args.seed_start + trial
                run_id = f"{variant}_r{ratio:.1f}_seed{seed}"
                result_path = results_dir / f"{run_id}.json"
                config_path = configs_dir / f"{run_id}.json"
                config = {
                    "seed": seed,
                    "A size": args.a_size,
                    "B size": b_size,
                    "A minus B size": 0,
                    "k": 8,
                    "d": counters,
                    "max num peels": round(args.delta_size * 2.5),
                    "counting": True,
                    "max comm rounds": 1,
                    "resolving round": -1,
                    "decoder": settings["decoder"],
                    "result filename": str(result_path),
                    **{
                        key: value
                        for key, value in settings.items()
                        if key != "ratios" and key != "decoder"
                    },
                }
                config_path.write_text(json.dumps(config, indent=2) + "\n")
                completed = subprocess.run(
                    [str(binary), str(config_path)],
                    cwd="/work",
                    text=True,
                    capture_output=True,
                    check=False,
                )
                (logs_dir / f"{run_id}.log").write_text(
                    completed.stdout + completed.stderr
                )
                if completed.returncode != 0:
                    raise RuntimeError(
                        f"{run_id} failed with exit code {completed.returncode}; "
                        f"see {logs_dir / f'{run_id}.log'}"
                    )
                if not result_path.exists():
                    raise RuntimeError(f"{run_id} did not write {result_path}")

                result = json.loads(result_path.read_text())
                cost = result["total communication cost"]
                times = result.get("time", [])
                records.append(
                    {
                        "variant": variant,
                        "seed": seed,
                        "counter_ratio": ratio,
                        "counters": counters,
                        "primary_decoder_success": result["primary decoder success"],
                        "end_to_end_success": result["success"],
                        "communication_bits": cost,
                        "communication_over_lower_bound": cost / bound,
                        "protocol_seconds": times[-1] if times else "",
                        "communication_rounds": result["communication rounds"],
                        "peak_memory_kib": parse_peak_memory_kib(
                            result.get("peak memory", "")
                        ),
                    }
                )
                print(
                    f"{run_id}: primary="
                    f"{records[-1]['primary_decoder_success']} "
                    f"bits={cost:.0f}",
                    flush=True,
                )

    summary = []
    for variant, settings in VARIANTS.items():
        for ratio in settings["ratios"]:
            group = [
                row
                for row in records
                if row["variant"] == variant and row["counter_ratio"] == ratio
            ]
            primary_rate = sum(row["primary_decoder_success"] for row in group) / len(group)
            end_to_end_rate = sum(row["end_to_end_success"] for row in group) / len(group)
            successful_costs = [
                row["communication_bits"]
                for row in group
                if row["primary_decoder_success"]
            ]
            timed = [row["protocol_seconds"] for row in group if row["protocol_seconds"] != ""]
            memory = [row["peak_memory_kib"] for row in group if row["peak_memory_kib"] != ""]
            summary.append(
                {
                    "variant": variant,
                    "counter_ratio": ratio,
                    "counters": group[0]["counters"],
                    "trials": len(group),
                    "primary_decoder_success_rate": primary_rate,
                    "end_to_end_success_rate": end_to_end_rate,
                    "mean_communication_bits": statistics.mean(
                        row["communication_bits"] for row in group
                    ),
                    "mean_bits_on_primary_success": (
                        statistics.mean(successful_costs) if successful_costs else ""
                    ),
                    "mean_communication_over_lower_bound": (
                        statistics.mean(successful_costs) / bound if successful_costs else ""
                    ),
                    "median_protocol_seconds": statistics.median(timed) if timed else "",
                    "median_communication_rounds": statistics.median(
                        row["communication_rounds"] for row in group
                    ),
                    "median_peak_memory_kib": statistics.median(memory) if memory else "",
                    "minimum_ratio_at_95pct_primary_success": "",
                }
            )

    for variant in VARIANTS:
        qualifying = [
            row["counter_ratio"]
            for row in summary
            if row["variant"] == variant
            and row["primary_decoder_success_rate"] >= 0.95
        ]
        threshold = min(qualifying) if qualifying else ""
        for row in summary:
            if row["variant"] == variant:
                row["minimum_ratio_at_95pct_primary_success"] = threshold

    write_csv(output_dir / "runs.csv", records, RUN_FIELDS)
    write_csv(output_dir / "summary.csv", summary, SUMMARY_FIELDS)
    print(f"\nWrote {output_dir / 'runs.csv'} and {output_dir / 'summary.csv'}")
    print(f"SetX counting lower bound for this workload: {bound:.1f} bits")
    threshold_rows = {}
    for variant in VARIANTS:
        threshold = next(
            row["minimum_ratio_at_95pct_primary_success"]
            for row in summary
            if row["variant"] == variant
        )
        print(
            f"{variant}: minimum tested d/|B\\A| at >=95% primary success: "
            f"{threshold if threshold != '' else 'not reached'}"
        )
        if threshold != "":
            threshold_rows[variant] = next(
                row
                for row in summary
                if row["variant"] == variant and row["counter_ratio"] == threshold
            )

    if len(threshold_rows) == len(VARIANTS):
        print("\nMetrics at each variant's minimum tested ratio meeting the target:")
        print("variant       d/|B\\A|  counters  primary success  mean bits on primary success")
        for variant, row in threshold_rows.items():
            print(
                f"{variant:<13} {row['counter_ratio']:>7.1f} "
                f"{row['counters']:>9} "
                f"{row['primary_decoder_success_rate']:>15.1%} "
                f"{row['mean_bits_on_primary_success']:>28.0f}"
            )
        for control in ("original_mp", "uniform_bp"):
            reference = threshold_rows[control]
            coupled = threshold_rows["spatial_bp"]
            counter_reduction = 1 - coupled["counters"] / reference["counters"]
            bit_reduction = (
                1
                - coupled["mean_bits_on_primary_success"]
                / reference["mean_bits_on_primary_success"]
            )
            counter_change = (
                f"{counter_reduction:.1%} fewer"
                if counter_reduction >= 0
                else f"{-counter_reduction:.1%} more"
            )
            bit_change = (
                f"{bit_reduction:.1%} fewer"
                if bit_reduction >= 0
                else f"{-bit_reduction:.1%} more"
            )
            print(
                f"spatial_bp vs {control}: "
                f"{counter_change} counters, {bit_change} communication bits"
            )


if __name__ == "__main__":
    main()
