#!/usr/bin/env python3
"""
Unified Motif Map Test Suite

Usage:
    python tests/test_motif.py            # Test all event types
    python tests/test_motif.py --event se # Test specific event type
    python tests/test_motif.py --list     # List available event types
"""

import argparse
import subprocess
import sys
from pathlib import Path
from datetime import datetime


# Motif event type configurations.
# Uses bundled event-specific rMATS files under data/test/clip/.
MOTIF_EVENTS = {
    "se": {
        "name": "Skipped Exon (SE)",
            "rMATS": "data/test/clip/SE/se.rMATS.txt",
        "output": "results/test_motif_se",
    },
    "a5ss": {
        "name": "Alternative 5' Splice Site (A5SS)",
        "rMATS": "data/test/clip/A5SS/a5ss.rMATS.txt",
        "output": "results/test_motif_a5ss",
    },
    "a3ss": {
        "name": "Alternative 3' Splice Site (A3SS)",
        "rMATS": "data/test/clip/A3SS/a3ss.rMATS.txt",
        "output": "results/test_motif_a3ss",
    },
    "mxe": {
        "name": "Mutually Exclusive Exons (MXE)",
        "rMATS": "data/test/clip/MXE/mxe.rMATS.txt",
        "output": "results/test_motif_mxe",
    },
    "ri": {
        "name": "Retained Intron (RI)",
        "rMATS": "data/test/clip/RI/ri.rMATS.txt",
        "output": "results/test_motif_ri",
    },
}

COMMON_PARAMS = {
    "genome": "hg19",
    "known_motifs": "data/testMotifs.txt",
    "motifs": "NA",
    "miso": "NA",
    "up": "NA",
    "down": "NA",
    "background": "NA",
}

DEFAULT_COMPARE_METHODS = [
    "fisher",
    "mannwhitney_greater",
    "brunnermunzel_greater",
]


def run_motif_test(event_type, config, fasta_root, verbose=False):
    """Run motif-map test for a specific event type."""
    print(f"\n{'='*70}")
    print(f"Testing Motif Map: {config['name']}")
    print(f"{'='*70}")

    cmd = [
        sys.executable, "cli.py", "motif-map", event_type,
        "--known-motifs", COMMON_PARAMS["known_motifs"],
        "--motifs", COMMON_PARAMS["motifs"],
        "--fasta-root", fasta_root,
        "--genome", COMMON_PARAMS["genome"],
        "--output", config["output"],
        "--rMATS", config["rMATS"],
        "--miso", COMMON_PARAMS["miso"],
        "--up", COMMON_PARAMS["up"],
        "--down", COMMON_PARAMS["down"],
        "--background", COMMON_PARAMS["background"],
    ]

    if verbose:
        print(f"Command: {' '.join(cmd)}")

    start_time = datetime.now()
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=Path(__file__).parent.parent,
        )

        elapsed = (datetime.now() - start_time).total_seconds()

        if result.returncode == 0:
            print(f"PASS in {elapsed:.1f}s")
            if verbose and result.stdout:
                print(f"Output: {result.stdout}")
            return True

        print(f"FAIL (exit code {result.returncode})")
        if result.stderr:
            print(f"Error: {result.stderr}")
        if result.stdout:
            print(f"Output: {result.stdout}")
        return False
    except Exception as e:
        print(f"FAIL with exception: {e}")
        return False


def count_significant_windows(pval_file: Path, threshold: float = 0.05) -> int:
    if not pval_file.exists():
        return -1
    count = 0
    with pval_file.open("r", encoding="utf-8", errors="replace") as handle:
        next(handle, None)
        for line in handle:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            try:
                if float(parts[2]) < threshold:
                    count += 1
            except ValueError:
                continue
    return count


def run_motif_method_compare(event_type, config, fasta_root, methods, verbose=False):
    print(f"\n{'='*70}")
    print(f"Comparing Motif methods: {config['name']}")
    print(f"{'='*70}")

    base_out = Path(config["output"])
    method_outputs = {method: f"{base_out}_{method}" for method in methods}

    for method, out_dir in method_outputs.items():
        cmd = [
            sys.executable, "cli.py", "motif-map", event_type,
            "--known-motifs", COMMON_PARAMS["known_motifs"],
            "--motifs", COMMON_PARAMS["motifs"],
            "--fasta-root", fasta_root,
            "--genome", COMMON_PARAMS["genome"],
            "--output", out_dir,
            "--rMATS", config["rMATS"],
            "--miso", COMMON_PARAMS["miso"],
            "--up", COMMON_PARAMS["up"],
            "--down", COMMON_PARAMS["down"],
            "--background", COMMON_PARAMS["background"],
            "--stat-method", method,
        ]
        if method == "permutation_one_sided":
            cmd.extend(["--stat-permutations", "200", "--stat-seed", "1337"])
        if verbose:
            print(f"Command: {' '.join(cmd)}")
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=Path(__file__).parent.parent,
        )
        if result.returncode != 0:
            print(f"FAIL ({method}) exit code {result.returncode}")
            if result.stderr:
                print(result.stderr)
            return False

    summary = []
    for method, out_dir in method_outputs.items():
        up_file = Path(out_dir) / "pVal.up.vs.bg.RNAmap.txt"
        dn_file = Path(out_dir) / "pVal.dn.vs.bg.RNAmap.txt"
        up_sig = count_significant_windows(up_file)
        dn_sig = count_significant_windows(dn_file)
        if min(up_sig, dn_sig) < 0:
            print(f"FAIL (missing motif p-value files for method={method})")
            return False
        summary.append(
            (method, up_sig, dn_sig)
        )

    for method, up_sig, dn_sig in summary:
        print(f"Significant windows (p<0.05): method={method}, up={up_sig}, dn={dn_sig}")

    print("PASS method comparison")
    return True


def list_events():
    """List all available event types."""
    print("\nAvailable Motif event types:")
    print("-" * 50)
    for event_type, config in MOTIF_EVENTS.items():
        print(f"  {event_type:6s} - {config['name']}")
    print()


def main():
    parser = argparse.ArgumentParser(
        description="Unified Motif Map Test Suite",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python tests/test_motif.py                    # Test all event types
  python tests/test_motif.py --event se         # Test SE only
  python tests/test_motif.py --event se a5ss    # Test SE and A5SS
  python tests/test_motif.py --verbose          # Show detailed output
  python tests/test_motif.py --list             # List event types
        """,
    )

    parser.add_argument(
        "--event",
        nargs="+",
        choices=list(MOTIF_EVENTS.keys()),
        help="Specific event type(s) to test (default: all)",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List available event types and exit",
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Show detailed output",
    )
    parser.add_argument(
        "--fasta-root",
        default="genomedata",
        help="Path to FASTA root directory (default: genomedata)",
    )
    parser.add_argument(
        "--genome",
        default="hg19",
        help="Genome build to use (default: hg19)",
    )
    parser.add_argument(
        "--compare-methods",
        action="store_true",
        help="Run selected method set and verify all method paths produce valid outputs.",
    )
    parser.add_argument(
        "--include-permutation",
        action="store_true",
        help="Include permutation_one_sided in --compare-methods runs (slower).",
    )
    args = parser.parse_args()

    if args.list:
        list_events()
        return 0

    COMMON_PARAMS["genome"] = args.genome
    events_to_test = args.event if args.event else list(MOTIF_EVENTS.keys())
    compare_methods = list(DEFAULT_COMPARE_METHODS)
    if args.include_permutation:
        compare_methods.append("permutation_one_sided")

    print("\n" + "=" * 70)
    print("MOTIF MAP TEST SUITE")
    print("=" * 70)
    print(f"Testing {len(events_to_test)} event type(s): {', '.join(events_to_test)}")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"FASTA root: {args.fasta_root}")
    print(f"Genome: {args.genome}")

    results = {}
    for event_type in events_to_test:
        config = MOTIF_EVENTS[event_type]
        if args.compare_methods:
            success = run_motif_method_compare(event_type, config, args.fasta_root, compare_methods, args.verbose)
        else:
            success = run_motif_test(event_type, config, args.fasta_root, args.verbose)
        results[event_type] = success

    print("\n" + "=" * 70)
    print("TEST SUMMARY")
    print("=" * 70)

    passed = sum(1 for success in results.values() if success)
    failed = len(results) - passed

    for event_type, success in results.items():
        status = "PASS" if success else "FAIL"
        print(f"  {event_type:6s} - {status}")

    print("-" * 70)
    print(f"Total: {len(results)} tests, {passed} passed, {failed} failed")
    print("=" * 70)

    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())


