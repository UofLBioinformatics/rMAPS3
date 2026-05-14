from pathlib import Path
import sys
from unittest.mock import patch

import typer

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import cli
from rmaps_core.clip_core import run_clip_map
from rmaps_core.motif_map_core import run_motif_map
from rmaps_core.path_utils import resolve_user_path


def expected_path(base_cwd: Path, relative_path: str) -> str:
    return str((base_cwd / relative_path).resolve())


def test_resolve_user_path_preserves_na() -> None:
    assert resolve_user_path("NA", Path.cwd()) == "NA"


def test_resolve_user_path_uses_invocation_cwd() -> None:
    base_cwd = Path.cwd() / "scratch_user_run"
    assert resolve_user_path("results/out", base_cwd) == expected_path(base_cwd, "results/out")


def test_motif_map_resolves_user_paths_and_runs_from_invocation_cwd() -> None:
    base_cwd = Path.cwd() / "scratch_user_run"
    with patch("rmaps_core.motif_map_core.maybe_prepare_rmats_input", side_effect=lambda rmats, output: rmats), patch(
        "rmaps_core.motif_map_core.subprocess.run"
    ) as run_mock:
        run_mock.return_value.returncode = 0

        run_motif_map(
            "se",
            known_motifs=Path("data/known.txt"),
            motifs="motifs.txt",
            fasta_root=Path("genomedata"),
            genome="mm39",
            output=Path("results/motif_se"),
            rmats="events.tsv",
            miso="NA",
            up="NA",
            down="NA",
            background="NA",
            label="RBP",
            intron=250,
            exon=50,
            window=50,
            step=1,
            sig_fdr=0.05,
            sig_delta_psi=0.05,
            separate=False,
            base_cwd=base_cwd,
        )

    cmd = run_mock.call_args.args[0]
    assert run_mock.call_args.kwargs["cwd"] == base_cwd
    assert cmd[cmd.index("-k") + 1] == expected_path(base_cwd, "data/known.txt")
    assert cmd[cmd.index("-m") + 1] == expected_path(base_cwd, "motifs.txt")
    assert cmd[cmd.index("--fasta-root") + 1] == expected_path(base_cwd, "genomedata")
    assert cmd[cmd.index("-o") + 1] == expected_path(base_cwd, "results/motif_se")
    assert cmd[cmd.index("-r") + 1] == expected_path(base_cwd, "events.tsv")
    assert cmd[cmd.index("-mi") + 1] == "NA"


def test_clip_map_resolves_user_paths_and_runs_from_invocation_cwd() -> None:
    base_cwd = Path.cwd() / "scratch_user_run"
    with patch("rmaps_core.clip_core.maybe_prepare_rmats_input", side_effect=lambda rmats, output: rmats), patch(
        "rmaps_core.clip_core.subprocess.run"
    ) as run_mock:
        run_mock.return_value.returncode = 0

        run_clip_map(
            "se",
            peak=Path("peaks.bed"),
            output=Path("results/clip_se"),
            rmats="events.tsv",
            miso="NA",
            up="NA",
            down="NA",
            background="NA",
            label="RBP",
            intron=250,
            exon=50,
            window=10,
            step=1,
            sig_fdr=0.05,
            sig_delta_psi=0.05,
            separate=False,
            base_cwd=base_cwd,
        )

    cmd = run_mock.call_args.args[0]
    assert run_mock.call_args.kwargs["cwd"] == base_cwd
    assert cmd[cmd.index("--peak") + 1] == expected_path(base_cwd, "peaks.bed")
    assert cmd[cmd.index("--output") + 1] == expected_path(base_cwd, "results/clip_se")
    assert cmd[cmd.index("--rMATS") + 1] == expected_path(base_cwd, "events.tsv")
    assert cmd[cmd.index("--miso") + 1] == "NA"


def test_convert_miso_resolves_input_and_output_from_invocation_cwd() -> None:
    base_cwd = Path.cwd() / "scratch_user_run"
    with patch("cli.Path.cwd", return_value=base_cwd), patch("cli.run_subprocess") as run_mock:
        run_mock.return_value = 0

        try:
            cli.convert_miso(
                event="se",
                miso=Path("inputs/events.miso"),
                output=Path("results/events.rmats.txt"),
                bayes_low=1,
                bayes_high=100,
            )
        except typer.Exit as exc:
            exit_code = exc.exit_code
        else:
            raise AssertionError("convert_miso did not exit")

    assert exit_code == 0
    cmd = run_mock.call_args.args[0]
    assert run_mock.call_args.kwargs["base_cwd"] == base_cwd
    assert cmd[-2:] == [
        expected_path(base_cwd, "inputs/events.miso"),
        expected_path(base_cwd, "results/events.rmats.txt"),
    ]


def test_exon_sets_resolves_input_and_output_from_invocation_cwd() -> None:
    base_cwd = Path.cwd() / "scratch_user_run"
    with patch("cli.Path.cwd", return_value=base_cwd), patch("cli.run_subprocess") as run_mock:
        run_mock.return_value = 0

        try:
            cli.exon_sets_se(
                input_file=Path("inputs/common.txt"),
                sample1="sample_a",
                sample2="sample_b",
                out_dir=Path("results/exons"),
            )
        except typer.Exit as exc:
            exit_code = exc.exit_code
        else:
            raise AssertionError("exon_sets_se did not exit")

    assert exit_code == 0
    cmd = run_mock.call_args.args[0]
    assert run_mock.call_args.kwargs["base_cwd"] == base_cwd
    assert cmd[-4:] == [
        expected_path(base_cwd, "inputs/common.txt"),
        "sample_a",
        "sample_b",
        expected_path(base_cwd, "results/exons"),
    ]


def main() -> None:
    test_resolve_user_path_preserves_na()
    test_resolve_user_path_uses_invocation_cwd()
    test_motif_map_resolves_user_paths_and_runs_from_invocation_cwd()
    test_clip_map_resolves_user_paths_and_runs_from_invocation_cwd()
    test_convert_miso_resolves_input_and_output_from_invocation_cwd()
    test_exon_sets_resolves_input_and_output_from_invocation_cwd()
    print("Path handling tests completed.")


if __name__ == "__main__":
    main()
