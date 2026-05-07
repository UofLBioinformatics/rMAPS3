# FAQ

Original web FAQ: <http://rmaps.cecsresearch.org/Help/FAQ>

## What changed between rMAPS3 and rMAPS2?

rMAPS3 keeps the five-event support from rMAPS2 (SE, MXE, A5SS, A3SS, and RI)
and modernizes the project for local use from this Git repository.

Major rMAPS3 updates include:

- Python 3 compatibility for the maintained CLI and local web UI workflows.
- A unified command-line entry point in `cli.py`.
- Local web UI support through `run_web.py`.
- Genome download/setup scripts and a documented genome-data layout.
- Selectable p-value methods for Motif Map and CLIP Map analyses.
- Version reporting with `python cli.py --version`.
- Event-specific CLIP RNA-map renderers and fallback PDF/PNG export when
  optional plotting dependencies are unavailable.

## Does rMAPS support frog genomes?

The original rMAPS2 web server documented support for two frog assemblies:

- `xenLae2`: Xenopus laevis, African clawed frog
- `xenTro9`: Xenopus tropicalis

In this repository, genome availability depends on the FASTA files installed in
the configured genome data directory. See [INSTALL.md](INSTALL.md) for genome
setup details.

## Do I need to use the provided genome fetch scripts?

No. The fetch scripts are convenience helpers. You may also download or prepare
FASTA files yourself, as long as they follow the expected directory and filename
layout.

## What input should I use: rMATS, MISO, or coordinate files?

Use rMATS output when possible. It is the simplest and most directly supported
workflow. MISO output is supported through conversion scripts. Coordinate files
are useful when you already have curated upregulated, downregulated, and
background event sets.

## Troubleshooting

### I get `FastaNotFoundError`. What should I check?

Verify that `--genome` and `--fasta-root` point to a FASTA file at:

```text
<fasta-root>/<build>/<build>.fa
```

For example, `--genome hg38 --fasta-root genomedata` expects:

```text
genomedata/hg38/hg38.fa
```

### MISO conversion fails. What should I check?

Make sure Perl is installed and available on `PATH`, and verify the converter
scripts exist under `bin/`:

```text
bin/miso2rMATS.SE.pl
bin/miso2rMATS.A3SS.pl
bin/miso2rMATS.A5SS.pl
bin/miso2rMATS.RI.pl
bin/miso2rMATS.MXE.pl
```

### CLIP plots are missing. What should I check?

Confirm the input produces non-empty upregulated, downregulated, and background
event groups. If those groups are empty, rMAPS can write partial logs or
intermediate files but cannot produce meaningful CLIP RNA maps.

Also check optional plotting dependencies such as TeX and Ghostscript. rMAPS3
includes fallback export behavior, but full plotting support is best with those
tools installed.
