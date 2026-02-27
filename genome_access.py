from pathlib import Path
from functools import lru_cache

from pyfaidx import Fasta


_RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(_RC_TABLE)[::-1]


@lru_cache(maxsize=None)
def load_genome(build: str, base_dir: str) -> Fasta:
    """
    Load a genome FASTA for the given build from a genomedata-style layout:
    base_dir/<build>/<build>.fa
    """
    base = Path(base_dir)
    fasta_path = base / build / f"{build}.fa"
    return Fasta(str(fasta_path), as_raw=True, sequence_always_upper=True)


def fetch_seq(fasta: Fasta, strand: str, chrom: str, start: int, end: int) -> str:
    """
    Fetch sequence from a pyfaidx Fasta, matching legacy sequence behavior:
    - 0-based, end-exclusive slicing
    - reverse-complement on '-' strand
    - on any error, return a string of Ns of the requested length
    """
    length = end - start
    if length <= 0:
        return ""
    try:
        record = fasta[chrom][start:end]
        seq = str(record)
        if strand == "-":
            seq = revcomp(seq)
        if len(seq) != length:
            if len(seq) < length:
                seq = (seq + ("N" * length))[:length]
            else:
                seq = seq[:length]
        return seq.upper()
    except Exception:
        return "N" * length

