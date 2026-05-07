#!/usr/bin/env python3
"""
CLIP-seq RNA map generation for Alternative 3' Splice Site (A3SS) events.

This program processes rMATS A3SS output and CLIP-seq peaks to generate RNA maps.
"""

import sys
import os
import logging
import time
import argparse
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from rmaps_core.clip_utils import (
    setup_logging,
    copy_file,
    process_rmats_file,
    convert_miso_to_rmats
)


def main():
    """Main execution function for A3SS CLIP-seq analysis."""
    
    parser = argparse.ArgumentParser(
        description='Making RNA map from CLIP-seq peaks and A3SS exon coordinates'
    )
    parser.add_argument('-p', '--peak', dest='peakFile', required=True, help='CLIP-seq peak file')
    parser.add_argument('-o', '--output', dest='output', required=True, help='output directory')
    parser.add_argument('-r', '--rMATS', dest='rMATS', required=True, help='rMATS A3SS file')
    parser.add_argument('-mi', '--miso', dest='miso', required=True, help='MISO A3SS file')
    parser.add_argument('-u', '--up', dest='up', required=True, help='upregulated exons file')
    parser.add_argument('-d', '--down', dest='dn', required=True, help='downregulated exons file')
    parser.add_argument('-b', '--background', dest='bg', required=True, help='background exons file')
    parser.add_argument('--label', type=str, dest='label', default="RBP", help='RBP label')
    parser.add_argument('--intron', type=int, dest='intron', default=250, help='intron length')
    parser.add_argument('--exon', type=int, dest='exon', default=50, help='exon length')
    parser.add_argument('--window', type=int, dest='window', default=10, help='window size')
    parser.add_argument('--step', type=int, dest='step', default=1, help='step size')
    parser.add_argument('--sigFDR', type=float, dest='sigFDR', default=0.05, help='FDR cutoff')
    parser.add_argument('--sigDeltaPSI', type=float, dest='sigDeltaPSI', default=0.05, help='Delta PSI cutoff')
    parser.add_argument('--separate', dest='separate', default=False, action='store_true', help='Separate plots')
    
    args = parser.parse_args()
    
    out_path = Path(args.output).resolve()
    out_path.mkdir(parents=True, exist_ok=True)
    exon_path = out_path / 'exon'
    exon_path.mkdir(parents=True, exist_ok=True)
    temp_path = out_path / 'temp'
    temp_path.mkdir(parents=True, exist_ok=True)
    
    script_path = Path(__file__).resolve().parent.parent
    bin_path = script_path / 'bin'
    
    version = "3.0.0"
    setup_logging(out_path, version)
    logging.debug('Start: %s', ' '.join(sys.argv))
    start_time = time.time()
    
    rmats = args.rMATS
    miso = args.miso
    
    if rmats == "NA" and miso == "NA" and (args.up == "NA" or args.dn == "NA" or args.bg == "NA"):
        print("Error: Need rMATS/MISO or all coordinate files")
        sys.exit(-99)
    
    if '_' in args.label:
        print('Error: Label cannot contain underscore')
        sys.exit(-1)
    
    if not 0.0 <= args.sigFDR <= 1.0:
        print('Error: FDR must be 0.0-1.0')
        sys.exit(-1)
    
    if miso != "NA":
        try:
            rmats = convert_miso_to_rmats(miso, temp_path, bin_path, 'a3ss')
        except RuntimeError as e:
            logging.error(str(e))
            sys.exit(-101)
    
    logging.debug("=" * 50)
    logging.debug("PROCESSING INPUT FILES")
    
    try:
        if rmats == "NA":
            copy_file(args.up, str(exon_path / 'up.coord.txt'))
            copy_file(args.dn, str(exon_path / 'dn.coord.txt'))
            copy_file(args.bg, str(exon_path / 'bg.coord.txt'))
            nu = sum(1 for _ in open(exon_path / 'up.coord.txt')) - 1
            nd = sum(1 for _ in open(exon_path / 'dn.coord.txt')) - 1
            nb = sum(1 for _ in open(exon_path / 'bg.coord.txt')) - 1
        else:
            coord_indices = (3, 4, 5, 6, 7, 8, 9, 10)
            header = 'chr\tstrand\tlongExonStart_0base\tlongExonEnd\tshortES\tshortEE\tflankingES\tflankingEE'
            nu, nd, nb = process_rmats_file(rmats, exon_path, header, coord_indices, args.sigFDR, args.sigDeltaPSI)
        
        logging.debug("Counts - Up: %d, Down: %d, Bg: %d", nu, nd, nb)
    except Exception as e:
        logging.error("Exception: %s", str(e))
        sys.exit(-1)
    
    logging.debug("=" * 50)
    logging.debug("GENERATING RNA MAP")
    
    import subprocess
    cmd = [
        sys.executable, str(bin_path / 'RNA.map.noWiggle.A3SS.py'),
        str(exon_path), str(Path(args.peakFile).resolve()),
        str(args.intron), str(args.exon), str(args.window), str(args.step), str(args.sigFDR),
        args.label, str(out_path), str(nu), str(nd), str(nb), str(int(args.separate)), 'A3SS'
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logging.error("RNA map failed: %s", result.stderr)
        sys.exit(-9)
    
    logging.debug("RNA map output:\n%s", result.stdout)
    
    runtime = time.time() - start_time
    logging.debug("Runtime: %02d:%02d:%02d", int(runtime/3600), int((runtime%3600)/60), int(runtime%60))
    print(f"CLIP-seq A3SS analysis completed! Results: {out_path}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
