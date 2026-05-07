#!/usr/bin/env python3
"""
CLIP-seq RNA map generation for Retained Intron (RI) events.

This program processes rMATS RI output and CLIP-seq peaks to generate RNA maps.
"""

import sys
import logging
import time
import argparse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from rmaps_core.clip_utils import (
    setup_logging, copy_file, process_rmats_file, convert_miso_to_rmats
)


def main():
    parser = argparse.ArgumentParser(description='CLIP-seq RNA map for RI events')
    parser.add_argument('-p', '--peak', dest='peakFile', required=True, help='CLIP-seq peaks')
    parser.add_argument('-o', '--output', dest='output', required=True, help='output dir')
    parser.add_argument('-r', '--rMATS', dest='rMATS', required=True, help='rMATS RI')
    parser.add_argument('-mi', '--miso', dest='miso', required=True, help='MISO RI')
    parser.add_argument('-u', '--up', dest='up', required=True, help='upregulated')
    parser.add_argument('-d', '--down', dest='dn', required=True, help='downregulated')
    parser.add_argument('-b', '--background', dest='bg', required=True, help='background')
    parser.add_argument('--label', default="RBP", help='RBP label')
    parser.add_argument('--intron', type=int, default=250, help='intron length')
    parser.add_argument('--exon', type=int, default=50, help='exon length')
    parser.add_argument('--window', type=int, default=10, help='window size')
    parser.add_argument('--step', type=int, default=1, help='step size')
    parser.add_argument('--sigFDR', type=float, default=0.05, help='FDR cutoff')
    parser.add_argument('--sigDeltaPSI', type=float, default=0.05, help='Delta PSI cutoff')
    parser.add_argument('--separate', default=False, action='store_true', help='Separate plots')
    
    args = parser.parse_args()
    
    out_path = Path(args.output).resolve()
    out_path.mkdir(parents=True, exist_ok=True)
    exon_path = out_path / 'exon'
    exon_path.mkdir(parents=True, exist_ok=True)
    temp_path = out_path / 'temp'
    temp_path.mkdir(parents=True, exist_ok=True)
    
    script_path = Path(__file__).resolve().parent.parent
    bin_path = script_path / 'bin'
    
    setup_logging(out_path, "3.0.0")
    logging.debug('Start: %s', ' '.join(sys.argv))
    start_time = time.time()
    
    rmats = args.rMATS
    miso = args.miso
    
    if rmats == "NA" and miso == "NA" and (args.up == "NA" or args.dn == "NA" or args.bg == "NA"):
        print("Error: Need rMATS/MISO or coordinate files")
        sys.exit(-99)
    
    if '_' in args.label:
        print('Error: Label cannot contain underscore')
        sys.exit(-1)
    
    if miso != "NA":
        try:
            rmats = convert_miso_to_rmats(miso, temp_path, bin_path, 'ri')
        except RuntimeError as e:
            logging.error(str(e))
            sys.exit(-101)
    
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
            # RI coordinates: chr, strand, riExonStart, riExonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE
            coord_indices = (3, 4, 5, 6, 7, 8, 9, 10)
            header = 'chr\tstrand\triExonStart_0base\triExonEnd\tupstreamES\tupstreamEE\tdownstreamES\tdownstreamEE'
            nu, nd, nb = process_rmats_file(rmats, exon_path, header, coord_indices, args.sigFDR, args.sigDeltaPSI)
        
        logging.debug("Counts - Up: %d, Down: %d, Bg: %d", nu, nd, nb)
    except Exception as e:
        logging.error("Exception: %s", str(e))
        sys.exit(-1)
    
    logging.debug("GENERATING RNA MAP")
    
    import subprocess
    cmd = [
        sys.executable, str(bin_path / 'RNA.map.noWiggle.RI.py'),
        str(exon_path), str(Path(args.peakFile).resolve()),
        str(args.intron), str(args.exon), str(args.window), str(args.step), str(args.sigFDR),
        args.label, str(out_path), str(nu), str(nd), str(nb), str(int(args.separate)), 'RI'
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logging.error("RNA map failed: %s", result.stderr)
        sys.exit(-9)
    
    runtime = time.time() - start_time
    logging.debug("Runtime: %02d:%02d:%02d", int(runtime/3600), int((runtime%3600)/60), int(runtime%60))
    print(f"CLIP-seq RI analysis completed! Results: {out_path}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
