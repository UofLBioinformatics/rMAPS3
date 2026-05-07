#!/usr/bin/env python3
"""
CLIP-seq RNA map generation for Alternative 5' Splice Site (A5SS) events.

This program processes rMATS A5SS output and CLIP-seq peaks to generate RNA maps.
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
    run_perl_command,
    copy_file,
    process_rmats_file,
    convert_miso_to_rmats
)


def main():
    """Main execution function for A5SS CLIP-seq analysis."""
    
    # Argument parsing
    parser = argparse.ArgumentParser(
        description='Making RNA map from CLIP-seq peaks and A5SS exon coordinates'
    )
    parser.add_argument('-p', '--peak', dest='peakFile', required=True,
                        help='CLIP-seq peak file')
    parser.add_argument('-o', '--output', dest='output', required=True,
                        help='output directory')
    parser.add_argument('-r', '--rMATS', dest='rMATS', required=True,
                        help='an rMATS output file from A5SS event')
    parser.add_argument('-mi', '--miso', dest='miso', required=True,
                        help='an miso output file from A5SS event')
    parser.add_argument(
        '-u', '--up', dest='up', required=True,
        help='tab delimited file containing upregulated exons. '
             'Headers: GeneID geneSymbol chr strand longExonStart_0base longExonEnd '
             'shortES shortEE flankingES flankingEE'
    )
    parser.add_argument(
        '-d', '--down', dest='dn', required=True,
        help='tab delimited file containing downregulated exons'
    )
    parser.add_argument(
        '-b', '--background', dest='bg', required=True,
        help='tab delimited file containing background exons'
    )
    parser.add_argument('--label', type=str, dest='label', default="RBP",
                        help='Label of RBP, e.g., PTB, ESRP')
    parser.add_argument('--intron', type=int, dest='intron', default=250,
                        help='intron length to examine (NTs from exon junction)')
    parser.add_argument('--exon', type=int, dest='exon', default=50,
                        help='exon length to examine (NTs from each end)')
    parser.add_argument('--window', type=int, dest='window', default=10,
                        help='window size for analysis')
    parser.add_argument('--step', type=int, dest='step', default=1,
                        help='step size for sliding window')
    parser.add_argument('--sigFDR', type=float, dest='sigFDR', default=0.005,
                        help='FDR cutoff for significant events')
    parser.add_argument('--sigDeltaPSI', type=float, dest='sigDeltaPSI', default=0.01,
                        help='Delta PSI cutoff for significant events')
    parser.add_argument('--separate', dest='separate', default=False,
                        action='store_true',
                        help='Create separate p-value plots')
    
    args = parser.parse_args()
    
    # Setup output directories
    out_path = Path(args.output).resolve()
    out_path.mkdir(parents=True, exist_ok=True)
    
    exon_path = out_path / 'exon'
    exon_path.mkdir(parents=True, exist_ok=True)
    
    temp_path = out_path / 'temp'
    temp_path.mkdir(parents=True, exist_ok=True)
    
    # Setup paths
    script_path = Path(__file__).resolve().parent.parent
    bin_path = script_path / 'bin'
    
    # Setup logging
    setup_logging(out_path)
    logging.debug('Start the program with arguments: %s', ' '.join(sys.argv))
    start_time = time.time()
    
    # Validate inputs
    rmats = args.rMATS
    miso = args.miso
    up = args.up
    dn = args.dn
    bg = args.bg
    
    if rmats == "NA" and miso == "NA" and (up == "NA" or dn == "NA" or bg == "NA"):
        print("Error: Need either rMATS input, MISO input, or coordinates for all three exon types")
        logging.error("Incorrect input: missing required files")
        sys.exit(-99)
    
    # Validate protein name
    if '_' in args.label:
        print('Error: Protein name should not contain underscore (_)')
        logging.error("Invalid protein name: contains underscore")
        sys.exit(-1)
    
    # Validate FDR
    if not 0.0 <= args.sigFDR <= 1.0:
        print('Error: FDR cutoff must be between 0.0 and 1.0')
        logging.error("Invalid FDR cutoff: %f", args.sigFDR)
        sys.exit(-1)
    
    peak_file_path = Path(args.peakFile).resolve()
    
    # Convert MISO to rMATS if needed
    if miso != "NA":
        logging.debug("Converting MISO format to rMATS")
        try:
            rmats = convert_miso_to_rmats(miso, temp_path, bin_path, 'a5ss')
        except RuntimeError as e:
            logging.error(str(e))
            sys.exit(-101)
    
    # Process rMATS or copy coordinate files
    logging.debug("=" * 50)
    logging.debug("MAKING INPUT FILES FROM rMATS FILE")
    
    try:
        if rmats == "NA":
            # Copy provided coordinate files
            logging.debug("Copying provided coordinate files")
            
            status, output = copy_file(up, str(exon_path / 'up.coord.txt'))
            logging.debug("Copying up file: %s", output)
            if status != 0:
                raise RuntimeError("Failed to copy up file")
            
            status, output = copy_file(dn, str(exon_path / 'dn.coord.txt'))
            logging.debug("Copying dn file: %s", output)
            if status != 0:
                raise RuntimeError("Failed to copy dn file")
            
            status, output = copy_file(bg, str(exon_path / 'bg.coord.txt'))
            logging.debug("Copying bg file: %s", output)
            if status != 0:
                raise RuntimeError("Failed to copy bg file")
            
            # Count lines in files (excluding header)
            nu = sum(1 for line in open(exon_path / 'up.coord.txt')) - 1
            nd = sum(1 for line in open(exon_path / 'dn.coord.txt')) - 1
            nb = sum(1 for line in open(exon_path / 'bg.coord.txt')) - 1
            
        else:
            # Process rMATS file
            # A5SS coordinate indices: chr, strand, longExonStart, longExonEnd, shortES, shortEE, flankingES, flankingEE
            # Typically columns 3-10 in rMATS output
            coord_indices = (3, 4, 5, 6, 7, 8, 9, 10)
            header = 'chr\tstrand\tlongExonStart_0base\tlongExonEnd\tshortES\tshortEE\tflankingES\tflankingEE'
            
            nu, nd, nb = process_rmats_file(
                rmats,
                exon_path,
                header,
                coord_indices,
                args.sigFDR,
                args.sigDeltaPSI
            )
        
        logging.debug("Final counts - Up: %d, Down: %d, Background: %d", nu, nd, nb)
        
    except Exception as e:
        logging.error("Exception in making input from rMATS output: %s", str(e))
        logging.error("Exception type: %s", type(e).__name__)
        sys.exit(-1)
    
    logging.debug("DONE MAKING INPUT FILES FROM rMATS FILE")
    logging.debug("=" * 50)
    
    # Generate RNA map
    logging.debug("=" * 50)
    logging.debug("GENERATING RNA MAP")
    
    rna_map_script = bin_path / 'RNA.map.noWiggle.A5SS.py'
    
    # Build command
    cmd = [
        sys.executable,
        str(rna_map_script),
        str(exon_path),
        str(peak_file_path),
        str(args.intron),
        str(args.exon),
        str(args.window),
        str(args.step),
        str(args.sigFDR),
        args.label,
        str(out_path),
        str(nu),
        str(nd),
        str(nb),
        str(int(args.separate)),
        'A5SS'
    ]
    
    logging.debug("Running RNA map with command: %s", ' '.join(cmd))
    
    import subprocess
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    logging.debug("RNA map was generated with status %d", result.returncode)
    
    if result.returncode != 0:
        logging.error("Error in generating RNA map: %d", result.returncode)
        logging.error("Error detail: %s", result.stderr)
        sys.exit(-9)
    
    logging.debug("RNA map output:\n%s", result.stdout)
    logging.debug("DONE GENERATING RNA MAP")
    logging.debug("=" * 50)
    
    # Calculate runtime
    end_time = time.time()
    runtime = end_time - start_time
    logging.debug("Program ended")
    logging.debug("Program ran %02d:%02d:%02d",
                  int(runtime / 3600),
                  int((runtime % 3600) / 60),
                  int(runtime % 60))
    
    print(f"CLIP-seq A5SS analysis completed successfully!")
    print(f"Results written to: {out_path}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
