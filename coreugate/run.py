#!/usr/bin/env python3

import logging
import argparse
import pathlib
import sys

from coreugate.coreugate import RunCoreugate

VERSION = '2.0.0'

def run_pipeline(args):
    '''
    Run the pipeline
    '''
    C = RunCoreugate(args)
    return(C.run_pipeline())

def set_parsers():
    # setup the parser
  
    parser = argparse.ArgumentParser(description='Coreugate - a cgMLST pipeline implementing chewBACCA',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + VERSION)
    parser.add_argument('--input_file',
        '-i', 
        help='Input file tab-delimited file3 columns isolate_id path_to_input_file (contigs)',
        default = '')
    parser.add_argument('--schema_path',
        '-s',
        help='Path to species schema/allele db (or url if using chewie Nomenclature server)',
        default = '')
    
    parser.add_argument('--prodigal_training', 
        '-p', 
        help='Prodigal file to be used in allele calling. See https://github.com/B-UMMI/chewBBACA/tree/master/CHEWBBACA/prodigal_training_files for options', 
        default= '')
    parser.add_argument(
        "--workdir",
        "-w",
        default=f"{pathlib.Path.cwd().absolute()}",
        help="Working directory, default is current directory",
    )   
    parser.add_argument('--threads',
        '-t',
        help='Number of threads to run chewBACCA', 
        default=16
    )
    parser.add_argument(
        '--filter_samples_threshold',
        '-ft',
        default= 0.95,
        help = f"The proportion of loci present in a sample for an sample to be included in further analysis (0-1)"
        )
    parser.add_argument(
        '--cluster',
        '-c',
        action = "store_true", help='If you would like to cluster the pairwise distance matrix. If selected you must provide a list of thresholds.' 
    )
    parser.add_argument(
        '--cluster_thresholds',
        '-ct',
        default='',
        help="Provide a comma separate list (NO SPACES) eg 20,40,200"
    )
    parser.add_argument(
        '--force',
        '-f',
        action = "store_true", help='If you want to force chewBBACA to re-run.'
        )
    parser.add_argument(
        '--report',
        help = "Save nextflow reports.",
        action="store_true"
    )
    parser.set_defaults(func=run_pipeline)
    args = parser.parse_args()
    
    if vars(args) == {}:
        parser.print_help(sys.stderr)
    else:
        print(args)
        args.func(args)
	

def main():
    """
    run pipeline
    """

    args = set_parsers()
    

if __name__ == "__main__":
    main()
