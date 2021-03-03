#!/usr/bin/env python3

import logging
import argparse
import configargparse
import pathlib
import sys

from Coreugate import RunCoreugate


def run_pipeline(args):
    '''
    Run the pipeline
    '''
    C = RunCoreugate(args)
    return(C.run_pipeline())

def set_parsers():
    # setup the parser
  
    parser = configargparse.ArgumentParser(description='Coreugate - a cgMLST pipeline implementing chewBACCA',formatter_class=configargparse.ArgumentDefaultsHelpFormatter)
    # use add_argument
    parser.add_argument('--input_file',
        '-i', 
        help='Input file tab-delimited file 2 or 3 columns isolate_id path_to_input_file(s)',
        default = '')
    parser.add_argument('--schema_path',
        '-sp',
        help='Path to species schema/allele db',
        default = '')
    parser.add_argument('--min_contig_size', 
        '-c', 
        help='Minumum contig size required for QC', 
        default=500)
    parser.add_argument('--min_contigs', 
        '-m', 
        help='Minumum number of contigs required for QC', 
        default=0)
    parser.add_argument('--assembler', 
        '-a', 
        help='Assembler to be used (options are: shovill-spades, shovil-skesa, shovill-velvet, skesa, spades)', 
        default= 'shovill-spades')
    parser.add_argument('--prodigal_training', 
        '-p', 
        help='Prodigal file to be used in allele calling. See https://github.com/B-UMMI/chewBBACA/tree/master/CHEWBBACA/prodigal_training_files for options', 
        default= '')
    parser.add_argument(
        "--singularity",
        "-S",
        action="store_true",
        help="If using singularity container for chewBACCA"
    )
    parser.add_argument(
        "--singularity_path",
        "-s",
        default="",
        help="Path to the singularity container for chewBACCA"
    )
    parser.add_argument(
        "--workdir",
        "-w",
        default=f"{pathlib.Path.cwd().absolute()}",
        help="Working directory, default is current directory",
    )
    parser.add_argument(
        "--resources",
        "-r",
        default=f"{pathlib.Path(__file__).parent }",
        help="Directory where templates are stored",
    )
    
    parser.add_argument('--threads',
        '-t',
        help='Number of threads to run chewBACCA', 
        default=1
    )
    
    parser.set_defaults(func=run_pipeline)
    args = parser.parse_args()
    
    if vars(args) == {}:
        parser.print_help(sys.stderr)
    else:
        
        args.func(args)
	

def main():
    """
    run pipeline
    """

    args = set_parsers()
    
    # if vars(args) == {}:
    #     parser.print_help(sys.stderr)
    # else:
    # print(logging.__file__)
    # logger = logging.getLogger(__name__)
    # logger.setLevel(logging.INFO)
    # # # create file handler which logs even debug messages
    # # # create file handler which logs even debug messages
    # fh = logging.FileHandler('coreugate.log')
    # fh.setLevel(logging.INFO)
    # # create console handler with a higher log level
    # ch = logging.StreamHandler()
    # ch.setLevel(logging.INFO)
    # # create formatter and add it to the handlers
    # formatter = logging.Formatter('[%(levelname)s:%(asctime)s] %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    # fh.setFormatter(formatter)
    # ch.setFormatter(formatter)
    # # add the handlers to the logger
    # logger.addHandler(ch)
    # logger.addHandler(fh)
    # # logger.info(f"Ran Coreugate with : {' '.join(sys.argv)}")
    # # args.func(args)


if __name__ == "__main__":
    main()
