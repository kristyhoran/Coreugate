import pathlib
import os
import pandas
import jinja2
import logging
import filecmp
import datetime
import numpy
import itertools
import subprocess
import pathlib
import re
# from Bio import SeqIO, Phylo
from packaging.version import Version


class RunCoreugate:


    def __init__(self,args):
        
        # set up logger
        self.logger =logging.getLogger(__name__) 
        self.logger.setLevel(logging.INFO)
        fh = logging.FileHandler('coreugate.log')
        fh.setLevel(logging.INFO)
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        formatter = logging.Formatter('[%(levelname)s:%(asctime)s] %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p') 
        ch.setFormatter(formatter)
        fh.setFormatter(formatter)
        self.logger.addHandler(ch) 
        self.logger.addHandler(fh)

        # set up files
        if args.input_file == '':
            self.logger.warning(f"input_file can not be empty. Please check your inputs and try again.") 
            raise SystemExit
        elif args.input_file != '':
            self.input_file = self._check_file(args.input_file)
        if args.schema_path == '':
            self.logger.warning(f"schema_path can not be empty. Please check your inputs and try again.") 
            raise SystemExit
        else:
            self.schema_path = self._check_file(args.schema_path) 
        self.prodigal_training = self.check_prodigal(args.prodigal_training)
        self.workdir = self._check_file(args.workdir)
        self.threads = args.threads
        self.force = args.force
        self.cluster = args.cluster
        if self.cluster:
            # check that thresholds are applied
            self.cluster_thresholds = self.check_cluster_thresholds(args.cluster_thresholds)
        else: 
            self.cluster_thresholds = args.cluster_thresholds
        self.filter_threshold = self.check_filter_threshold(args.filter_samples_threshold)
        self.report = args.report
        
    def check_prodigal(self, ptf):

        if ptf == '':
            return ''
        elif pathlib.Path(ptf).exists():
            return pathlib.Path(ptf).absolute(  )
        else:
            self.logger.critical(f"There seems to be a problem with you prodigal training file. Please try again.")
            raise SystemExit

    def check_filter_threshold(self, threshold):
        try:
            t = float(threshold)
            if 0 <= t <=1.0:
                return t
            else: 
                self.logger.critical(f"Proportion must be between 0 and 1. Please try again.")
        except ValueError:
                self.logger.critical(f"There seems to be a problem with your proortion. Please try again.")
                raise SystemExit

    def check_cluster_thresholds(self, thresholds):
        if thresholds.split(',') == '':
            self.logger.critical(f"If you wish to cluster the pairwise allelic distance matrix you will need to supply at least on threshold.")
            raise SystemExit
        else:
            print(thresholds)
            for t in thresholds.split(','):
                try:
                    isinstance(int(t), int)
                    
                    
                except ValueError:
                    self.logger.critical(f"There seems to be a problem with your thresholds. Please try again.")
                    raise SystemExit
                except TypeError:
                    self.logger.critical(f"There seems to be a problem with your thresholds. Please try again.")
                    raise SystemExit
            return thresholds

    def _check_file(self, path):
        '''
        check version of software installed, if chewbbaca add string
        :software: name of software (str)
        '''
        if path == '':
            self.logger.warning(f"{path} can not be empty. Please check your inputs and try again.")
            raise SystemExit
        elif pathlib.Path(path).exists() and os.access(path, os.R_OK):
            return pathlib.Path(path)
        else:
            self.logger.warning(f"{path} can not be found or is not accessible. Please check your inputs and try again.")
            raise SystemExit

    def check_version(self, software):
        '''
        check version of software installed, if chewbbaca add string
        :software: name of software (str)
        '''
        # self.logger.info(f"Checking that {software} is installed and recording version.")
        version_pat = re.compile(r'\bv?(?P<major>[0-9]+)\.(?P<minor>[0-9]+)\.(?P<release>[0-9]+)(?:\.(?P<build>[0-9]+))?\b')
        try:
            cmd = f"{software} --version 2>&1"
            self.logger.info(f"Checking that {software} is installed and recording version : {cmd}")
            p = subprocess.run(cmd, shell = True, capture_output = True, encoding = 'utf-8')
            # print(sft)
            sft = p.stdout.split()[-1]
            # print(type(sft))
            # print(version_pat.match(sft))
            sft_version = version_pat.match(sft).group()
            # pritn()
            self.logger.info(f"{software} {sft_version} found. Good job!")
            return sft_version
        except FileNotFoundError:
            logger.info(f"{software} is not installed.")
            raise SystemExit
    


    def check_chewbbaca(self):
        '''
        check chewbbaca version -> needs to be 2.0.16
        '''
        # self.logger.info(f'Checking if chewBACCA is installed.')
        try:
            chewie_version = self.check_version('chewBBACA.py')
            base_version = Version("2.0.16")
            if Version(chewie_version) >= base_version:
                self.logger.info(f"chewBBACA version {chewie_version} has been found.")
            else:
                self.logger.warning(f"chewBACCA {chewie_version} has been found. This is not compatible with coreugate, please install chewBBACA version >= 2.0.16 before proceeding (https://github.com/B-UMMI/chewBBACA). Or use `-s Y`. Exiting....")
                raise SystemExit
        except FileNotFoundError:
            self.logger.warning(f"chewBBACA is not installed. please install chewBBACA <= version 2.0.16 before proceeding (https://github.com/B-UMMI/chewBBACA). Exiting....")
            raise SystemExit
    

    def run_checks(self):
        
        self.check_chewbbaca()
        # self.check_assemblers()

    def check_input_file(self, tab):
        
        if tab.shape[1] != 2:
            self.logger.critical(f"Your input file is not in the correct format. The input file should be a tab-delimited file, column 1 is isolate ID and column 2 is path to contigs, with no header line. Exiting....")
            raise SystemExit
        return True
        
    def link_inputs(self, data_source, isolate_id):
        '''
        check if read source exists if so check if target exists - if not create isolate dir and link. If already exists report that a dupilcation may have occured in input and procedd

        '''
        # check that job directory exists
        # logger.info(f"Checking that reads are present.")
        
        if data_source.exists() and os.access(data_source, os.R_OK):
            I = pathlib.Path(f"{isolate_id}") # the directory where contigs will be stored for the isolate
            if not I.exists():
                I.mkdir()
            data_target = I / f"{isolate_id}.fa"
            if not data_target.exists():
                subprocess.run(f"ln -sf {data_source} {data_target}", shell = True)
                # data_target.symlink_to(data_source)
        else:
            logger.warning(f"{data_source} does not seem to a valid path. Please check your input and try again.")
            raise SystemExit()
    
    def check_inputs_exists(self, tab):
        '''
        check that the correct paths have been given
        if reads not present path_exists will cause a FileNotFound error and warn user
        :input
            :tab: dataframe of the input file
        '''
        self.logger.info(f"Checking that all the data files exist.")
        # pos_1 = 'R1.fq.gz' if self.input_type == 'READS' else 'contigs.fa'
        paths = []
        for i in tab.iterrows():
            # print(i[1][0])
            if not '#' in i[1][0]:
                c = i[1][1]
                islt = f"{i[1][0].strip()}"
                self._check_file(c)
                self.link_inputs(pathlib.Path(c), isolate_id=islt)
                ctg = f"{islt}.fa"
                paths.append(f"{pathlib.Path(islt,ctg)}")
            
        return paths

    
    def prep_external_schema(self):
        '''
        check if the schema has been prepped
        :path:
            path is path to schema
        :cpu:
            cpu for running PrepExternalSchema (threads option)
        '''
        short_path= self.schema_path / 'short'
        # if no short directory run PrepExternalSchema
        if not os.path.exists(short_path):
            self.logger.info(f'Preparing external schema. This may take some time - feel free to go get coffee/tea/refreshment of choice.')
            subprocess.call(['chewBBACA.py', 'PrepExternalSchema', '-i', f"{self.schema_path}", '--cpu', f"{self.threads}", '-v' ])
        else:
            self.logger.info(f'Schema ({self.schema_path}) is already in the correct format. Congratulations.')
        # if short/fasta is present and empty then remove it to prevent problems with chewbbaca
        fasta = self.schema_path.glob("*.fasta")
        self.logger.info('Checking for abberent fasta files. Incorrectly formatted fasta files and empty fasta files will be removed.')
        for f in fasta:
            if os.path.exists(short_path / f / '_short.fasta') and os.path.getsize(short_path / f / '_short.fasta') == 0:
                subprocess.call(['rm', short_path / f / '*'])
                subprocess.call(['rm', self.schema_path / f / '*'])
        # if temp exists remove it
        if os.path.exists(self.schema_path / 'temp'):
            subprocess.call(['rm', '-r',  str(self.schema_path/ 'temp')])
        
   
    def link_schema(self):

        target = pathlib.Path(self.schema_path.name)
        source = self.schema_path
        self.prep_external_schema()

        if not target.exists():
            target.symlink_to(source)


    def setup_working_directory(self):
        '''
        ensure all required files are linked to the working directory, job_directory and data directory
        :input_path: input_file file path
        :workdir: working directory
        :day: day for log
        :job_id: unique identifer will be used to make the job_directory
        :schema: path to schema
        :cpu: cpu for running PrepExternalSchema if needed
        '''
        # make the job directory
        self.logger.info(f"Setting up schema for cgMLST")
        self.link_schema()
               
        # set up data directory                
        # link input data to input directory
        self.logger.info(f"Linking source data to working directory")
        tab = pandas.read_csv(self.input_file, sep = '\t', header = None)
        self.check_input_file(tab = tab)
        # print(tab)
        self.contigs = self.check_inputs_exists(tab = tab)
        # generate a generic input file for chewie
        self.logger.info(f"Writing an input file for chewie")
        o = pathlib.Path('contigs.txt')
        o.write_text('\n'.join(self.contigs))

    def write_workflow(self):
        '''
        write the Snakefile and config.yaml
        '''
        
        self.logger.info(f"Setting up specific workflow")

        cfg = { 'input': 'contigs.txt',
                'schema_path': f"{self.schema_path}",
                'ptf':f"{self.prodigal_training}",
                'chewie_threads': self.threads,
                'force': 'false',
                'filter_threshold': self.filter_threshold,
                'cluster': f"{self.cluster}".lower(),
                'cluster_threshold':self.cluster_thresholds
        }        
        # read the config file which is written with jinja2 placeholders (like django template language)
        template_dir = f"{pathlib.Path(__file__).parent / 'utils'}"
        config_template = jinja2.Template(pathlib.Path(template_dir, 'config.j2.yaml').read_text())
        config = self.workdir /'nextflow.config'
        
        config.write_text(config_template.render(cfg))
        self.logger.info(f'this is the schema path : {self.schema_path}')
        self.logger.info(f"Config file successfully created")
        
    def setup_pipeline(self):
        
        self.setup_working_directory()
        self.write_workflow()

    def run_workflow(self):
        '''
        run Coreugate
        set the current directory to working dir for correct running of pipeline if singularity = Y then run with singularity
        if the pipeline wroks, return True else False
        '''
        resume = f"" if self.force else "-resume"
        report = "-with-report coreugate_resource_report.html -with-trace"
        extra_args = f"{resume} {report}"
        cmd = f"nextflow {pathlib.Path(__file__).parent / 'utils' / 'main.nf'} {extra_args}"
        self.logger.info(f"Running {cmd} - patient you must be.")
        wkf = subprocess.run(cmd, shell = True)
        while True:
            if wkf.stdout != None:
                line = wkf.stdout.readline().strip()
                if not line:
                    break
            line = ''
            break
            self.self.logger.info(f"{line}")
        if wkf.returncode == 0:
            return True
        else:
            return False
        
    
    def run_pipeline(self):

        self.run_checks()
        self.setup_pipeline()
        self.run_workflow()
        
