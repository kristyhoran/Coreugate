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
from Bio import SeqIO, Phylo
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
        else:
            self.input_file = self._check_file(args.input_file) 
        self.input_type = self._input_type()
        if args.schema_path == '':
            self.logger.warning(f"schema_path can not be empty. Please check your inputs and try again.") 
            raise SystemExit
        else:
            self.schema_path = self._check_file(args.schema_path) 
        self.min_contig_size = args.min_contig_size
        self.min_contigs = args.min_contigs
        self.assembler = args.assembler
        self.prodigal_training = args.prodigal_training
        self.run_with_singularity = args.singularity
        self.singularity_path = args.singularity_path
        self.workdir = self._check_file(args.workdir)
        self.resources = args.resources
        self.threads = args.threads

    def _check_file(self, path):
        '''
        check version of software installed, if chewbbaca add string
        :software: name of software (str)
        '''
        if path == '':
            self.logger.warning(f"{path} can not be empty. Please check your inputs and try again.")
            raise SystemExit
        elif pathlib.Path(path).exists():
            return pathlib.Path(path)
        else:
            self.logger.warning(f"{path} can not be found. Please check your inputs and try again.")
            raise SystemExit

    def _input_type(self):
        """
        check the dimensions of the input file, if == 2 then inputs is contigs, if == 3 the reads
        """
        tab = pandas.read_csv(self.input_file, sep = '\t')
        if tab.shape[1] == 2:
            self.logger.info(f"You are running COREugate with from assemblies.")
            return "CONTIGS"
        elif tab.shape[1] == 3:
            self.logger.info(f"You are running COREugate with from reads.")
            return "READS"
        else:
            self.logger.warning(f"Sorry but your input file is incorrectly foramtted, delimiter should be a tab and it should contain only 2 or 3 columns (contigs or reads respoectively). Please check and try again.")
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
                self.logger.warning(f"chewBACCA {chewie_version} has been found. This is not compatible with COREugate, please install chewBBACA version >= 2.0.16 before proceeding (https://github.com/B-UMMI/chewBBACA). Or use `-s Y`. Exiting....")
                raise SystemExit
        except FileNotFoundError:
            self.logger.warning(f"chewBBACA is not installed. please install chewBBACA <= version 2.0.16 before proceeding (https://github.com/B-UMMI/chewBBACA). Or use `-s Y`. Exiting....")
            raise SystemExit
    

        
    def check_assemblers(self):
        '''
        check assemblers
        '''
        self.logger.info(f'Checking that {self.assembler} is installed.')
        tocheck = 'shovill' if 'shovill' in self.assembler else self.assembler
        asmb_version = self.check_version(tocheck)
        
    def run_checks(self):
        
        if self.run_with_singularity:
            self.logger.info(f"You are using singularity containers. Good luck")           
        else:
            self.check_chewbbaca()
        
        if not self.run_with_singularity and self.input_type == 'READS':
            self.check_assemblers()


    def link_inputs(self, data_source, isolate_id, data_name):
        '''
        check if read source exists if so check if target exists - if not create isolate dir and link. If already exists report that a dupilcation may have occured in input and procedd

        '''
        # check that job directory exists
        # logger.info(f"Checking that reads are present.")
        
        R = self.workdir / f"{self.input_type}"
        if not R.exists():
            R.mkdir()
        
        # if f"{read_source}"[0] != '/':
        #     read_source = self.workdir / read_source
        
        if data_source.exists():
            I = R / f"{isolate_id}" # the directory where reads will be stored for the isolate
            if not I.exists():
                I.mkdir()
            data_target = I / f"{data_name}"
            if not data_target.exists():
                subprocess.run(f"cp {data_source} {data_target}", shell = True)
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
        pos_1 = 'R1.fq.gz' if self.input_type == 'READS' else 'contigs.fa'
        for i in tab.iterrows():
            # print(i[1][0])
            if not '#' in i[1][0]:
                r1 = i[1][1]
                self._check_file(r1)
                self.link_inputs(pathlib.Path(r1), isolate_id=f"{i[1][0].strip()}", data_name=pos_1)
                if self.input_type == 'READS':
                    r2 = i[1][2]
                    self._check_file(pathlib.Path(r2))
                    self.link_reads(pathlib.Path(r2), isolate_id=f"{i[1][0].strip()}", data_name='R2.fq.gz')
        return True

    
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
        # print(tab)
        self.check_inputs_exists(tab = tab)
        self.isolates = ' '.join(list(tab.iloc[ : , 0]))

    def write_workflow(self):
        '''
        write the Snakefile and config.yaml
        '''
        if self.prodigal_training != '':
            ptf_string = f"--ptf {self.prodigal_training}"
        else:
            ptf_string = ''
        self.logger.info(f"Setting up specific workflow")        
        # read the config file which is written with jinja2 placeholders (like django template language)
        template_dir = f"{pathlib.Path(self.resources, 'templates')}"
        config_template = jinja2.Template(pathlib.Path(template_dir, 'config.yaml').read_text())
        config = self.workdir /'config.yaml'
        
        config.write_text(config_template.render(isolates = self.isolates, data_type = self.input_type, min_contig_size = self.min_contig_size, min_contigs = self.min_contigs,cpu=self.threads, schema_path = self.schema_path,workdir = self.workdir, ptf = ptf_string, assembler = self.assembler, script_path = template_dir, input_file = f"{self.input_file}"))
        self.logger.info(f'this is the schema path : {self.schema_path}')
        self.logger.info(f"Config file successfully created")
        # template_dir = f"{pathlib.Path(self.resources, 'templates')}"
        # rfile_template = jinja2.Template(pathlib.Path(template_dir, 'distances.R').read_text())
        # rfile = 'distances.R'
        # rfile.write_text(rfile_template.render(script_dir = template_dir))

    def setup_pipeline(self):
        
        self.setup_working_directory()
        self.write_workflow()

    def run_workflow(self):
        '''
        run Coreugate
        set the current directory to working dir for correct running of pipeline if singularity = Y then run with singularity
        if the pipeline wroks, return True else False
        '''
        
        if self.run_with_singularity:
            cmd = f"snakemake -s {pathlib.Path(self.resources, 'templates', 'Snakefile')} --use-singularity --singularity-args '--bind /home' -d {self.workdir} -j 1"
        else:
            cmd = f"snakemake -s {pathlib.Path(self.resources, 'templates', 'Snakefile')} -d {self.workdir} -j 1"
        self.logger.info(f"Running {cmd} - patient you must be.")
        wkf = subprocess.run(cmd, shell = True)
        if wkf.returncode == 0:
            return True
        else:
            return False
        
    def finish_workflow(self):
        '''
        final message at completion of workflow. If workflow goes to completion print 'thanks for coming message'
        '''
        output_stats = self.workdir / 'overall_statistics.tsv'
        output_alleles = self.workdir / 'overall_alleles.tsv'
        self.logger.info(f'Checking that output files are present.')
        if output_stats.exists() and output_alleles.exists():
            
            self.logger.info(f"COREugate has finished.")
            self.logger.info(f"May the force be with you.") 
        else:
            self.logger.warning(f"Something has gone wrong, please check logs and try again. If problem persists, contact developer.")
        
    def run_pipeline(self):

        self.run_checks()
        self.setup_pipeline()
        self.run_workflow()
        self.finish_workflow()
