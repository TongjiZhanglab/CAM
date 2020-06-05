#!/usr/bin/python
# Time-stamp: <2016-12-25 Shengen Hu>
"""
 <CAM: a quality control and analysis pipeline for MNase-seq data>
    Copyright (C) <2017>  <Shengen Hu>
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
# ------------------------------------
# Python Modual
# ------------------------------------

import os
import sys
import time
import string
import argparse
import subprocess

# -----------------------------------
# custom package
# -----------------------------------
import CAMpipe

### tool function
from CAMpipe.Utility      import          (sp,
                                   pdf_name,
                                   raise_error,
                                   wlog,
                                   readAnnotation,
                                   textformat,
                                   CMD,
                                   check_filelist
                                   )
### read and generate config file
from CAMpipe.parse_config import (gen_conf,
                                   read_conf,
                                   make_conf,
                                   )     


                                   
# -------------------
# main step
# -------------------
from CAMpipe.step0_integrate_data   import step0_integrate_data
from CAMpipe.step1_preprocess       import step1_preprocess
from CAMpipe.step2_QC    import step2_QC
from CAMpipe.step3_nucarray import step3_nucarray
from CAMpipe.step5_summary import step5_summary

# ------------------------------------
# Misc functions
# ------------------------------------

    
#### read options

class ChiLinParser(argparse.ArgumentParser):
    '''
    Ex version for argparse(parameter parser) , add raise error function .
    '''
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit()

def parse_args():
    '''
    Read parameter 
    '''
    description = "CAM -- a quality control and analysis pipeline for MNase-seq data"
    parser = ChiLinParser(description = description, version = "CAM 1.0.0 20170102")
    sub_parsers = parser.add_subparsers(help = "sub-command help", dest = "sub_command")

    ### generate config file
    template_parser = sub_parsers.add_parser("gen",  help = "generate a template of config file",
                                             description = "CAM config file generation. Usage: CAM.py gen  -n config_name.conf")
    template_parser.add_argument("-n","--name", dest="config_name",required = True,help="name of your config file : config_name.conf")

    ### run config file
    pipe_parser = sub_parsers.add_parser("run", help = "run pipeline with a config file input",
                                         description = "Run CAM pipeline with a config file input")
    pipe_parser.add_argument("-c","--config", required = True,
                             help = "specify the config file, -c config_name.conf" )
    pipe_parser.add_argument("-f","--force_overwrite",dest='fover',  default=False, action='store_true', 
                             help = "The pipeline will over write output result if the output folder is already exist " )
    pipe_parser.add_argument("--clean",dest='Clean' , default=False, action='store_true',
                             help = "remove intermediate result generated during CAM,default is NO" )
    ### simple mode
    simple_parser = sub_parsers.add_parser("simple", help = "run CAM using simple mode",
                                         description = "(Run CAM pipeline using simple mode/command line mode) \
                                                        Usage: CAM.py -a raw_Reads.fastq -t SE -n outname -s hg19 \
                                                        --fa /home/user/data/hg19.fa")
    simple_parser.add_argument("-a","--inputa", dest = 'inputa',required = True,
                             help = "[required] main input file of CAM, accept fastq (raw MNase-seq data), \
                                     sam/bam (bowtie alignment output) and bed (alignment output transformed from sam, \
                                     6 column, 5th column is bowtie map quality and 6th column is strand ) format , \
                                     format is decided by the extension of input file(.fastq, .sam or .bed) " )
    simple_parser.add_argument("-b","--inputb",dest='inputb',required = False,
                             help = "This pararmeter is for the 2nd part of paired-end fastq files. \
                                     Required only if input file is fastq file and sequencing type is paired-end sequencing, \
                                     only accept raw fastq input. The extension should be .fastq" )
    simple_parser.add_argument("-c","--customregion",dest='customregion',required = False,
                             help = "[optional, absolute path]user defined cutsom region to display nucleosome pattern. \
                                     CAM will add an additional nucleosome profile on input region if this parameter is specified. \
                                     Only accept bed format input (3 column or 6 column bed). The extension should be .bed. \
                                     if 6 column bed and the 6th column is +/-, CAM will display strand specific signal. \
                                     limit of custom region is 50k, if the number of custom region is greater than 50k, CAM will take top 50k regions for nucleosome profile")
    simple_parser.add_argument("-t","--seqtype", dest="seqtype",default = "SE",choices=["SE","PE"],
                             help="sequencing type of your MNase-seq sample, choose from SE(single end, default) and PE(paired end)")
    simple_parser.add_argument("-n","--name", dest="name",required = True,
                             help="[required] name of you config file and output dir, name only , no extension. \
                                   The output files will be named like name.pdf, name.txt ... ")
    simple_parser.add_argument("-s","--species",dest='species', required = False,
                             help = "[required] genome version name:  \
                                     CAM require genome version to find corresponded annotation file. by default we support hg19,hg38,mm9,mm10 genome version,\
                                     users can simply input, eg: [-s hg19]  to specify genome version. \
                                     users can also add other genome versions by adding annotation files in the package (see documents for detail. \
                                     NOTE: this parameter should correspond with --fa/--mapindex/--genome2bit, in other word, make sure you are using same genome version.")
    simple_parser.add_argument("--fa",dest='fa',required = False,
                             help = "[required, absolute path] aka genome_fasta: whole genome sequence in fasta(.fa) or 2bit(.fa) format,\
                                     input file can only be fasta and 2bit format, and file name should end with .fa or .2bit ,\
                                     CAM will build mapindex (if need mapping) before process. ")
    simple_parser.add_argument("--mapindex",dest='mapindex',required = False,
                             help = "mapping index folder, user can add input mapping index directly to skip bowtie-build step, which is time consuming\
                                     The index file should named like hg19.1.ebwt (should be end with .ebwt). \
                                     this parameter should be absolute path of index filename prefix \
                                     (without trailing .X.ebwt).  eg: --mapindex /Users/mapping_index/hg19.bowtie/hg19 \
                                     (then under your folder  /Users/mapping_index/hg19.bowtie/  \
                                     there should be hg19.1.ebwt, hg19.rev.1.ebwt ....  , see documents for detail)" )
    simple_parser.add_argument("-p","--thread",dest='P' ,default='8', 
                             help = "number of alignment threads for mapping step, ignored for sam/bed input" )
    simple_parser.add_argument("-f","--force_overwrite",dest='fover',  default=False, action='store_true', 
                             help = "specify the config file to create output folder , this cmd will remove existing result if set True ~!! " )
    simple_parser.add_argument("--clean",dest='Clean' , default=False, action='store_true',
                             help = "remove intermediate result generated during CAM, default is No" )
    simple_parser.add_argument("--task",dest='task' , default=4,choices=['1','2','3','4'], required = False,
                             help="[required] to proacess which process")
    # simple_parser.add_argument("--start",dest='start' , default=0, required = False,
                             # help="[required] to starta from which process")
    args = parser.parse_args()
    ## generate config file template 
    if args.sub_command == "gen":
        gen_conf(args.config_name)
        sys.exit(0)

    ## run CAM with config file input
    if args.sub_command == "run":
        if os.path.isfile(args.config):
            return args
        else:
            print ('ERROR : -c input is not a config file\n')
            print pipe_parser.print_help()
            sys.exit()
    
    ## run CAM with a simple mode, input parameter in command line
    if args.sub_command == "simple":
        if args.name.endswith('.conf'):
            args.name = args.name[:-5]
        make_conf(args.inputa,args.inputb,args.seqtype,args.name,args.customregion,args.fover,args.species,args.P,args.fa,args.mapindex)
        args.config = args.name + '.conf'
        return args

# ------------------------------------
# Main function
# ------------------------------------

def main():

    args = parse_args()
    # print int(args.task)
    # print int(args.task) >=2
    # print int(args.task) >=3
    # print int(args.task) ==4
    # sys.exit(1)

    # print 'Run only %s tasks' % args.task
    # sys.exit(1)

    conf_dict = read_conf(args.config)
    ### read raw path of output dir, the startdir will be used when the input file is not in absolute path
    conf_dict['General']['startdir'] = os.getcwd()+'/'

    ### check output name and dir from input parameter
    if conf_dict['General']['outname'] == "":
        print 'your outname cannot be left blank,exit'
        sys.exit(1)
    if "." in conf_dict['General']['outname']:
        oldname = conf_dict['General']['outname']
        newname = oldname.replace(".","-")
        conf_dict['General']['outname'] = newname
        print 'replace outname from %s to %s for latex summary'%(oldname,newname)
    if conf_dict['General']['outputdirectory'] == "":
        conf_dict['General']['outputdirectory'] = conf_dict['General']['outname']
        print 'output directory is blank, use outname as directory name and set output directory in current folder'
    if "~" in conf_dict['General']['outname']:
        print 'ERROR: ~ cannot appeared in outname, current outname is %s'%(conf_dict['General']['outname'])
        sys.exit(1)
    if "~" in conf_dict['General']['outputdirectory']:
        print 'ERROR: require absolute path for outputdirectory'
        sys.exit(1)
    if not conf_dict['General']['outputdirectory'].endswith('/'):
        conf_dict['General']['outputdirectory'] += '/'
    if not conf_dict['General']['outputdirectory'].startswith('/'):
        conf_dict['General']['outputdirectory'] = conf_dict['General']['startdir'] + conf_dict['General']['outputdirectory']
    
    ### creat output dir
    if os.path.isfile(conf_dict['General']['outputdirectory'].rstrip("/")):
        print 'ERROR: name of your output dir %s is exist as a file, cannot create a dir,CAM exit'%(conf_dict['General']['outputdirectory'].rstrip("/"))
        sys.exit(1)
    elif os.path.isdir(conf_dict['General']['outputdirectory']):
        if not args.fover:
            print 'name of your output dir %s is exist as a dir, overwrite function (-f) is turned off, write output result in existing dir'%(conf_dict['General']['outputdirectory'])
        else: 
            print 'name of your output dir %s is exist as a dir, overwrite function (-f) is turned on, write output result in existing dir'%(conf_dict['General']['outputdirectory'])
    else:
        os.system("mkdir %s"%(conf_dict['General']['outputdirectory']))
    
    ### move to output dir
    os.chdir(conf_dict['General']['outputdirectory'])
    ## cp config file to output folder
    if args.config.startswith('~') or args.config.startswith('/'):
        cmd = 'cp %s .'%(args.config)
    else:
        cmd = 'cp %s .'%(conf_dict['General']['startdir']+args.config)
    CMD(cmd)
    ### specify the main progress log file
    logfile = conf_dict['General']['outputdirectory']+'progress_log.txt'
    filelist = conf_dict['General']['outputdirectory']+'progress_filelist.txt'
    ## remove existing log file. 
    if os.path.isfile(logfile) and args.fover:
        CMD('rm %s'%logfile)
    if os.path.isfile(filelist) and args.fover:
        CMD('rm %s'%filelist) 
    wlog('',filelist)   
    ### Rscript location 
    conf_dict['rscript'] = os.path.join(CAMpipe.__path__[0], "Rscript/")    
    conf_dict['clean'] = args.Clean
    
    ### default annotation location
    conf_dict['default_anno_dir'] = os.path.join(CAMpipe.__path__[0], "annotation/")
    
    ### main step for CAM, see individual script for detail note.
    # preparing step, integrate parameter, prepare for following step
    # wlog(args.task,logfile)
    wlog("Step0: preparation",logfile)
    t = time.time()
    step0_integrate_data(conf_dict,logfile)
    step0time = time.time() -t
    wlog("running time for Step0: %s"%(step0time),logfile)

    
    # main data processing step, including mapping, generate bigwig file
    wlog("Step1: data-preprocess",logfile)
    final_file = conf_dict['General']['outputdirectory'] + 'preprocess/' + conf_dict['General']['outname'] + '.bb'
    final_file_bp = conf_dict['General']['outputdirectory'] + 'summary/' + conf_dict['General']['outname'] + '.bb'
    # print not args.fover and (os.path.isfile(final_file) or os.path.isfile(final_file_bp))
    # sys.exit(1)
    if not args.fover and (os.path.isfile(final_file) or os.path.isfile(final_file_bp)) :
        print 'Step1 has been finished.'
    else:
        t = time.time()
        step1_preprocess(conf_dict,logfile,filelist)
        step1time = time.time() -t
        wlog("running time for Step1: %s"%(step1time),logfile)

    # nucarray step, including  nucarray detection, annotation of detected nucarray
    if int(args.task) >=2:
        wlog("Step2: call well-positioned nucleosome array (analysis)",logfile)
        final_file = conf_dict['General']['outputdirectory'] + 'nucarray/' + conf_dict['General']['outname'] + '_center_position.bw'
        final_file_bp = conf_dict['General']['outputdirectory'] + 'summary/' + conf_dict['General']['outname'] + '_center_position.bw'
        if not args.fover and (os.path.isfile(final_file) or os.path.isfile(final_file_bp)):
            print 'Step2 has been finished.'
        else:
            t = time.time()
            step3_nucarray(conf_dict,logfile,filelist)
            step3time = time.time() -t
            wlog("running time for Step2: %s"%(step3time),logfile)

    # QC and analysis step, including QC and analysis
    if int(args.task) >=3:
        wlog("Step3: QC and nucleosome profile",logfile)
        final_file = conf_dict['General']['outputdirectory'] + 'QC/' + conf_dict['General']['outname'] + '_Features.txt'
        file_count_cmd = 'ls %s'%(conf_dict['General']['outputdirectory'] + 'QC/')
        files = sp(file_count_cmd)[0].split('\n')
        files_num = 0
        for i in range(len(files)):
            if files[i].endswith('.pdf'):
                files_num += 1
        # print 'Total pdf files:%s' % files_num
        # print os.path.isfile(final_file) and files_num==8
        # sys.exit(1)
        # print os.path.isfile(final_file)
        if not args.fover and ((os.path.isfile(final_file) and files_num==8)):
            print 'Step3 has been finished.'
        else:
            t = time.time()
            step2_QC(conf_dict,logfile,filelist)
            step2time = time.time()-t
            wlog("running time for Step3: %s"%(step2time),logfile)

    if int(args.task) ==4:
        wlog("Step4: summary and output",logfile)
        final_file = conf_dict['General']['outputdirectory'] + 'summary/' + conf_dict['General']['outname'] + '_summary.pdf'
        if os.path.isfile(final_file):
            print 'Step4 has been finished.'
        else:
            t = time.time()
            step5_summary(conf_dict,logfile,filelist)
            step5time = time.time()-t
            wlog("running time for Step4: %s"%(step5time),logfile)

    
if __name__== '__main__':
    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt CAM\n")
        sys.exit(1)

