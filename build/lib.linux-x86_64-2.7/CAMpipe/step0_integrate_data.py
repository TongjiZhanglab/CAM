#!/usr/bin/env python

# ------------------------------------
# Python Modual
# ------------------------------------

import os
import sys
from string import *

# --------------------------
# custom package
# --------------------------


### tool function
from CAMpipe.Utility      import (sp,
                                   sperr,
                                   pdf_name,
                                   raise_error,
                                   wlog,
                                   ewlog,
                                   CMD,
                                   checkfa,
                                   checkbed
                                   )
# --------------------------
# main 
# --------------------------

def step0_integrate_data(conf_dict,logfile):
    '''
    step0 integrate data 
    check and complement parameter
    '''
    wlog("Start CAM",logfile)
    wlog("Step0: Data integrate",logfile)
    
    ### check output name
    if "/" in conf_dict['General']['outname']:
        ewlog("outname is the name of all your output result, cannot contain "/", current outname is  %s"%(conf_dict['General']['outname']),logfile)

    ### check data path of inputa
    if "~" in conf_dict['General']['inputa']:
        ewlog('require absolute path for input file, input file cannot contain "~", current input file is %s'%(conf_dict['General']['inputa']),logfile)
    if not conf_dict['General']['inputa'].startswith('/'):
        conf_dict['General']['inputa'] = conf_dict['General']['startdir'] + conf_dict['General']['inputa']
    if not os.path.isfile(conf_dict['General']['inputa']):
        ewlog("input file %s not found"%(conf_dict['General']['inputa']),logfile)

    ### check data format
    if conf_dict['General']['seqtype'] == "PE":
        wlog('sequencing type is set as PE(paired end)',logfile)
    elif conf_dict['General']['seqtype'] == "SE":
        wlog('sequencing type is set as SE(single end)',logfile)
    else:
        ewlog('seqtype can only be SE/PE, current value is %s'%(conf_dict['General']['seqtype']),logfile)

    if conf_dict['General']['inputa'].endswith('.fastq') :
        conf_dict['General']['format'] = 'fastq' 
        wlog('input file is raw sequencing file (fastq)',logfile)
        if conf_dict['General']['seqtype'] == "PE":
            if conf_dict['General']['inputb'].strip() == "" :
                ewlog('2nd part of input file (inputb, fastq) is left blank while seqtype is PE (paired end). Please make sure you correctly input your data and specify the sequencing type. CAM exit',logfile)
            elif not conf_dict['General']['inputb'].endswith('.fastq'):
                ewlog('2nd part of input file (inputb) does not endwith .fastq. Please make sure you correctly input your data. CAM exit',logfile)
            ### check data path of inputb
            if "~" in conf_dict['General']['inputb']:
                ewlog('require absolute path for input file(part2), input file cannot contain "~", current input file(part2) is %s. CAM exit'%(conf_dict['General']['inputb']),logfile)
            if not conf_dict['General']['inputb'].startswith('/'):
                conf_dict['General']['inputb'] = conf_dict['General']['startdir'] + conf_dict['General']['inputb']
            if not os.path.isfile(conf_dict['General']['inputb']):
                ewlog("input file(part2) %s not found. CAM exit"%(conf_dict['General']['inputb']),logfile)
        elif conf_dict['General']['seqtype'] == "SE":
            if conf_dict['General']['inputb'].strip() != "" :
                wlog('2nd part of input file(inputb) is not left blank while seqtype is SE (single end), ignore inputb parameter',logfile)
                
    elif conf_dict['General']['inputa'].endswith('.sam'):
        conf_dict['General']['format'] = 'sam' 
        wlog('input file is aligned sequencing file (sam)',logfile)
        if conf_dict['General']['inputb'].strip() != "" :
            wlog('2nd part of input file(inputb) is not left blank while input file is aligned sam file, ignore inputb parameter',logfile)
    elif conf_dict['General']['inputa'].endswith('.bed'):
        conf_dict['General']['format'] = 'bed' 
        wlog('input file is aligned sequencing file (bed)',logfile)
        if conf_dict['General']['inputb'].strip() != "" :
            wlog('2nd part of input file(inputb) is not left blank while input file is aligned bed file, ignore inputb parameter',logfile)
    else: 
        ewlog('input file is not in a proper format (fastq/sam/bed), current input file is %s. CAM exit'%(conf_dict['General']['inputa']),logfile)                

    ### check custom region
    conf_dict['Step2_QC']['plotcustom'] = 1
    if conf_dict['General']['customregion'].strip() == "":
        wlog('no custom region input, custom region profile will be skipped',logfile)
        conf_dict['Step2_QC']['plotcustom'] = 0
    if "~" in conf_dict['General']['customregion']:
        wlog('require absolute path for custom region, custom region cannot contain "~", current input file is %s, custom region profile will be skipped'%(conf_dict['General']['customregion']),logfile)
        conf_dict['Step2_QC']['plotcustom'] = 0
    if not conf_dict['General']['customregion'].startswith('/'):
        conf_dict['General']['customregion'] = conf_dict['General']['startdir'] + conf_dict['General']['customregion']
    if not os.path.isfile(conf_dict['General']['customregion']):
        wlog("custom region %s not found, custom region profile will be skipped"%(conf_dict['General']['customregion']),logfile)
        conf_dict['Step2_QC']['plotcustom'] = 0
    if not conf_dict['General']['customregion'].endswith('.bed') or checkbed(conf_dict['General']['customregion']) == 0:
        wlog('custom region is not in bed format or not endswith .bed, current name of custom region is %s, custom region profile will be skipped'%(conf_dict['General']['customregion']),logfile)
        conf_dict['Step2_QC']['plotcustom'] = 0
                
    if conf_dict['Step2_QC']['plotcustom'] == 1  :
        if int(conf_dict['Step2_QC']['customregion_dis']) > 5000 or int(conf_dict['Step2_QC']['customregion_dis']) < 200:
            wlog("distance of custom region length (customregion_dis parameter) should be set between 200 ~ 5000, current value is %s, set customregion_dis to default 1000bp"%(conf_dict['Step2_QC']['customregion_dis'],logfile))
            conf_dict['Step2_QC']['customregion_dis'] = 1000
    
    ### check species parameter
    if conf_dict['Step1_preprocess']['species'] == "":
        ewlog("species is not given, CAM exit",logfile) 
    elif not conf_dict['Step1_preprocess']['species'] in ['hg38','hg19','mm10','mm9']:
        ewlog("species should be chose from [hg38,hg19,mm10,mm9](case sensitive), current species is %s, CAM exit"%(conf_dict['Step1_preprocess']['species']),logfile)         
    else:
        if not os.path.isfile(conf_dict['default_anno_dir'] + conf_dict['Step1_preprocess']['species'] + '_refgenes.txt') :
            ewlog("gene annotation file for genome version: %s (file should be %s_refgenes.txt) is not detected, make sure you already have (don't remove) corresponded genome version installed, CAM exit"%(conf_dict['Step1_preprocess']['species'], conf_dict['Step1_preprocess']['species']))
        if not os.path.isfile(conf_dict['default_anno_dir'] + conf_dict['Step1_preprocess']['species'] + '.genome'):
            ewlog("genome length file for genome version: %s (file should be %s.genome) is not detected, make sure you already have corresponded genome version installed, CAM exit"%(conf_dict['Step1_preprocess']['species'], conf_dict['Step1_preprocess']['species']))
        if not os.path.isfile(conf_dict['default_anno_dir'] + 'DHS_' + conf_dict['Step1_preprocess']['species'] + '.bed' ):
            ewlog("union DHS file for genome version: %s (file should be DHS_%s.bed) is not detected, make sure you already have corresponded genome version installed, CAM exit"%(conf_dict['Step1_preprocess']['species'], conf_dict['Step1_preprocess']['species']))
        else:
            wlog("corresponded annotation files for %s it detected"%(conf_dict['Step1_preprocess']['species']),logfile)     
            
    ### use default annotation file, if species is correctly inputted
    conf_dict['Step1_preprocess']['gene_annotation'] = conf_dict['default_anno_dir'] + conf_dict['Step1_preprocess']['species'] + '_refgenes.txt'
    conf_dict['Step1_preprocess']['genome_length'] = conf_dict['default_anno_dir'] + conf_dict['Step1_preprocess']['species'] + '.genome'
    conf_dict['Step1_preprocess']['union_dhs'] = conf_dict['default_anno_dir'] + 'DHS_' + conf_dict['Step1_preprocess']['species'] + '.bed'
                
    ### mapping index and 2bit
    if conf_dict['Step1_preprocess']['genome_fasta'] == "":
        ewlog("genome_fasta is not given, CAM exit",logfile) 
    else:
        if "~" in conf_dict['Step1_preprocess']['genome_fasta']:
            ewlog("genome_fasta file: %s is not absolute path, CAM exit"%(conf_dict['Step1_preprocess']['genome_fasta']),logfile)        
        if not conf_dict['Step1_preprocess']['genome_fasta'].startswith("/") : 
            conf_dict['Step1_preprocess']['genome_fasta'] = conf_dict['General']['startdir'] + conf_dict['Step1_preprocess']['genome_fasta']
        if os.path.isfile(conf_dict['Step1_preprocess']['genome_fasta']) :
            if conf_dict['Step1_preprocess']['genome_fasta'].endswith('.fa') and checkfa(conf_dict['Step1_preprocess']['genome_fasta']) == 1:
                wlog("genome_fasta file: %s is detected (fasta format)"%(conf_dict['Step1_preprocess']['genome_fasta']),logfile)
                conf_dict['Step1_preprocess']['usefa'] = 1
            elif conf_dict['Step1_preprocess']['genome_fasta'].endswith('.2bit') :
                wlog("genome_fasta file: %s is detected (2bit format)"%(conf_dict['Step1_preprocess']['genome_fasta']),logfile)
                conf_dict['Step1_preprocess']['usefa'] = 0
            else:
                ewlog("genome_fasta file: %s is not in fasta(.fa) or 2bit(.2bit) format, CAM exit"%(conf_dict['Step1_preprocess']['genome_fasta']),logfile)

        else:
            ewlog("genome_fasta file: %s is not a regular file, CAM exit"%(conf_dict['Step1_preprocess']['genome_fasta']),logfile)
            
 
    ### check options
    try:
        int(conf_dict['Step1_preprocess']['mapping_p'])
        #wlog('mapping thread is %s'%(str(int(conf_dict['Step1_preprocess']['mapping_p']))),logfile)
    except:
        ewlog('mapping_p should be int, current value is %s'%(conf_dict['Step1_preprocess']['mapping_p']),logfile)
    if int(conf_dict['Step1_preprocess']['trim3end']) < 0:
        ewlog('trim3end should be greater/equal than 0, current value is %s'%(conf_dict['Step1_preprocess']['trim3end']),logfile)
    
    if not int(conf_dict['Step1_preprocess']['q30filter']) in [0,1]:
        ewlog('q30filter measurement can only be 0/1, current value is %s'%(conf_dict['Step1_preprocess']['q30filter']),logfile)
    if not int(conf_dict['Step1_preprocess']['trim3end']) in [0,1]:
        ewlog('q30filter measurement can only be 0/1, current value is %s'%(conf_dict['Step1_preprocess']['q30filter']),logfile)
    if not int(conf_dict['Step1_preprocess']['rpm']) in [0,1]:
        ewlog('rpm function can only be 0/1, current value is %s'%(conf_dict['Step1_preprocess']['rpm']),logfile)
    
    if not int(conf_dict['Step2_QC']['sample_down_reads']) > 1000000:
        ewlog('sample_down_reads should greater than 1000000 to make sure enough reads for QC step, current value is %s'%(conf_dict['Step2_QC']['sample_down_reads']),logfile)

    if	not int(conf_dict['Step2_QC']['upstreamtss']) > 500:
        wlog('upstreamtss should be greater than 500 , current value is %s, CAM adjust upstreamtss to 500bp'%(conf_dict['Step2_QC']['upstreamtss']),logfile)
        conf_dict['Step2_QC']['upstreamtss'] = 500
    
    if	not int(conf_dict['Step2_QC']['downstreamtss']) > 500:
        wlog('downstreamtss should be greater than 500 , current value is %s, CAM adjust downstreamtss to 500bp'%(conf_dict['Step2_QC']['downstreamtss']),logfile)
        conf_dict['Step2_QC']['downstreamtss'] = 500

    if not int(conf_dict['Step2_QC']['reads_distance_range']) > 150:
        wlog('reads_distance_range should be greater than 150 , current value is %s, CAM adjust reads_distance_range to 150bp'%(conf_dict['Step2_QC']['reads_distance_range']),logfile)
        conf_dict['Step2_QC']['reads_distance_range'] = 150
    wlog('check options: DONE ',logfile)
    

    ### check Rscript
    #if not 'Usage' in sperr('Rscript')[1] and not 'version' in sperr('Rscript')[1]:
    #    ewlog('require Rscript',logfile)
    
    ### check pdflatex
    if sp('pdflatex --help')[0] == "":
        wlog('pdflatex was not installed, CAM is still processing but no summary QC report generated',logfile)
        conf_dict['General']['latex'] = 0
    else:
        conf_dict['General']['latex'] = 1

    wlog('Step0 Data integrate DONE',logfile)

    return conf_dict
    
    
    
