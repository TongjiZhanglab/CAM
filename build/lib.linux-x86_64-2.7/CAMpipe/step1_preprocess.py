#!/usr/bin/env python

# ------------------------------------
# Python Modual
# ------------------------------------

import os
import sys
import string
import subprocess
import time
# --------------------------
# custom package
# --------------------------

### tool function
from CAMpipe.Utility      import (sp,
                                   pdf_name,
                                   raise_error,
                                   wlog,
                                   ewlog,
                                   rwlog,
                                   flog,   
                                   CMD,
                                   createDIR,
                                   correct_genome_length,
                                   fastq_reads_length,
                                   transform_sam2bed,
                                   transform_bed2bed,
                                   bed2allbw,check_filelist
                                   )
# --------------------------
# main 
# --------------------------

def step1_preprocess(conf_dict,logfile,filelist):
    '''
    main data processing step, including mapping, transformation
    for fastq format : 
        bowtie mapping,
        q30 filter (optional),
        transform to bed 
    for sam format:
        q30 filter (optional),
        transform to bed    
    ''' 
    ### create annotation dir
    t0 = time.time()

    wlog('generate annotation files',logfile)
    genome_annotation_dir = conf_dict['General']['outputdirectory'] + 'annotation/'
    createDIR(genome_annotation_dir)
    os.chdir(genome_annotation_dir)
    
    SPNAME_raw = conf_dict['Step1_preprocess']['genome_fasta'].split('/')[-1].split('.')
    if len(SPNAME_raw) > 1:
        SPNAME = ".".join(SPNAME_raw[:-1])
    else:
        SPNAME = SPNAME_raw[0]
    if conf_dict['Step1_preprocess']['usefa'] == 1:

        genome2bit_gen_cmd = 'faToTwoBit %s %s_genome2bit.2bit'%(conf_dict['Step1_preprocess']['genome_fasta'],SPNAME)
        rwlog(genome2bit_gen_cmd,logfile)
        conf_dict['Step1_preprocess']['genome2bit'] = genome_annotation_dir + SPNAME + '_genome2bit.2bit'

    else:
        conf_dict['Step1_preprocess']['genome2bit'] = conf_dict['Step1_preprocess']['genome_fasta']
        fa_gen_cmd = 'twoBitToFa %s %s_genomefa.fa'%(conf_dict['Step1_preprocess']['genome2bit'],SPNAME)
        rwlog(fa_gen_cmd,logfile)
        conf_dict['Step1_preprocess']['genome_fasta'] = genome_annotation_dir + SPNAME + '_genomefa.fa'
            
    if conf_dict['General']['format'] == 'fastq':
        wlog('check bowtie index',logfile)
        if conf_dict['Step1_preprocess']['mapindex'] == "":
            wlog("no bowtie mapindex inputted",logfile)
            generate_index = 1
        elif not os.path.isfile( conf_dict['Step1_preprocess']['mapindex']+'.1.ebwt' ):
            wlog("cannot find inputted bowtie index file : %s "%(conf_dict['Step1_preprocess']['mapindex']+'.1.ebwt'),logfile)
            generate_index = 1
        else:
            generate_index = 0
        if generate_index == 1:
            wlog('generate bowtie index',logfile)
            if conf_dict['Step1_preprocess']['usefa'] == 1:
                genome2bit_gen_cmd = 'faToTwoBit %s %s_genome2bit.2bit'%(conf_dict['Step1_preprocess']['genome_fasta'],SPNAME)
                rwlog(genome2bit_gen_cmd,logfile)
                conf_dict['Step1_preprocess']['genome2bit'] = genome_annotation_dir + SPNAME + '_genome2bit.2bit'            
            else:
                conf_dict['Step1_preprocess']['genome2bit'] = conf_dict['Step1_preprocess']['genome_fasta']
                fa_gen_cmd = 'twoBitToFa %s %s_genomefa.fa'%(conf_dict['Step1_preprocess']['genome2bit'],SPNAME)
                rwlog(fa_gen_cmd,logfile)
                conf_dict['Step1_preprocess']['genome_fasta'] = genome_annotation_dir + SPNAME + '_genomefa.fa'
            mapindex_gen_cmd = 'bowtie-build %s %s_bowtie_index'%(conf_dict['Step1_preprocess']['genome_fasta'],SPNAME)
            rwlog(mapindex_gen_cmd,logfile)
            conf_dict['Step1_preprocess']['mapindex'] = genome_annotation_dir + SPNAME + '_bowtie_index'
            wlog('generate bowtie index done: %s'%(conf_dict['Step1_preprocess']['mapindex']),logfile)
        else:
            wlog('use user inputted bowtie index',logfile)
        
    wlog("prepare annotation done, time: %s"%(time.time()-t0),logfile)
    ### create mapping dir 
    mapping_dir = conf_dict['General']['outputdirectory'] + 'preprocess/'
    createDIR(mapping_dir)
    
    ### check reads file format , 
    ## start mapping step if format is fastq
    os.chdir(mapping_dir)

    if conf_dict['General']['format'] == 'fastq':
        ### sam file name
        conf_dict['Step1_preprocess']['sam'] = mapping_dir + conf_dict['General']['outname'] + '.sam'
        conf_dict['Step1_preprocess']['maplog'] = mapping_dir + conf_dict['General']['outname'] + '.bowtieout'
        # check the progress for mapping
        if check_filelist(filelist,conf_dict['Step1_preprocess']['sam'])==0 or not os.path.isfile(conf_dict['Step1_preprocess']['sam']):
            wlog('Now start mapping in %s , all mapping result will be here'%(mapping_dir),logfile)
            t1 = time.time()
            ### judge fastq reads length
            ReadLengthA = fastq_reads_length(conf_dict['General']['inputa'])
            if ReadLengthA[1] == 'difflen':
                wlog('WARNING: reads length in input fastq %s is not consistant'%(conf_dict['General']['inputa']),logfile)
            if conf_dict['General']['seqtype'] == 'PE':
                ReadLengthB = fastq_reads_length(conf_dict['General']['inputb'])
                if ReadLengthB[1] == 'difflen':
                    wlog('WARNING: reads length in input fastq (part2) %s is not consistant, nucpipe exit'%(conf_dict['General']['inputb']),logfile)
                if ReadLengthA[0] != ReadLengthB[0]:
                    wlog('WARNING: read length in 2 part of input fastq file is different',logfile)
                ReadLength = min(ReadLengthA[0] ,ReadLengthB[0])
            else:
                ReadLength = ReadLengthA[0]
            wlog('reads length is detected as %s'%(str(ReadLength)),logfile)
            
            ### check : read_length - trim3end > 18 
            if int(conf_dict['Step1_preprocess']['trim3end']) == 0:
                pass
            else:
                read_left_length = ReadLength - int(conf_dict['Step1_preprocess']['trim3end'])
                if read_left_length < 18:
                    ewlog('user set trim3end length is %s, the left reads length is %sbp, less than 18bp, CAM exit'%(conf_dict['Step1_preprocess']['trim3end'],str(read_left_length)),logfile)
            
            ### check bowtie software
            if sp('which bowtie')[0].strip() == "":
                ewlog('bowtie is not detected in default PATH, make sure you installed bowtie and export it into default PATH',logfile)
            # SE / PE bowtie mapping
            if conf_dict['General']['seqtype'] == 'PE':
                wlog('seqtype is PE (paired end), mapping with bowtie paired end mode',logfile)
                mapping_cmd = 'bowtie -X %s -3 %s --chunkmbs 256 -m 1 -p %s -S %s -1 %s -2 %s  %s 2>&1 >>/dev/null |tee -a %s'%( \
                               conf_dict['Step1_preprocess']['fragment_length_limit'],\
                               conf_dict['Step1_preprocess']['trim3end'], \
                               conf_dict['Step1_preprocess']['mapping_p'], \
                               conf_dict['Step1_preprocess']['mapindex'], \
                               conf_dict['General']['inputa'], \
                               conf_dict['General']['inputb'], \
                               conf_dict['Step1_preprocess']['sam'], \
                               conf_dict['Step1_preprocess']['maplog'])
                rwlog(mapping_cmd,logfile)
            elif conf_dict['General']['seqtype'] == 'SE':
                wlog('seqtype is SE (single end), mapping with bowtie single end mode',logfile)
                mapping_cmd = 'bowtie -3 %s --chunkmbs 256 -m 1 -p %s -S %s  %s  %s   >>/dev/null |tee -a %s'%( \
                               conf_dict['Step1_preprocess']['trim3end'], \
                               conf_dict['Step1_preprocess']['mapping_p'], \
                               conf_dict['Step1_preprocess']['mapindex'], \
                               conf_dict['General']['inputa'], \
                               conf_dict['Step1_preprocess']['sam'], \
                               conf_dict['Step1_preprocess']['maplog'])
                rwlog(mapping_cmd,logfile)
                # write the finished file into progress_filelist.txt
                flog(conf_dict['Step1_preprocess']['sam'],filelist)
            else:
                ewlog('wrong seqtype, current seqtype is %s'%(conf_dict['General']['seqtype']))
            wlog('mapping done, time: %s'%(time.time()-t1),logfile)
        else:
            wlog('mapping has been finished.',logfile)
    ### for sam/bed file, skip mapping step        
    elif conf_dict['General']['format'] == 'sam':
        wlog('input file format is sam, skip mapping step',logfile)
        conf_dict['Step1_preprocess']['sam'] = conf_dict['General']['inputa']
    elif conf_dict['General']['format'] == 'bed':
        wlog('input file format is bed, skip mapping step',logfile)
    else: 
        ewlog('input file is not in a proper format (fastq/sam/bed), current input file is %s. CAM exit'%(conf_dict['General']['inputa']),logfile)                

    ### transform sam to bed
    
    t2 = time.time()
    if conf_dict['General']['seqtype']=='PE':
        conf_dict['Step1_preprocess']['bed'] = mapping_dir + conf_dict['General']['outname'] + '_PEtoSE.bed'
    else:
        conf_dict['Step1_preprocess']['bed'] = mapping_dir + conf_dict['General']['outname'] + '.bed'
    if conf_dict['General']['format'] == 'fastq' or conf_dict['General']['format'] == 'sam':
        # check the progress for samtobed
        if check_filelist(filelist,conf_dict['Step1_preprocess']['bed'])==0 or not os.path.isfile(conf_dict['Step1_preprocess']['bed']):
            map_reads, fraglen_Dict = transform_sam2bed(conf_dict['Step1_preprocess']['sam'],\
                                                        conf_dict['Step1_preprocess']['bed'],\
                                                        conf_dict['General']['seqtype'],\
                                                        conf_dict['Step1_preprocess']['q30filter'],\
                                                        conf_dict['Step1_preprocess']['fragment_length_limit'])
            # write the finished file into progress_filelist.txt
            flog(conf_dict['Step1_preprocess']['bed'],filelist)

    else:
        ## for PEbed input, convert to SEbed, output only one strand for the fragments with only 1bp; for SEbed input, check format
        # check the progress for samtobed
        # finished_files = open(filelist).readlines()
        # finished_files = [i.strip() for i in finished_files]
        # status = 1
        # if conf_dict['Step1_preprocess']['bed'] not in finished_files:
        #     status = 0
        # print status
        # print check_filelist(filelist,conf_dict['Step1_preprocess']['bed'])==0
        # sys.exit(1)
        # print check_filelist(filelist,conf_dict['Step1_preprocess']['bed'])==0
        # sys.exit(1)
        if check_filelist(filelist,conf_dict['Step1_preprocess']['bed'])==0 or not os.path.isfile(conf_dict['Step1_preprocess']['bed']):
            map_reads, fraglen_Dict = transform_bed2bed(conf_dict['General']['inputa'],\
                                                        conf_dict['Step1_preprocess']['bed'],\
                                                        conf_dict['General']['seqtype'],\
                                                        conf_dict['Step1_preprocess']['q30filter'],\
                                                        conf_dict['Step1_preprocess']['fragment_length_limit'])
            flog(conf_dict['Step1_preprocess']['bed'],filelist)
            conf_dict['Step2_QC']['map_reads'] = map_reads
            conf_dict['Step2_QC']['freglen_dict'] = fraglen_Dict
            outf3 = open(mapping_dir+'freglen_tmp.txt','w')
            for i in sorted(conf_dict['Step2_QC']['freglen_dict'].keys()):
                newll = [i,conf_dict['Step2_QC']['freglen_dict'][i]]
                outf3.write("\t".join(map(str,newll))+"\n")
            outf3.close()
            flog(mapping_dir+'freglen_tmp.txt',filelist)
    # fraglen_Dict should be write into files

    transform_time = time.time() -t2
    wlog('transforming done, time: %s'%(transform_time),logfile)
    
    t3 = time.time()
    
    ### generate extbw, +1bpbw, -1bpbw, 10M +1bp bed, 10M -1bp bed
    conf_dict['Step1_preprocess']['profilebw'] = mapping_dir + conf_dict['General']['outname'] + '_profile.bw'
    conf_dict['Step1_preprocess']['profilebdg'] = mapping_dir + conf_dict['General']['outname'] + '_profile.bdg'
    conf_dict['Step1_preprocess']['centerwig'] = mapping_dir + conf_dict['General']['outname'] + '_center.wig'
    conf_dict['Step1_preprocess']['sdplus1bpbed'] = mapping_dir + conf_dict['General']['outname'] + '_sdplus1bp.bed'
    conf_dict['Step1_preprocess']['minus1bpbed'] = mapping_dir + conf_dict['General']['outname'] + '_minus1bp.bed'

    conf_dict['Step1_preprocess']['genome_length_use'] = mapping_dir + '%s_GemomeLengthTmp.genome'%(conf_dict['General']['outname'])
    correct_genome_length(conf_dict['Step1_preprocess']['genome_length'],conf_dict['Step1_preprocess']['genome_length_use'])
    # check the progress for bed2bw
    if check_filelist(filelist,conf_dict['Step1_preprocess']['centerwig'])==0 or not os.path.isfile(conf_dict['Step1_preprocess']['centerwig']):
        if not conf_dict['Step2_QC'].has_key('map_reads'):
            coverage_count_cmd = 'wc -l %s' % conf_dict['Step1_preprocess']['bed']
            # print coverage_count_cmd
            # print os.system(coverage_count_cmd)
            # print os.popen(coverage_count_cmd).read()
            conf_dict['Step2_QC']['map_reads'] = int(sp(coverage_count_cmd)[0].split()[0])
        print conf_dict['Step2_QC']['map_reads']            
        bed2allbw(conf_dict['Step1_preprocess']['bed'],\
                  conf_dict['Step1_preprocess']['genome_length_use'],\
                  conf_dict['Step1_preprocess']['profilebdg'],\
                  conf_dict['Step1_preprocess']['profilebw'],\
                  conf_dict['Step1_preprocess']['centerwig'],\
                  conf_dict['Step1_preprocess']['minus1bpbed'],\
                  conf_dict['Step1_preprocess']['sdplus1bpbed'],\
                  conf_dict['Step2_QC']['map_reads'],\
                  conf_dict['Step2_QC']['sample_down_reads'],\
                  conf_dict['Step1_preprocess']['rpm'])
        flog(conf_dict['Step1_preprocess']['centerwig'],filelist)
    else:
        wlog('pileup has been finished.',logfile)
    bed2bw_time = time.time() -t3
    wlog('pileup done, time: %s'%(bed2bw_time),logfile)
    # transform bed to bb
    # for fastq input, samtool view to transform sam to bed
    # than sort the bed file, and transform to bb
    if conf_dict['General']['seqtype']=='PE':
        conf_dict['Step1_preprocess']['mapped_bed'] = mapping_dir + conf_dict['General']['outname'] + '.bed'
        if check_filelist(filelist,conf_dict['Step1_preprocess']['mapped_bed'])==0  or not os.path.isfile(conf_dict['Step1_preprocess']['mapped_bed']):
            if conf_dict['General']['format']!='bed':
                cmd = 'samtools view %s -Sb | bamToBed -i stdin > %s' % (conf_dict['Step1_preprocess']['sam'],conf_dict['Step1_preprocess']['mapped_bed'])
                rwlog(cmd,logfile)
                flog(conf_dict['Step1_preprocess']['mapped_bed'],filelist)
            else:
                cmd = 'cp %s %s' % (conf_dict['General']['inputa'],conf_dict['Step1_preprocess']['mapped_bed'])
                rwlog(cmd,logfile)
                flog(conf_dict['Step1_preprocess']['mapped_bed'],filelist)
    else:
        conf_dict['Step1_preprocess']['mapped_bed'] = conf_dict['Step1_preprocess']['bed']
    
    conf_dict['Step1_preprocess']['mapped_bb'] = mapping_dir + conf_dict['General']['outname'] + '.bb'

    # cmd = 'awk \'{FS="\\t";OFS="\\t";if (NF==6 && $2>=0 && $3>=0) print $0;}\' %s | sort -k1,1 -k2,2n > %s.tmp' % (conf_dict['Step1_preprocess']['mapped_bed'],conf_dict['Step1_preprocess']['mapped_bed'])
    # try :
    if check_filelist(filelist,'%s.tmp' % conf_dict['Step1_preprocess']['mapped_bed'])==0 or not os.path.isfile(conf_dict['Step1_preprocess']['mapped_bed']):
        cmd = 'awk %s{FS="\\t";OFS="\\t";if (NF==6 && $2>=0 && $3>=0) print $0;}%s %s | sort -k1,1 -k2,2n > %s.tmp' % ("'","'",conf_dict['Step1_preprocess']['mapped_bed'],conf_dict['Step1_preprocess']['mapped_bed'])
        rwlog(cmd,logfile)
        flog('%s.tmp' % conf_dict['Step1_preprocess']['mapped_bed'],filelist)
    if sp('which bedToBigBed')[0].strip() == "":
        ewlog('bedToBigBed is not detected in default PATH, make sure you installed bowtie and export it into default PATH',logfile)
    else:
        if check_filelist(filelist,conf_dict['Step1_preprocess']['mapped_bb'])==0 or not os.path.isfile(conf_dict['Step1_preprocess']['mapped_bb']):
            cmd = 'bedToBigBed %s.tmp %s %s' % (conf_dict['Step1_preprocess']['mapped_bed'],conf_dict['Step1_preprocess']['genome_length_use'],conf_dict['Step1_preprocess']['mapped_bb'])
            rwlog(cmd,logfile)
            flog(conf_dict['Step1_preprocess']['mapped_bb'],filelist)
    # cmd = 'rm %s.tmp ' % conf_dict['Step1_preprocess']['mapped_bed']
    # rwlog(cmd,logfile)
    # except:
        # print 'Failed to transform bed to bb.'

    return conf_dict
