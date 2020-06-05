#!/usr/bin/env python

# ------------------------------------
# Python Modual
# ------------------------------------

import os
import sys
import string
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
                                   calculate_genome_length,
                                   tssprofile,
                                   ATprofile,
                                   reads_reads_distance,
                                   generate_utr_overlap,
                                   generate_DHS_overlap,
                                   customprofile,
                                   check_filelist
                                   )
# --------------------------
# main 
# --------------------------
def step2_QC(conf_dict,logfile,filelist):
    '''
    QC step
    mapping stat
    4 QC plots
    '''
    # start    
    ### create  QC dir and conduct QC
    qcdir = conf_dict['General']['outputdirectory'] + 'QC/'
    mapping_dir = conf_dict['General']['outputdirectory'] + 'preprocess/'
    arraydir = conf_dict['General']['outputdirectory'] + 'nucarray/'
    createDIR(qcdir)
    os.chdir(qcdir)

    ### start QC plots

    ### TSS profile

    t1 = time.time()
    # resstore dictionary information that may skipped due to process interruption and resart
    conf_dict['Step2_QC']['ave_tssprofile'] = qcdir + conf_dict['General']['outname'] + '_Tss_profile.txt'
    conf_dict['Step1_preprocess']['profilebw'] = mapping_dir + conf_dict['General']['outname'] + '_profile.bw'
    print conf_dict['Step2_QC']['ave_tssprofile']
    if check_filelist(filelist,conf_dict['Step2_QC']['ave_tssprofile'])==0 or not os.path.isfile(conf_dict['Step2_QC']['ave_tssprofile']):
        tssprofile(conf_dict['Step1_preprocess']['profilebw'],\
               conf_dict['Step2_QC']['ave_tssprofile'],\
               conf_dict['Step1_preprocess']['gene_annotation'],\
               conf_dict['Step2_QC']['upstreamtss'],\
               conf_dict['Step2_QC']['downstreamtss'])
        flog(conf_dict['Step2_QC']['ave_tssprofile'],filelist)
        wlog('tss profile done, time: %s'%(time.time()-t1),logfile)
    else:
      wlog('tss profile has been finished.',logfile)
    
    ### custom region profile
    if int(conf_dict['Step2_QC']['plotcustom']) == 1:
        t2 = time.time()
        # resstore dictionary information that may skipped due to process interruption and resart
        conf_dict['Step2_QC']['custom_profile'] = qcdir + conf_dict['General']['outname'] + '_CustomRegion_profile.txt'
        if check_filelist(filelist,conf_dict['Step2_QC']['custom_profile'])==0 or not os.path.isfile(conf_dict['Step2_QC']['custom_profile']):
            customprofile(conf_dict['Step1_preprocess']['profilebw'],\
                      conf_dict['Step2_QC']['custom_profile'],\
                      conf_dict['General']['customregion'],\
                      conf_dict['Step2_QC']['customregion_dis'])
            wlog('custom region profile done, time: %s'%(time.time()-t2),logfile)
            flog(conf_dict['Step2_QC']['custom_profile'],filelist)
        else:
            wlog('custom region profile has been finished.')
            
    ### AT frac
    t3 = time.time()
    # resstore dictionary information that may skipped due to process interruption and resart
    conf_dict['Step2_QC']['at_frac'] = qcdir + conf_dict['General']['outname'] + '_ATfrac.txt'
    conf_dict['Step1_preprocess']['sdplus1bpbed'] = mapping_dir + conf_dict['General']['outname'] + '_sdplus1bp.bed' 
    conf_dict['Step1_preprocess']['minus1bpbed'] = mapping_dir + conf_dict['General']['outname'] + '_minus1bp.bed'
    conf_dict['Step1_preprocess']['genome_length_use'] = mapping_dir + '%s_GemomeLengthTmp.genome'%(conf_dict['General']['outname'])
    genome_annotation_dir = conf_dict['General']['outputdirectory'] + 'annotation/'
    SPNAME_raw = conf_dict['Step1_preprocess']['genome_fasta'].split('/')[-1].split('.')
    if len(SPNAME_raw) > 1:
        SPNAME = ".".join(SPNAME_raw[:-1])
    else:
        SPNAME = SPNAME_raw[0]
    if conf_dict['Step1_preprocess']['usefa'] == 1:
        conf_dict['Step1_preprocess']['genome2bit'] = genome_annotation_dir + SPNAME + '_genome2bit.2bit'
    else:
        conf_dict['Step1_preprocess']['genome2bit'] = conf_dict['Step1_preprocess']['genome_fasta']

    if check_filelist(filelist,conf_dict['Step2_QC']['at_frac'])==0  or not os.path.isfile(conf_dict['Step2_QC']['at_frac']):
        ATfrac = ATprofile(conf_dict['Step1_preprocess']['sdplus1bpbed'],conf_dict['Step1_preprocess']['genome2bit'])
        outf2 = open(conf_dict['Step2_QC']['at_frac'],'w')
        for i in range(len(ATfrac)):
            newll = [i+4, ATfrac[i]]        
            outf2.write("\t".join(map(str,newll))+"\n")
        outf2.close()
        wlog('AT profile done, time: %s'%(time.time()-t3),logfile)
        flog(conf_dict['Step2_QC']['at_frac'],filelist)
    else:
        wlog('AT profile has finished.',logfile)
    
    ### RR distance /frag len
    t4 = time.time()
    conf_dict['Step2_QC']['fraglen'] = qcdir + conf_dict['General']['outname'] + '_fraglen.txt'

    if conf_dict['General']['seqtype'] == "SE":
        if check_filelist(filelist,conf_dict['Step2_QC']['fraglen'])==0  or not os.path.isfile(conf_dict['Step2_QC']['fraglen']):
            nucleosome_length = reads_reads_distance(conf_dict['Step1_preprocess']['sdplus1bpbed'],\
                                             conf_dict['Step1_preprocess']['minus1bpbed'],\
                                             conf_dict['Step1_preprocess']['genome_length_use'],\
                                             conf_dict['Step2_QC']['reads_distance_range'])
            outf3 = open(conf_dict['Step2_QC']['fraglen'],'w')    
            for i in range(len(nucleosome_length)):
                newll = [i,nucleosome_length[i]]
                outf3.write("\t".join(map(str,newll))+"\n")
            outf3.close()
            flog(conf_dict['Step2_QC']['fraglen'],filelist)
            wlog('nucleosomal DNA length done, time: %s'%(time.time()-t4),logfile)
    else:
        if check_filelist(filelist,conf_dict['Step2_QC']['fraglen'])==0  or not os.path.isfile(conf_dict['Step2_QC']['fraglen']):
            rwlog('cp %s %s'%(mapping_dir + 'freglen_tmp.txt',qcdir+conf_dict['General']['outname'] + '_fraglen.txt' ),logfile)
            flog(conf_dict['Step2_QC']['fraglen'],filelist)
        # storege in conf_dic
        # outf3 = open(conf_dict['Step2_QC']['fraglen'],'w')
        # for i in sorted(conf_dict['Step2_QC']['freglen_dict'].keys()):
        #     newll = [i,conf_dict['Step2_QC']['freglen_dict'][i]]
        #     outf3.write("\t".join(map(str,newll))+"\n")
        # outf3.close()    
            wlog('nucleosomal DNA length done, time: %s'%(time.time()-t4),logfile)
    


    ###  feature
    t5 = time.time()

    if conf_dict['Step1_preprocess']['species'] in ['hg38','hg19']:
        effective_gs = 2.7e9
    elif conf_dict['Step1_preprocess']['species'] in ['mm10','mm9']:
        effective_gs = 1.87e9
    else:
        effective_gs = int( calculate_genome_length(conf_dict['Step1_preprocess']['genome_length_use']) )

    # seq coverage
    try:
        seq_coverage = float(conf_dict['Step2_QC']['map_reads']) * 194 / effective_gs
    except:
        coverage_count_cmd = 'wc -l %s' % mapping_dir + conf_dict['General']['outname']+'.bed'
        seq_coverage = float(int(sp(coverage_count_cmd)[0].split()[0])) * 194 / effective_gs
        # seq_coverage = 1000000
        # count bed file
    # array number
    conf_dict['Step3_nucarray']['arrayselect'] = arraydir + conf_dict['General']['outname'] + '_Nucleosome_Array.bed'
    array_count_cmd = 'wc -l %s'%(conf_dict['Step3_nucarray']['arrayselect'])
    # print array_count_cmd
    array_num = int(sp(array_count_cmd)[0].split()[0])

    # array on UTR 3kb
    array_on_utr, total_utr_length = generate_utr_overlap(conf_dict['Step3_nucarray']['arrayselect'], conf_dict['Step1_preprocess']['gene_annotation'], 3000)
    array_on_DHS, total_DHS_length = generate_DHS_overlap(conf_dict['Step3_nucarray']['arrayselect'], conf_dict['Step1_preprocess']['union_dhs'])    
    utr_fold = round( (float(array_on_utr)/array_num) / (float(total_utr_length)/effective_gs) ,4)
    DHS_fold = round( (float(array_on_DHS)/array_num) / (float(total_DHS_length)/effective_gs) ,4)
    # gene level array annotation +-3kb
    # print conf_dict['rscript']
    ### generate QC plots
    cmd = 'Rscript %s %s %s %s %s %s %s %s %s'%(conf_dict['rscript'] + 'QCplots.r',\
                                    conf_dict['General']['outname'],\
                                    conf_dict['Step2_QC']['upstreamtss'],\
                                    conf_dict['Step2_QC']['downstreamtss'],\
                                    conf_dict['Step2_QC']['plotcustom'],\
                                    seq_coverage,\
                                    conf_dict['Step2_QC']['reads_distance_range'],\
                                    utr_fold,\
                                    DHS_fold)
    wlog(cmd,logfile)
    rot_score, nuclen, NFRscore, PSarray = sp(cmd)[0].split()[-4:]
    
    outf4 = open(qcdir + conf_dict['General']['outname'] + '_Features.txt','w')

    outf4.write("\t".join(map(str,['seq_coverage', seq_coverage]))+"\n")
    outf4.write("\t".join(map(str,['rot_score', rot_score]))+"\n")
    outf4.write("\t".join(map(str,['nuclen', nuclen]))+"\n")
    outf4.write("\t".join(map(str,['NFRscore', NFRscore]))+"\n")
    outf4.write("\t".join(map(str,['PSarray', PSarray]))+"\n")
    outf4.write("\t".join(map(str,['array_on_utr', array_on_utr]))+"\n")
    outf4.write("\t".join(map(str,['array_on_DHS', array_on_DHS]))+"\n")
    outf4.write("\t".join(map(str,['array_num', array_num]))+"\n")
    outf4.write("\t".join(map(str,['total_utr_length', total_utr_length]))+"\n")
    outf4.write("\t".join(map(str,['total_DHS_length', total_DHS_length]))+"\n")
    outf4.write("\t".join(map(str,['effective_gs', effective_gs]))+"\n")
    outf4.write("\t".join(map(str,['enrichment_on_UTR', utr_fold]))+"\n")
    outf4.write("\t".join(map(str,['enrichment_on_DHS', DHS_fold]))+"\n")
    outf4.close()
    flog(qcdir + conf_dict['General']['outname'] + '_Features.txt',filelist)

    wlog("QC plots and feature summary done, time: %s"%(time.time()-t5),logfile)
    return conf_dict









