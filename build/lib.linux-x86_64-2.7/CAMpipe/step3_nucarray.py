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
                                   modify_wig_signal,
                                   generate_position_signal,
                                   make_array,
                                   generate_geneLevel_arrayAnnotation,
                                   signal_on_aray,check_filelist)
# --------------------------
# main 
# --------------------------
def step3_nucarray(conf_dict,logfile,filelist):
    '''
    analysis part
    mainly Rscript
    dimentional reduction + clustering
    '''   

    arraydir = conf_dict['General']['outputdirectory'] + 'nucarray/'
    mapping_dir = conf_dict['General']['outputdirectory'] + 'preprocess/'
    createDIR(arraydir)
    os.chdir(arraydir)
 
    win_size = 7
    band_width = 3
    
    conf_dict['Step3_nucarray']['modifywig_center'] = arraydir + conf_dict['General']['outname'] + '_center_modify.wig'
    conf_dict['Step3_nucarray']['positionwig_center'] = arraydir + conf_dict['General']['outname'] + '_center_position.wig'
    conf_dict['Step3_nucarray']['positionbw_center'] = arraydir + conf_dict['General']['outname'] + '_center_position.bw'
    conf_dict['Step3_nucarray']['arrayall'] = arraydir + conf_dict['General']['outname'] + '_Nucleosome_Array_all.bed'
    conf_dict['Step3_nucarray']['arrayselect'] = arraydir + conf_dict['General']['outname'] + '_Nucleosome_Array.bed'
    conf_dict['Step3_nucarray']['gene_array_anno'] = arraydir + conf_dict['General']['outname'] + '_geneLevel_nucarrayAnnotation.bed'
    conf_dict['Step3_nucarray']['profilebw_onarray'] = arraydir + conf_dict['General']['outname'] + "_profile_on_Nucleosome_Array.bw"
    
    if check_filelist(filelist,conf_dict['Step3_nucarray']['modifywig_center'])==0  or not os.path.isfile(conf_dict['Step3_nucarray']['modifywig_center']):
        conf_dict['Step1_preprocess']['centerwig'] = mapping_dir + conf_dict['General']['outname'] + '_center.wig'
        modify_wig_signal(conf_dict['Step1_preprocess']['centerwig'],\
                      conf_dict['Step3_nucarray']['modifywig_center'],\
                      conf_dict['Step3_nucarray']['window_size'],\
                      conf_dict['Step3_nucarray']['smooth_bandwidth'])
        flog(conf_dict['Step3_nucarray']['modifywig_center'],filelist)
    else:
        wlog('modifywig_center wig exists.',logfile)

    bg_center = generate_position_signal(conf_dict['Step3_nucarray']['modifywig_center'],\
                             conf_dict['Step3_nucarray']['positionwig_center'],\
                             conf_dict['Step3_nucarray']['window_size'])
    conf_dict['Step3_nucarray']['bg_value'] = bg_center
    if check_filelist(filelist,conf_dict['Step3_nucarray']['positionwig_center'])==0 or not os.path.isfile(conf_dict['Step3_nucarray']['positionwig_center']):
        make_array(conf_dict['Step3_nucarray']['positionwig_center'],\
               conf_dict['Step3_nucarray']['arrayall'],\
               bg_center)
        flog(conf_dict['Step3_nucarray']['arrayall'],filelist)
        array_cmd = """awk '{if ($3 - $2 > %s && $5 > %s) print $0;}' %s > %s """%(conf_dict['Step3_nucarray']['array_length'],\
                                                                         conf_dict['Step3_nucarray']['array_fold'],\
                                                                         conf_dict['Step3_nucarray']['arrayall'],\
                                                                         conf_dict['Step3_nucarray']['arrayselect'])
        rwlog(array_cmd,logfile)
        flog(conf_dict['Step3_nucarray']['arrayselect'],filelist)
        wlog('Nucleosome array background value:%s'%bg_center,logfile)
    if check_filelist(filelist,conf_dict['Step3_nucarray']['gene_array_anno'])==0 or not os.path.isfile(conf_dict['Step3_nucarray']['gene_array_anno']):
        generate_geneLevel_arrayAnnotation(conf_dict['Step3_nucarray']['arrayselect'],conf_dict['Step1_preprocess']['gene_annotation'],conf_dict['Step3_nucarray']['gene_array_anno'])
        flog(conf_dict['Step3_nucarray']['gene_array_anno'],filelist)
    if check_filelist(filelist,conf_dict['Step3_nucarray']['profilebw_onarray'])==0 or not os.path.isfile(conf_dict['Step3_nucarray']['profilebw_onarray']):
        conf_dict['Step1_preprocess']['genome_length_use'] = mapping_dir + '%s_GemomeLengthTmp.genome'%(conf_dict['General']['outname'])
        conf_dict['Step1_preprocess']['profilebdg'] = mapping_dir + conf_dict['General']['outname'] + '_profile.bdg'
        signal_on_aray( conf_dict['Step1_preprocess']['profilebdg'],\
                    conf_dict['Step3_nucarray']['arrayselect'], \
                    conf_dict['Step3_nucarray']['profilebw_onarray'],\
                    conf_dict['Step1_preprocess']['genome_length_use'])
        flog(conf_dict['Step3_nucarray']['profilebw_onarray'],filelist)
    if check_filelist(filelist,conf_dict['Step3_nucarray']['positionbw_center'])==0 or not os.path.isfile(conf_dict['Step3_nucarray']['positionbw_center']):
        conf_dict['Step1_preprocess']['genome_length_use'] = mapping_dir + '%s_GemomeLengthTmp.genome'%(conf_dict['General']['outname'])
        cmd = 'wigToBigWig %s %s %s' % (conf_dict['Step3_nucarray']['positionwig_center'],conf_dict['Step1_preprocess']['genome_length_use'],conf_dict['Step3_nucarray']['positionbw_center'])
        rwlog(cmd,logfile)
        flog(conf_dict['Step3_nucarray']['positionbw_center'],filelist)
    return conf_dict


