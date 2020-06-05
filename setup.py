#!/usr/bin/env python
"""Description
Setup script for CAM  -- QC and analysis pipeline for MNase-seq data
Copyright (c) 2017 Shengen Hu <tarelahu@gmail.com>
This code is free software; you can redistribute it and/or modify it
under the terms of the Artistic License (see the file COPYING included
with the distribution).
"""
import os
import sys
import subprocess
from distutils.core import setup, Extension



def sp(cmd):
    '''
    Call shell cmd or software and return its stdout
    '''
    a=subprocess.Popen(cmd, stdout=subprocess.PIPE, shell='TRUE')
    ac = a.communicate()
    return ac

def compile_bedtools():
    curdir = os.getcwd()
    os.chdir('refpackage/bedtools')
    sp('make 1>/dev/null 2>&1 ')
    sp('chmod 755 *')
    os.chdir(curdir)
    
def check_bedtools():
    checkhandle = sp('which bedtools')
    if checkhandle[0].strip() == "":
        return 0
    else:
        return 1
def check_R():
    checkhandle = sp('which Rscript')
    if checkhandle[0].strip() == "":
        return 0
    else:
        return 1       
def main(): 
    if sys.version_info[0] != 2 or sys.version_info[1] < 7:
	    print >> sys.stderr, "ERROR: CAM requires Python 2.7"
	    sys.exit()
    has_R = check_R()
    if has_R == 0:
	    print >> sys.stderr, "ERROR: CAM requires R & Rscript under default PATH"
	    sys.exit()
        
    has_bedtools = check_bedtools()
    print 'Intalling CAM, may take serval minutes'
    
    if has_bedtools == 0:
        print 'bedtools is not detected under default PATH, bedtools is also installed'
        compile_bedtools()
        SCRIPT = ['bin/CAM.py','refpackage/bedtools/bin/bedtools']
    else:
        SCRIPT = ['bin/CAM.py']
        
    for Tools in ['bedClip','bedGraphToBigWig','bigWigSummary','faToTwoBit','twoBitToFa']:
        if sp('which %s'%Tools)[0].strip() == "":
            sp('chmod 755 refpackage/%s'%Tools)
            SCRIPT.append('refpackage/%s'%Tools)
    
    if sp('which bowtie')[0].strip() == "":
        sp('chmod 755 refpackage/bowtie*')
        SCRIPT.extend(['refpackage/bowtie-align-l','refpackage/bowtie','refpackage/bowtie-align-s','refpackage/bowtie-align-l-debug','refpackage/bowtie-build-l-debug','refpackage/bowtie-build-l','refpackage/bowtie-build','refpackage/bowtie-align-s-debug','refpackage/bowtie-inspect-l','refpackage/bowtie-inspect','refpackage/bowtie-build-s-debug','refpackage/bowtie-build-s','refpackage/bowtie-inspect-s-debug','refpackage/bowtie-inspect-s'])

    PKGdir = {'CAMpipe' : 'lib'}
    PKGs = ['CAMpipe']
    
    try:
        import twobitreader
    except:
        print 'twobitreader is not detected, install twobitreader'
        PKGdir['twobitreader'] = 'twobit_src'
        PKGs.append('twobitreader')
        
    setup(name="CAMpipe",
          version="1.0.0",
          description="CAM: QC and analysis pipelie for MNase-seq",
          author='Shengen Hu',
          author_email='Tarelahu@gmail.com',
          url='https://Tarela@bitbucket.org/tarela/cam',
          package_dir=PKGdir,
          packages=PKGs,
          package_data={'CAMpipe': ['Config/CAM_template.conf',
                                  'Rscript/QCplots.r',
                                  'annotation/DHS_hg19.bed',
                                  'annotation/DHS_hg38.bed',
                                  'annotation/DHS_mm10.bed',
                                  'annotation/DHS_mm9.bed',
                                  'annotation/hg19.genome',
                                  'annotation/hg19_refgenes.txt',
                                  'annotation/hg38.genome',
                                  'annotation/hg38_refgenes.txt',
                                  'annotation/mm10.genome',
                                  'annotation/mm10_refgenes.txt',
                                  'annotation/mm9.genome',
                                  'annotation/mm9_refgenes.txt' ]},
          scripts=SCRIPT,
                    
          classifiers=[
        'Development Status :: version1.0 finish',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Artistic License',
        'Operating System :: POSIX',
        'Programming Language :: Python',
        'Topic :: pipeline',
        ],
          requires=[],
      )

    print 'Installation of CAM is DONE'
    


if __name__ == '__main__':
    main()

