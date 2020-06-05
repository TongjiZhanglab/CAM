#!/usr/bin/env python
# ------------------------------------
"""
Function declare:

def gen_conf (name)
def read_conf(configfile)

"""
# -----------------------------------


import os,sys
import ConfigParser

import CAMpipe

### config template
CONFIG_TEMPLATE = os.path.join(CAMpipe.__path__[0], "Config/CAM_template.conf")
### generate a config
def gen_conf(CONF_name):
    '''
    Generate a config file
    ''' 
    inf = open(CONFIG_TEMPLATE)
    if not CONF_name.endswith('.conf'):
        CONF_name += '.conf'
    outf = open(CONF_name,'w')
    for line in inf:
        outf.write(line)
    outf.close()
    inf.close()

### read config
def read_conf(conf_file):
    '''
    Read config file and return a dict containing all infomation
    '''
    conf_dict = {}
    cf = ConfigParser.SafeConfigParser()
    cf.read(conf_file)
    #section_list = sorted( cf.sections() , key=lambda x:int(x.lstrip('Step').split("_")[0]))
    for st in cf.sections():
        conf_dict[st]={}
        for item in cf.items(st):
            conf_dict[st][item[0]]=item[1]
    return conf_dict

### generate a config file in simple mode
#        make_conf(args.barcode,args.reads,args.name,args.fover,args.CBL,args.UMIL,args.RF,args.P)

#def make_conf(barcode_file,reads_file,outname,fover,cellbarcodeL,umiL,geneanno,P,mapindex,checkmem,maptool,select_cell,remove_lowdup):
def make_conf(inputa,inputb,seqtype,outname,customregion,fover,Species,P,genome_fasta,mapindex):  
    name = outname
    # if os.path.isfile(name+'.conf') and  not fover :
    #     print 'config file "%s" using same name exist , choose other name or add -f to overwrite'%(name+".conf")
    #     sys.exit(1)

    # not exist or for over write:
    if not os.path.isfile(name+'.conf') or fover :
        inf = open(CONFIG_TEMPLATE)
        outf = open(name+".conf",'w')
        #print inputb is None
        for line in inf:
            if line.startswith('inputa ='):
                newline = 'inputa = ' + inputa + '\n'
            elif line.startswith('inputb ='):
                if inputb is None:
                    newline = line
                else:
                    newline = 'inputb = ' + inputb + '\n'
            elif line.startswith('seqtype ='):
                newline = 'seqtype = ' + seqtype + '\n'
            elif line.startswith('outputdirectory ='):
                newline = 'outputdirectory = ' + name + '\n'
            elif line.startswith('outname ='):
                newline = 'outname = ' + name + '\n'
            elif line.startswith('species ='):
                if Species:
                    newline = 'species = ' + Species + '\n'
                else:
                    newline = line	                
            elif line.startswith('customregion ='):
                if customregion:
                    newline = 'customregion = ' + customregion + '\n'
                else:
                    newline = line	
            elif line.startswith('genome_fasta ='):
                if genome_fasta:
                    newline = 'genome_fasta = ' + genome_fasta + '\n'
                else:
                    newline = line
            elif line.startswith('mapindex ='):
                if mapindex:
                    newline = 'mapindex = ' + mapindex + '\n'
                else:
                    newline = line
            elif line.startswith('mapping_p ='):
                newline = 'mapping_p = ' + str(P) + '\n'
            else:
                newline = line
            outf.write(newline)
        inf.close()
        outf.close()
    else:
        conf_dict = read_conf(name+'.conf')
        print name+'.conf'
        status = 1
        # default: should not change parameters
        # if conf_dict['General']['inputa']!= inputa or check_inputb!= inputb or conf_dict['General']['seqtype']!= seqtype or conf_dict['General']['outputdirectory']!= name or conf_dict['General']['outname']!= name or conf_dict['Step1_preprocess']['species']!= species or conf_dict['General']['customregion']!= customregion or conf_dict['Step1_preprocess']['genome_fasta']!= genome_fasta or conf_dict['Step1_preprocess']['mapindex']!= mapindex:
        if conf_dict['General']['outputdirectory']!= name:
            print 'Please keep the previous output directory'
            sys.exit(1)
            if conf_dict['Step1_preprocess']['species']!= Species:
                print 'Please keep the previous species'
                sys.exit(1)
                if conf_dict['General']['customregion']!= customregion:
                    print 'Please keep the previous customregion'
                    sys.exit(1)
                    if conf_dict['Step1_preprocess']['genome_fasta']!= genome_fasta:
                        print 'Please keep the previous genome_fasta'
                        sys.exit(1)
                        if conf_dict['Step1_preprocess']['mapindex']!= mapindex:
                            print 'Please keep the previous mapindex'
                            sys.exit(1)
            # print 'Please keep the previous parameters, or choose the other name'
            # sys.exit(1)
        else:
            print 'Correct parameters'
            # sys.exit(1)

        # the conf file has already exists, check the information in the exist conf file with 
    return  

def check_filelist(filename,checklist):
    name = filename
    finished_files = open(name).readlines()
    finished_files = [i.strip() for i in finished_files]
    status = 1
    for each in checklist:
        if each not in finished_files:
            status = 0
    #  status = 1 means the process has been finished.
    return status

