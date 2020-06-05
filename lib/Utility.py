#!/usr/bin/env python
"""

Function declare:
 
def CMD                  (cmd)
def sp                   (cmd)
def sperr                (cmd)
def raise_error          ()
def pdf_name             (input_name)
def wlog                 (message,logfile)
def ewlog                (message,logfile)
def rwlog                (message,logfile)
def readAnnotation       (annotation)
def textformat           (inp)
def createDIR            (dirname)
def strlatexformat       (instr)
def transform_sam        (samfile,full_ext_bed,plus_1bp_bed,minus_1bp_bed,SDplus_1bp_bed,sample_reads,q30filter)
def transform_refgene    (refgene,ttsdis,outname)
"""
# import rpm
import subprocess
import sys
import os
import math
import random
import string
import twobitreader

def CMD(cmd):
    os.system(cmd)

def sp(cmd):
    '''
    Call shell cmd or software and return its stdout
    '''
    a=subprocess.Popen(cmd, stdout=subprocess.PIPE, shell='TRUE')
    ac = a.communicate()
    return ac
def sperr(cmd):
    '''
    Call shell cmd or software and return its stdout
    '''
    a=subprocess.Popen(cmd, stderr=subprocess.PIPE, shell='TRUE')
    ac = a.communicate()
    return ac
   
def raise_error():
    '''
    Raise an error messgae and exit
    '''
    print 'error occurs, check log file~!'
    sys.exit(1)


def detect_memory():
    meminfo={}#OrderedDict()
    try:
        with open('/proc/meminfo') as f:
            for line in f:
                meminfo[line.split(':')[0].strip()] = line.split(':')[1].strip()
        totalM = meminfo['MemTotal'].split()
        #freeM = meminfo['MemFree'].split()
        if totalM[1].lower() == "kb":
            try:
                totalM_G = int(totalM[0])/1e6
                return totalM_G
            except:
                return 'NA'
        else:
            return 'NA'    
    except:
        return 'NA'


def pdf_name(input_name):
    '''
    Change filename to pdf file name
    '''
    outputname = "_".join(input_name.split('.')[:-1])+".pdf"
    return outputname
    
def wlog(message,logfile):
    '''
    print a message and write the message to logfile
    '''
    print message
    os.system('echo "[LOG] %s " >> %s'%(message,logfile))
 
def flog(message,logfile):
    '''
    print a finished file to logfile
    '''
    print message
    os.system('echo "%s" >> %s'%(message,logfile))

def ewlog(message,logfile):
    '''
    print an error message and write the error message to logfile
    then exit Dr.seq
    error messages start with [ERROR]
    '''
    print "[ERROR] %s "%(message)
    os.system('echo "[ERROR] %s " >> %s'%(message,logfile))
    raise_error()
    
def rwlog(cmd,logfile) :
    '''
    print an (shell) command line and write the command line to logfile
    then conduct the command line
    command lines start with [CMD]
    '''
    print "[CMD] %s "%(cmd)
    os.system('echo "[CMD] %s " >> %s'%(cmd,logfile))
    CMD(cmd)
   
def check_filelist(filename,checklist):
    name = filename
    finished_files = open(name).readlines()
    finished_files = [i.strip() for i in finished_files]
    status = 1
    if checklist not in finished_files:
        status = 0
    #  status = 1 means the process has been finished.
    return status

    
def readAnnotation(annotation):
    '''
    read full annotation file and output as a dictionary 
    file format is fixed to UCSC full annotation format
    '''
    inf = open(annotation)
    outdict = {}
    for line in inf:
        ll = line.split()
        outdict[ll[1]] = ll[12]
    return outdict
    
     
def textformat(inp):
    '''
    transfer 1000000 to  1,000,000 for better visualization in output report
    '''
    o = ''
    comma = 0
    for i in (inp[::-1]):
        comma += 1
        o += i
        if comma%3 == 0:
            o += ','
    return o[::-1].strip(',')

def createDIR(dirname):
    '''
    check dir name and create new dir
    '''
    if not os.path.isdir(dirname):
        os.system('mkdir %s'%(dirname))

def strlatexformat(instr):
    outstr = instr.replace('_','\_')
    return(outstr)

def checkfa(fastafile):
    inf = open(fastafile)
    l1 = inf.readline()
    l2 = inf.readline().strip().upper()
    if not l1.startswith('>') :
        return 0
    if not set(l2).issubset(['A', 'C', 'G', 'N', 'T']):
        return 0
    return 1

def checkbed(bedfile):
    inf =open(bedfile)
    count = 0
    lendict = []
    for line in inf:
        count += 1
        if count == 10:
            break
        ll = line.strip().split('\t')
        if len(ll) < 3:
            return 0
        lendict.append(len(ll))
    inf.close()
    if len(set(lendict)) > 1:
        return 0
    else:
        return 1
    
def correct_genome_length(inGL,outGL):
    inf = open(inGL)
    outf = open(outGL,'w')
    for line in inf:
        ll = line.strip().split("\t")
        try :
            int(ll[1])
        except:
            continue
        outf.write(line)
    inf.close()
    outf.close()
    return

def fastq_reads_length(FQfile1):
    inf1 = open(FQfile1)
    read_length = []
    count = 0
    while 1:
        a=inf1.readline()
        if a.strip() == "" or count > 400:
            break
        elif count %4 == 1:
            RL = len(a.strip())
            read_length.append(RL)
        else:
            pass
        count += 1
    if len(set(read_length)) != 1 :
        return [read_length[0],"difflen"]
    else:
        return [read_length[0],"samelen"]
        

def transform_sam2bed(samfile,reads_bed,seqtype,q30filter,fraglen_range_raw):
    '''
    transform aligned samfile to bed file, only storage the fragment distribution for PE
    '''
    q30 = int(q30filter)
    inf = open(samfile)
    outbed = open(reads_bed,'w')
    fraglen_dict = {}
    totalN = 0
    fraglen_range = int(fraglen_range_raw)
    ### single end sam
    ## step 1 , transform to bed
    if seqtype == "SE":
        for line in inf:
            if line.startswith('@') or line.strip() == "":
                continue
            ll = line.strip().split("\t")
            chrom = ll[2]
            start = int(ll[3])-1
            seqlen = len(ll[9])
            end = start + seqlen
            txname = ll[0]
            mapQ = ll[4]
            flag = bin(int(ll[1]))[2:]
            if len(flag) < 5:
                if len(flag) < 3:
                    strand = "+"
                else:
                    if flag[-3] == "1":
                        continue
                    else:
                        strand = "+"
            else:
                if flag[-3] == "1":
                    continue
                elif flag[-5] == "1":
                    strand = "-"
                else:
                    strand = "+"
            if q30 == 1 and int(mapQ) < 30:
                continue
            newll = [chrom,start,end,txname,mapQ,strand]
            totalN += 1
            outbed.write("\t".join(map(str,newll))+"\n")
    ### paired end sam
    elif seqtype == 'PE':
        ### create fragment length dictionary
        for i in range(fraglen_range+1):
            fraglen_dict[i] = 0
            
        random.seed(1228)
        for line in inf:
            if line.startswith('@') or line.strip() == "":
                continue
            ll = line.strip().split("\t")
            pairlen = int(ll[8])
            ### ll[8] = 0 means its single end sam file or two pair in different chromosome
            if pairlen == 0:
                continue
            ### chrom2 != "=" , means two pair of reads don't map on same chromsome
            chrom2 = ll[6]
            if  chrom2 != '=':
                continue
            flag = bin(int(ll[1]))[2:]
            if len(flag) < 5:
                if len(flag) < 3:
                    strand = "+"
                else:
                    if flag[-3] == "1":
                        continue
                    else:
                        strand = "+"
            else:
                if flag[-3] == "1":
                    continue
                elif flag[-5] == "1":
                    strand = "-"
                else:
                    strand = "+"
            ### in paired end sam file we don't need - strand, cuz we can infer whole fragment from +strand reads
            if strand == '-':
                continue
            if int(pairlen) < 0 and strand == "+":
                continue #print 'WARNING: +reads but pairlen < 0',ll
            chrom = ll[2]
            start1 = int(ll[3])-1
            seqlen1 = len(ll[9])
            end1 = start1 + seqlen1
            strand1 = "+"
            start2 = int(ll[7])-1
            end2 = start2 + seqlen1
            strand2 = "-"
            txname = ll[0]
            mapQ = ll[4]
            
            if q30 == 1 and int(mapQ) < 30:
                continue
            if start1 > end2:
                continue
            if pairlen <= fraglen_range:
                fraglen_dict[pairlen] +=  1
            frag_center = int((start1+end2)/2)
            newstart1 = frag_center - 74
            newend2 = frag_center + 74
            
            RDnum = int(random.randint(1,2))
            if RDnum == 1:
                newll = [chrom,newstart1,newstart1+1,txname,mapQ,strand1]
            else:
                newll = [chrom,newend2-1,newend2,txname,mapQ,strand2]
            totalN += 1
            outbed.write("\t".join(map(str,newll))+"\n")
 
    outbed.close()
    return [totalN, fraglen_dict]


def transform_bed2bed(inputbed,reads_bed,seqtype,q30filter,fraglen_range_raw):
    '''
    transform bed file to bed file, only storage the fragment distribution for PE
    '''
    q30 = int(q30filter)
    outbed = open(reads_bed,'w')
    fraglen_dict = {}
    totalN = 0
    fraglen_range = int(fraglen_range_raw)
    
    if seqtype == "SE":
        inf = open(inputbed)
        for line in inf:
            ll = line.strip().split("\t")
            if len(ll) < 6:
                continue
            if not ll[5] in ['+','-']:
                continue
            mapQ = ll[4]
            if q30 == 1 and int(mapQ) < 30:
                continue
            totalN += 1
            outbed.write(line)
    ### paired end sam
    elif seqtype == 'PE':
        ### create fragment length dictionary
        for i in range(int(fraglen_range)+1):
            fraglen_dict[i] = 0
            
        random.seed(1228)
        inf = open(inputbed)
        ispair2 = 0
        for line in inf:
            ll = line.strip().split("\t")
            
            if ispair2 == 0:
                line1 = line
                chrom1 = ll[0]
                start1 = int(ll[1])
                end1 = int(ll[2])
                fullname1 = ll[3]
                mapQ1 = ll[4]
                strand1 = ll[5]
                ispair2 = 1
            else:
                line2 = line
                chrom2 = ll[0]
                start2 = int(ll[1])
                end2 = int(ll[2])
                fullname2 = ll[3]
                mapQ2 = ll[4]
                strand2 = ll[5]
                ispair2 = 0
                if chrom1 != chrom2:
                    continue
                if q30 == 1 :
                    if int(mapQ1) < 30 or int(mapQ2) < 30:
                        continue
                if strand1 == "+":
                    if strand2 != '-':
                        print 'WARNING: two pair of reads have same strand: ',fullname1,strand1,fullname2,strand2
                        continue
                    pairlen = end2 - start1
                    frag_center = int((start1+end2)/2)
                elif strand1 == '-':
                    if strand2 != '+':
                        print 'WARNING: two pair of reads have same strand: ',fullname1,strand1,fullname2,strand2
                        continue
                    pairlen = end1 - start2
                    frag_center = int((start2+end1)/2)
                else:
                    print 'WARNING: unexpected strand value: ',fullname1,strand1,fullname2,strand2
                    continue
                if pairlen <= 0:
                    continue
                if pairlen <= fraglen_range:
                    fraglen_dict[pairlen] +=  1

                newstart1 = frag_center - 74
                newend2 = frag_center + 74
                    
                RDnum = int(random.randint(1,2))
                if RDnum == 1:
                    newll = [chrom1,newstart1,newstart1+1,fullname1,mapQ1,strand1]
                else:
                    newll = [chrom1,newend2-1,newend2,fullname2,mapQ2,strand2]
#                RDnum = int(random.randint(1,2))
#                if RDnum == 1:
#                    #newll = [chrom1,start1,end1,fullname1,mapQ1,strand1]
#                    outbed.write(line1)
#                else:
#                    outbed.write(line2)
#                    #newll = [chrom2,start2,end2,fullname2,mapQ2,strand2]
                totalN += 1
                outbed.write("\t".join(map(str,newll))+"\n")
                # write only 1bp for the start or end for the shifted and extended 147bp fragments
 
    outbed.close()
    inf.close()
    return [totalN, fraglen_dict]


def bed2allbw(inputbed, genome_length, extbdg, extbw, centerwig, m1bed, sdp1bed, totalN, sample_reads, rpm):
    centerbw = centerwig[:-3]+'bw'
    ### calculate percentage of reads to be sample down
    p = 10000 * (int(sample_reads)*1.0/(totalN/2))
    # for PE data, inputbed storage only 1bp, for SE data ,inputbed storage the whole reads length
    ### generate different format of bed files
    # sort inputbed by chromsome 
    tmpsortbed = 'rawbed_sortbychrm_tmp.bed'
    tmpextbed = 'extbed_tmp.bed'
    tmpcenterbed = 'centerbed_tmp.bed'
    
    # read genome length to dict: GLdict
    inf = open(genome_length)
    GLdict = {}
    for line in inf:
        ll = line.strip().split("\t")
        GLdict[ll[0]] = int(ll[1])
    inf.close()
    
    sp('sort -k 1,1 %s > %s'%(inputbed,tmpsortbed))
    inf = open(tmpsortbed)
    # center 74bp fragments bed
    outf_fullext = open(tmpextbed,'w')
    # center 1bp
    outf_center = open(tmpcenterbed,'w')
    outf_minus1bp = open(m1bed,'w')
    # sample down plus reads (1bp), 10M
    outf_SDplus1bp = open(sdp1bed,'w')
    for line in inf:
        ll = line.strip().split("\t")
        chrom = ll[0]
        start = int(ll[1])
        end = int(ll[2])
        readsname = ll[3]
        mapQ = ll[4]
        strand = ll[5]
        
        if not GLdict.has_key(chrom):
            continue
        if start < 0 or end > GLdict[chrom]:
            continue
        
        if strand == '+':
            newll = [chrom,start+37,start+37*3,readsname,mapQ,strand]
            outf_fullext.write("\t".join(map(str,newll))+"\n")

            newll = [chrom,start+37*2,start+37*2+1,readsname,mapQ,strand]
            outf_center.write("\t".join(map(str,newll))+"\n")

            newll = [chrom,start,start+1,readsname,mapQ,strand]
            #outf_plus1bp.write("\t".join(map(str,newll))+"\n")
            # 
            if random.randint(1,10000) <= p:
                outf_SDplus1bp.write("\t".join(map(str,newll))+"\n")
        else:
            if (end-37*2) < 0:
                continue
            newll = [chrom,max(0,end-37*3),end-37,readsname,mapQ,strand]
            outf_fullext.write("\t".join(map(str,newll))+"\n")

            newll = [chrom,max(0,end-37*2-1),end-37*2,readsname,mapQ,strand]
            outf_center.write("\t".join(map(str,newll))+"\n")

            newll = [chrom,end-1,end,readsname,mapQ,strand]
            outf_minus1bp.write("\t".join(map(str,newll))+"\n")
            #if random.randint(1,10000) <= p:
            #    outf_SDminus1bp.write("\t".join(map(str,newll))+"\n")
    
    inf.close()
    outf_fullext.close()
    outf_center.close()
    #outf_plus1bp.close()
    outf_minus1bp.close()
    outf_SDplus1bp.close()
    #outf_SDminus1bp.close()
    
    ### generate bw
    # ext bw
    if int(rpm) == 1:
        # normalized profile by coverage
        # output: extbdg, extbw
        bed2bw(tmpextbed, genome_length, extbdg, extbw, totalN, "NA") 
        # bed2bw(tmpcenterbed, genome_length, "NA", centerbw, totalN, centerwig)
    else:
        bed2bw(tmpextbed, genome_length, extbdg, extbw, "NA", "NA")
        
    bed2bw(tmpcenterbed, genome_length,"NA", "NA", "NA", centerwig)
    bed2bw(tmpcenterbed, genome_length,"NA", centerbw, totalN, 'NA')
    
    # 1bp bw    
    #bed2bw(tmpplusbed, genome_length, p1bw, "NA")
    #bed2bw(tmpminusbed, genome_length, m1bw, "NA")

    ### clean tmp files
    #sp("rm %s"%(tmpsortbed)) 
    #sp("rm %s"%(tmpextbed))
    #sp("rm %s"%(tmpplusbed))
    #sp("rm %s"%(tmpminusbed))
    
    return




def bed2bw(inputbed, genome_length, outbdg, outbw, total_reads, outwig):

    ### bed -> bdg    
    if total_reads == "NA":
        # not scale
        sp("bedtools genomecov -bga -i %s -g %s > rawbdg_TMPbed2bw.bdg"%(inputbed, genome_length ))
    else:
        # scale by coverage(to 10M)
        CONSTANT = 1000000.0/total_reads
        sp("bedtools genomecov -bga -scale %s -i %s -g %s > rawbdg_TMPbed2bw.bdg"%( CONSTANT, inputbed, genome_length ))
        
    ### bdg clip
    sp('bedClip rawbdg_TMPbed2bw.bdg %s clipbdg_TMPbed2bw.bdg'%(genome_length))
    ### bdg sort
    sp('sort -k1,1 -k2,2n clipbdg_TMPbed2bw.bdg > clipsortbdg_TMPbed2bw.bdg')
    ### bdg -> bw
    if outbw != "NA":
        sp('bedGraphToBigWig clipsortbdg_TMPbed2bw.bdg %s %s'%(genome_length, outbw))
    if outwig != "NA":
        bdg2wig("clipsortbdg_TMPbed2bw.bdg", outwig, genome_length)
    ### clean tmp file
    sp('rm rawbdg_TMPbed2bw.bdg')
    sp('rm clipbdg_TMPbed2bw.bdg')
    if not outbdg == "NA":
        sp('mv clipsortbdg_TMPbed2bw.bdg %s'%outbdg)
    return 
 
def signal_on_aray(inputbdg, inputarray, outbw, genome_length):
    sp('bedtools intersect -a %s -b %s -u > bdgonarray_TMPsignalonarray.bdg'%(inputbdg, inputarray))
    sp('bedGraphToBigWig bdgonarray_TMPsignalonarray.bdg %s %s'%(genome_length, outbw))
    sp('rm bdgonarray_TMPsignalonarray.bdg')
    return
 
def bdg2wig(inputbdg, outwig, genome_length):
    inf = open(genome_length)
    GLdict = {}
    for line in inf:
        ll = line.strip().split("\t")
        GLdict[ll[0]] = int(ll[1])
    inf.close()
        
    inf = open(inputbdg)
    outf = open(outwig,'w')
    usedchr = []
    for line in inf:
        ll = line.strip().split("\t")
        chrm = ll[0]
        start = int(ll[1])
        end = int(ll[2])
        value = float(ll[3])
        
        if start == 0 and end == GLdict[chrm] and value == 0:
            continue
        
        region_len = end-start
        if not chrm in usedchr:
            newline = 'fixedStep chrom=%s start=1 step=10\n'%(chrm)
            outf.write(newline)
            usedchr.append(chrm)
            poscount = 0
            valuecount = 0
        while 1:
            if region_len >= (10-poscount):
                valuecount += (10-poscount)*value
                outf.write(str(1.0 * valuecount / 10)+"\n")
                region_len -= (10-poscount)
                poscount = 0
                valuecount = 0
            else:
                valuecount += region_len * value
                poscount += region_len           
                region_len = 0
            
            if region_len == 0:
                break
    inf.close()
    outf.close()
            
        
def fetch_total_reads_fromLog(bowtielog):
    useline = open(bowtielog).readline()
    return useline.split()[0]

def calculate_genome_length(genome_length_file):
    total_len = 0
    inf = open(genome_length_file)     
    for line in inf:
        ll = line.strip().split("\t")
        if len(ll) > 0:
            total_len += int(ll[1])
    inf.close()
    return total_len

def applyMAT(inline,COF,NORM,N):
    # gaussian smooth
    sig_sum = 0
    for i in range(N):
        sig_sum += inline[i]*COF[i]
    return (sig_sum/NORM)    

def modify_wig_signal(raw_wig, modify_wig, win_size, band_width):
    inf = open(raw_wig)
    outf=  open(modify_wig,'w')
    win_size = int(win_size) / 10
    band_width = int(band_width) / 10
    
    buf = [0.0]*win_size
    coeff = [0.0]*(2*win_size+1)
    nvc = len(coeff)
    norm = 0
    for i in range(2*win_size+1):
        coeff[i] =  pow( math.e ,  -0.5 * pow((i-win_size + 0.0) / band_width , 2) ) 
        norm += coeff[i]

    tmplist = [0.0]*(2*win_size+1)  
    # modifylist = [0.0]*(2*win_size+1)

    INIT = 0
    count= 0
    while 1:
        line = inf.readline()
        if line.strip() == "":
            break
        if line.startswith('fixedStep'):
            for i in range(win_size):
                # write for the last 7 bin of the previous chromosome
                if INIT == 1:
                    tmplist = tmplist[1:]+[0]
                    gs_sig = applyMAT(tmplist,coeff,norm,nvc)
                    count += 1
                    modify_sig = abs(float(gs_sig) - float(last_sig))
                    last_sig = gs_sig
                    outf.write(str(modify_sig)+"\n")
            INIT = 1
            outf.write(line)
            tmplist = [0.0]*(2*win_size+1)
            count = 0
        else:
            raw = float(line.strip())
            tmplist = tmplist[1:]+[raw]
            gs_sig = applyMAT(tmplist,coeff,norm,nvc)
            count += 1
            # for the first 7 bin, gs_signal = raw signal
            if count > win_size:
                if count == win_size +1:
                    modify_sig = 0
                else:
                    modify_sig = abs(float(gs_sig) - float(last_sig))
                    # modif_sig = difference the current gs_sig and previous base pair
                last_sig = gs_sig
                outf.write(str(modify_sig)+"\n")
    for i in range(win_size):
        # write for the last 7 bin of the last chromosome
        tmplist = tmplist[1:]+[0]
        gs_sig = applyMAT(tmplist,coeff,norm,nvc)
        count += 1
        modify_sig = abs(float(gs_sig) - float(last_sig))
        last_sig = gs_sig
        outf.write(str(modify_sig)+"\n")
    
    outf.close()
    inf.close()
    
def generate_position_signal(modify_wig, position_wig, win_size):
    inf = open(modify_wig)
    outf=  open(position_wig,'w') 
    win_size = int(win_size) / 10   
    tmplist = [0.0]*(2*win_size+1)
    total_sum_value = 0
    total_bin_number = 0  
    INIT = 0
    count= 0
    while 1:
        line = inf.readline()
        if line.strip() == "":
            break
        if line.startswith('fixedStep'):
            # the start of each chromosome
            for i in range(win_size):
                if INIT == 1:
                    tmplist = tmplist[1:]+[0]
                    count += 1
                    if tmplist[ win_size ] == max(tmplist):
                        thismax_x = count
                        thismax_y = tmplist[ win_size ] 
                        for j in range(thismax_x - lastmax_x):
                            pos_sig = 1.0*(j)*(thismax_y-lastmax_y)/(thismax_x-lastmax_x) + lastmax_y
                            total_sum_value += pos_sig
                            total_bin_number += 1
                            outf.write(str(pos_sig)+"\n")
                        lastmax_x = thismax_x
                        lastmax_y = thismax_y
            # assign number for the first step in each chromosome
            INIT = 1
            outf.write(line)
            tmplist = [0.0]*(2*win_size+1)
            count = 0
        else:
            raw = float(line.strip())
            tmplist = tmplist[1:]+[raw]
            count += 1
            if count > win_size:
                if count == win_size +1:
                    lastmax_x = count
                    lastmax_y = tmplist[ -1*win_size -1]
                else:
                    if tmplist[ -1*win_size -1] == max(tmplist):
                        thismax_x = count
                        thismax_y = tmplist[ -1*win_size -1]
                        for j in range(thismax_x - lastmax_x):
                            pos_sig = 1.0*(j)*(thismax_y-lastmax_y)/(thismax_x-lastmax_x) + lastmax_y
                            total_sum_value += pos_sig
                            total_bin_number += 1
                            outf.write(str(pos_sig)+"\n")
                        lastmax_x = thismax_x
                        lastmax_y = thismax_y
    
    for i in range(win_size):
        tmplist = tmplist[1:]+[0]
        count += 1
        if tmplist[ win_size ] == max(tmplist):
            thismax_x = count
            thismax_y = tmplist[ win_size ] 
            for j in range(thismax_x - lastmax_x):
                pos_sig = 1.0*(j)*(thismax_y-lastmax_y)/(thismax_x-lastmax_x) + lastmax_y
                total_sum_value += pos_sig
                total_bin_number += 1
                outf.write(str(pos_sig)+"\n")
            lastmax_x = thismax_x
            lastmax_y = thismax_y
    
    outf.close()
    inf.close()
    bgvalue =  1.0 * total_sum_value / total_bin_number 
    return bgvalue    
         
def make_array(position_wig, array_bed, bgvalue):
    inf = open(position_wig)
    outf = open(array_bed,'w')
    BG = float(bgvalue)
    bin_count = 0
    onArray = 0
    nameN = 0
    for line in inf:
        if line.startswith('fixedStep'):
            if onArray == 1:
                onArray = 0
                end = bin_count
                if end != start:
                    foldchange = sumValue/(end-start)/BG
                    nameN += 1
                    newll =  [chrm,(start-1)*10,(end-1)*10,'array'+str(nameN),foldchange]
                    outf.write("\t".join(map(str,newll))+"\n")
            bin_count = 0
            chrm = line.split()[1].strip().split("=")[1]
            continue
        bin_count += 1
        value = float(line.strip())
        if value > BG:
            if onArray == 0:
                onArray = 1
                start = bin_count
                sumValue = value
            else:
                sumValue += value
        else:
            if onArray == 0:
                pass
            else:
                onArray = 0
                end = bin_count 
                foldchange = sumValue/(end-start)/BG
                nameN += 1
                newll = [chrm,(start-1)*10,(end-1)*10,'array'+str(nameN),foldchange]
                outf.write("\t".join(map(str,newll))+"\n")
    inf.close()
    outf.close()


def sum_bed_length(infile):
    inf = open(infile)
    count = 0
    for line in inf:
        ll   = line.split()
        count += int(ll[2]) - int(ll[1])
    inf.close()
    return (count)

def generate_utr_overlap(inbed, gene_annotation, promoter_dis):
    inf = open(gene_annotation)
    outf = open('tmp_utr_region.TMPbed','w')
    for line in inf:
        if line.startswith('#'):
            continue
        ll = line.strip().split('\t')
        txname = ll[1]
        chrm  = ll[2]
        strand = ll[3]
        start = ll[4]
        end = ll[5]
        if strand == "+":
            tss = int(start)
            newll = [chrm, tss , tss + int(promoter_dis)]
        else:
            tss = int(end)
            newll = [chrm, tss - int(promoter_dis), tss ]

        if tss < int(promoter_dis):
            continue
        #newll = [chrm, tss - int(promoter_dis), tss + int(promoter_dis)]
        outf.write("\t".join(map(str,newll))+"\n")
    inf.close()
    outf.close()
    array_utr_ov = sp('bedtools intersect -a %s -b %s -u | wc -l'%(inbed, 'tmp_utr_region.TMPbed'))[0].split()[0]
    sp('sort -k 1,1 -k 2,2g tmp_utr_region.TMPbed | bedtools merge -i - > tmp_utr_regionMerge.TMPbed')
    total_utr_len = sum_bed_length('tmp_utr_regionMerge.TMPbed')
    sp('rm tmp_utr_region.TMPbed')
    sp('rm tmp_utr_regionMerge.TMPbed')
    
    return [array_utr_ov,total_utr_len]

def generate_DHS_overlap(inbed, DHS):
    total_DHS_len = sum_bed_length(DHS)
    array_DHS_ov = sp('bedtools intersect -a %s -b %s -u | wc -l'%(inbed, DHS))[0].split()[0]
    return [array_DHS_ov, total_DHS_len]


def integrate_overlap(raw_gene_bed, raw_ov_array, summary_ov_array):
    inf = open(raw_gene_bed)
    raw_dict = {}
    for line in inf:
        raw_dict[line.strip()] = []
    inf.close()
    inf = open(raw_ov_array)
    for line in inf:
        ll = line.strip().split("\t")
        anno = "\t".join(ll[:5])
        arraylen = int(ll[7]) - int(ll[6])
        raw_dict[anno].append(arraylen)
    inf.close()
    outf = open(summary_ov_array,'w')
    for ANNO in sorted(raw_dict.keys()):
        if len(raw_dict[ANNO]) == 0:
            newll = ANNO.split("\t")[:4] + ["-1", ANNO.split("\t")[4]]
        else: 
            #if len(raw_dict[ANNO]) > 1:
            #    print ANNO
            newll = ANNO.split("\t")[:4] + [str( max(raw_dict[ANNO]) ), ANNO.split("\t")[4]]
        outf.write("\t".join(newll)+"\n")
    outf.close()

def generate_geneLevel_arrayAnnotation(inbed, gene_annotation, summary_ov_array):
    inf = open(gene_annotation)
    outf = open('tmp_promoter_region.TMPbed','w')
    for line in inf:
        if line.startswith('#'):
            continue
        ll = line.strip().split('\t')
        txname = ll[1]
        chrm  = ll[2]
        strand = ll[3]
        start = ll[4]
        end = ll[5]
        if strand == "+":
            tss = int(start)
            newll = [chrm, tss - 3000 , tss + 3000, txname, strand]
        else:
            tss = int(end)
            newll = [chrm, tss - 3000 , tss + 3000, txname, strand]

        if tss < 3000:
            continue

        outf.write("\t".join(map(str,newll))+"\n")
    inf.close()
    outf.close()
    sp('bedtools intersect -a %s -b %s -wo > %s'%( 'tmp_promoter_region.TMPbed', inbed, 'tmp_promoter_overlap.TMPbed'))
    integrate_overlap('tmp_promoter_region.TMPbed', 'tmp_promoter_overlap.TMPbed', summary_ov_array)

    sp('rm tmp_promoter_region.TMPbed')
    sp('rm tmp_promoter_overlap.TMPbed')



def tssprofile(inputbw,outmatrix,tssannotation,UP,DOWN):
    inf = open(tssannotation)
    outf = open(outmatrix,'w')
    Ndatapoints = ( int(UP)+int(DOWN) ) /10
#    Gene_count = 0
    sum_result = [0]*Ndatapoints
    ave_result = [0]*Ndatapoints
    for line in inf:
        if line.startswith('#') or line.strip() == "":
            continue
        ll = line.strip().split('\t')
        txname = ll[1]
        chrm  = ll[2]
        strand = ll[3]
        start = ll[4]
        end = ll[5]
        if strand == "+":
            tss = int(start)
            if tss < int(UP):
                continue
            cmd = 'bigWigSummary %s %s %s %s %s'%(inputbw,chrm,tss-int(UP),tss+int(DOWN),Ndatapoints)
            
            result =  sp(cmd)[0].strip().split("\t")
        else:
            tss = int(end)
            if tss < int(DOWN):
                continue
            cmd = 'bigWigSummary %s %s %s %s %s'%(inputbw,chrm,tss-int(DOWN),tss+int(UP),Ndatapoints)
            result = sp(cmd)[0].strip().split("\t")[::-1]
        if len(result) == 0 or "n/a" in result:
            pass
        elif len(result) == 1 and result[0] == "":
            pass 
        else:
            outf.write("\t".join(map(str, [chrm, start, end, txname, '.', strand]+ result))+"\n")

    inf.close()
    outf.close()

def customprofile(inputbw,outmatrix,customRegion,plotdis):
    inf = open(customRegion)
    outf = open(outmatrix,'w')
    Ndatapoints = ( 2* int(plotdis) ) /10
    count = 0
    for line in inf:
        if line.strip() == "":
            continue
        count += 1
        if count > 50000:
            break
        ll = line.strip().split('\t')
        chrm  = ll[0]
        start = ll[1]
        end = ll[2]
        center = (int(start) + int(end) )/2 
        if center < int(plotdis) :
            continue
        
        cmd = 'bigWigSummary %s %s %s %s %s'%(inputbw,chrm,center-int(plotdis),center+int(plotdis),Ndatapoints)
        raw_result =  sp(cmd)[0].strip().split("\t")
        if len(ll) >= 6 and ll[5] == "-":
            result = raw_result[::-1]
        else:
            result = raw_result
        if len(ll) >= 6:
            anno = ll[:6]
        else:
            anno = ll + ["."]*(6-len(ll))
        if len(result) == 0 or "n/a" in result:
            pass
        elif len(result) == 1 and result[0] == "":
            pass 
        else:
            outf.write("\t".join(map(str, anno + result))+"\n")

    inf.close()
    outf.close()



def ATprofile(inputbed,genome2bit):
    inf = open(inputbed)
    genome = twobitreader.TwoBitFile(genome2bit) 
    ATfrac = [0]*140
    ave_ATfrac = [0]*140
    reads_count = 0
    for line in inf:
        ll = line.strip().split('\t')
        #chrom = ll[0]
        if ll[5] == "+":
            seq = genome[ll[0]][(int(ll[1])+4):(int(ll[1]) + 145)].upper()
        else:
            #transtab = string.maketrans("ACGTNX","TGCANX")
            seq = genome[ll[0]][(int(ll[2]) - 145):(int(ll[2])-4)].upper()[::-1]#.translate(transtab)[::-1]
        if len(seq) != 141:
            continue        
        reads_count += 1
        for i in range(140):
            if seq[i:(i+2)] in ["AA","TT","AT","TA"]:
                ATfrac[i] += 1
    
    for i in range(140):
        ave_ATfrac[i] = (ATfrac[i]*1.0) / (reads_count*1.0)
    inf.close()
    return ave_ATfrac

def judgeY(indata):
    if int(indata) > 0:
        return 1
    else:
        return 0
    
def reads_reads_distance(SDplus1bp, minus1bp, genome_length, reads_distance_range):
    inf = open(genome_length)
    GLdict = {}
    for line in inf:
        ll = line.strip().split("\t")
        GLdict[ll[0]] = int(ll[1])
    inf.close()
    infplus = open(SDplus1bp)
    infminus = open(minus1bp)
    reads_count = 0
    occurance = [0] * int(reads_distance_range)
    ave_occurance = [0] * int(reads_distance_range)
    usedchr = []
    for line in infplus:
        ll = line.strip().split("\t")
        if not ll[0] in usedchr:
            usedchr.append(ll[0])
            tmpCHR = ll[0]
            tmpVector = [0]*(GLdict[tmpCHR]+1)
            while 1:
                LL = infminus.readline().strip().split("\t")
                if len(LL) == 0:
                    break
                if LL[0] == tmpCHR:
                    tmpVector[int(LL[2])] = 1
                else:
                    break
        start = int(ll[1])
        if start + int(reads_distance_range) < len(tmpVector):
            rawdata = tmpVector[start:(start+ int(reads_distance_range))]
            reads_count += 1
            for i in range(int(reads_distance_range)):
                if int(rawdata[i]) > 0:
                    occurance[i] += 1
    inf.close()
    for i in range(int(reads_distance_range)):
        ave_occurance[i] = (occurance[i]*1.0) / (reads_count*1.0)
    return ave_occurance
     

