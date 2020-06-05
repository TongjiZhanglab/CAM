#!/usr/bin/env python
# coding: utf-8
# ------------------------------------
# Python Modual
# ------------------------------------

import os
import sys
import time
import string

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
                                   textformat,
                                   strlatexformat,check_filelist)
# --------------------------
# main 
# --------------------------
def step5_summary(conf_dict,logfile,filelist):
    '''
    analysis part
    mainly Rscript
    dimentional reduction + clustering
    '''
    # start
    # create section for 
    
    summarydir = conf_dict['General']['outputdirectory'] + 'summary/'
    mapping_dir = conf_dict['General']['outputdirectory'] + 'preprocess/'
    arraydir = conf_dict['General']['outputdirectory'] + 'nucarray/'
    qcdir = conf_dict['General']['outputdirectory'] + 'QC/'
    createDIR(summarydir)
    os.chdir(summarydir)
    
    result_folder = conf_dict['General']['outputdirectory'] + 'summary/plots/'
    createDIR(result_folder)
    # confirm result path
    conf_dict['Step1_preprocess']['profilebw'] = mapping_dir + conf_dict['General']['outname'] + '_profile.bw'
    if conf_dict['General']['seqtype']=='PE':
        conf_dict['Step1_preprocess']['mapped_bed'] = mapping_dir + conf_dict['General']['outname'] + '.bed'
    else:
        conf_dict['Step1_preprocess']['mapped_bed'] = conf_dict['Step1_preprocess']['bed']
    conf_dict['Step1_preprocess']['mapped_bb'] = mapping_dir + conf_dict['General']['outname'] + '.bb'
    conf_dict['Step2_QC']['ave_tssprofile'] = qcdir + conf_dict['General']['outname'] + '_Tss_profile.txt'
    conf_dict['Step2_QC']['fraglen'] = qcdir + conf_dict['General']['outname'] + '_fraglen.txt'
    conf_dict['Step3_nucarray']['arrayselect'] = arraydir + conf_dict['General']['outname'] + '_Nucleosome_Array.bed'
    conf_dict['Step3_nucarray']['gene_array_anno'] = arraydir + conf_dict['General']['outname'] + '_geneLevel_nucarrayAnnotation.bed'
    conf_dict['Step3_nucarray']['profilebw_onarray'] = arraydir + conf_dict['General']['outname'] + "_profile_on_Nucleosome_Array.bw"
    conf_dict['Step3_nucarray']['positionbw_center'] = arraydir + conf_dict['General']['outname'] + '_center_position.bw'
    conf_dict['Step1_preprocess']['centerbw'] = mapping_dir + conf_dict['General']['outname'] + '_center.bw'
    # conf_dict['Step1_preprocess']['bowtieout'] = mapping_dir + conf_dict['General']['outname'] + '.bowtieout'
    ## collect results

    wlog('collect output files',logfile)
    rwlog('mv %s %s'%(conf_dict['Step1_preprocess']['profilebw']  , summarydir),logfile)    
    rwlog('mv %s %s'%(conf_dict['Step1_preprocess']['mapped_bed']  , summarydir),logfile)
    if os.path.exists(conf_dict['Step1_preprocess']['mapped_bb']) :
        rwlog('mv %s %s'%(conf_dict['Step1_preprocess']['mapped_bb']  , summarydir),logfile)    
    rwlog('mv %s %s'%(conf_dict['Step2_QC']['ave_tssprofile']  , summarydir),logfile)
    rwlog('mv %s %s'%(conf_dict['Step2_QC']['fraglen']  , summarydir),logfile)
    if int(conf_dict['Step2_QC']['plotcustom']) == 1:
        rwlog('mv %s %s'%(conf_dict['Step2_QC']['custom_profile']  , summarydir),logfile)
    rwlog('mv %s %s'%(conf_dict['Step3_nucarray']['arrayselect']  , summarydir),logfile)
    rwlog('mv %s %s'%(conf_dict['Step3_nucarray']['gene_array_anno']  , summarydir),logfile)
    rwlog('mv %s %s'%(conf_dict['Step3_nucarray']['profilebw_onarray']  , summarydir),logfile) 
    rwlog('mv %s %s'%(conf_dict['Step3_nucarray']['positionbw_center']  , summarydir),logfile) 
    rwlog('mv %s %s'%(conf_dict['Step1_preprocess']['centerbw']  , summarydir),logfile) 
    

    # plot material
    rwlog('cp %s %s'%(qcdir + conf_dict['General']['outname'] + '_Features.txt'  , result_folder),logfile)
    rwlog('cp %s %s'%(qcdir + conf_dict['General']['outname'] +"_seqCov.pdf"  , result_folder),logfile)
    rwlog('cp %s %s'%(qcdir + conf_dict['General']['outname'] +"_profile.pdf"  , result_folder),logfile)
    rwlog('cp %s %s'%(qcdir + conf_dict['General']['outname'] +"_ATfrac.pdf"  , result_folder),logfile)
    rwlog('cp %s %s'%(qcdir + conf_dict['General']['outname'] +"_fraglen.pdf"  , result_folder),logfile)
    rwlog('cp %s %s'%(qcdir + conf_dict['General']['outname'] +"_Nucdep.pdf"  , result_folder),logfile)
    rwlog('cp %s %s'%(qcdir + conf_dict['General']['outname'] +"_Nucfuzz.pdf"  , result_folder),logfile)
    rwlog('cp %s %s'%(qcdir + conf_dict['General']['outname'] +"_dhs.pdf"  , result_folder),logfile)
    rwlog('cp %s %s'%(qcdir + conf_dict['General']['outname'] +"_utr.pdf"  , result_folder),logfile)

    os.chdir(result_folder)


    # read feature
    inf = open(conf_dict['General']['outname'] + '_Features.txt')
    for line in inf:
        if line.startswith('seq_coverage'):
            seq_coverage = round(float(line.strip().split()[1]),2)
        if line.startswith('rot_score'):
            rot_score = round(float(line.strip().split()[1]),4)
        if line.startswith('nuclen'):
            nuclen = int(float(line.strip().split()[1]))
        if line.startswith('NFRscore'):
            NFRscore = round(float(line.strip().split()[1]),4)
        if line.startswith('PSarray'):
            PSarray = round(float(line.strip().split()[1]),4)
        if line.startswith('array_on_utr'):
            array_on_utr = float(line.strip().split()[1])
        if line.startswith('array_on_DHS'):
            array_on_DHS = float(line.strip().split()[1])
        if line.startswith('array_num'):
            array_num = float(line.strip().split()[1])
        if line.startswith('total_utr_length'):
            total_utr_length = float(line.strip().split()[1])
        if line.startswith('total_DHS_length'):
            total_DHS_length = float(line.strip().split()[1])
        if line.startswith('effective_gs'):
            effective_gs = float(line.strip().split()[1])
        if line.startswith('enrichment_on_UTR'):
            UTR_fold = float(line.strip().split()[1])
        if line.startswith('enrichment_on_DHS'):
            DHS_fold = float(line.strip().split()[1])
    inf.close()
    if rot_score < 0.08:
        rot_judge = "Fail"
    else:
        rot_judge = "Pass"
    if nuclen < 140 or nuclen > 155:
        nuclen_judge = "Fail"
    else:
        nuclen_judge = "Pass"
    if NFRscore >= 0.4:
        NFR_judge = "Pass"
    else:
        NFR_judge = "Fail"
    if PSarray <= 0.4:
        PSarray_judge = "Pass"
    else:
        PSarray_judge = "Fail"       
    if UTR_fold<1:
        UTR_judge = "Fail"
    else:
        UTR_judge = "Pass"
    if DHS_fold<2:
        DHS_judge = "Fail"
    else:
        DHS_judge = "Pass"

    wlog('generate qc documents',logfile)
    ### initiate 
    QCdoc = """\documentclass[11pt,a4paper]{article}
\usepackage{tabularx}
\usepackage[english]{babel}
\usepackage{array}
\usepackage{graphicx}
\usepackage{color}
\DeclareGraphicsExtensions{.eps,.png,.pdf,.ps}
\\begin{document}
\\title{QC and analysis reports for MNase-seq data : %s}

\\vspace{-1cm}
\maketitle
\\tableofcontents
\\newpage
\\newpage
\section{Data description}
\\begin{quotation}
Table 1 mainly describe the input file and mapping and analysis parameters.
\end{quotation}
\\begin{table}[h]
\caption{Data description}\label{bstable}
\\begin{tabularx}{\\textwidth}{ |X|l| }

"""%(strlatexformat(conf_dict['General']['outname']))
    ### table1 prepare parameter
    if conf_dict['General']['seqtype'] == "PE":
        Seqtype = "Paired end"
        if conf_dict['General']['format'].upper() == "FASTQ":
            inputb = conf_dict['General']['inputb'].split("/")[-1]
        else:
            inputb = "NA"
    else:
        inputb = "NA"
        Seqtype = "Single end"
    if int(conf_dict['Step1_preprocess']['q30filter']) == 1:
        q30filter = "True"
    if int(conf_dict['Step2_QC']['plotcustom']) == 1:
        customRegion = conf_dict['General']['customregion'].split("/")[-1]
    else:
        customRegion = "NA"

    QCdoc += """      
\hline
parameter & value  \\\\
\hline
output name & %s \\\\
\hline
input file A & %s \\\\
\hline
input file B & %s \\\\
\hline
input format & %s  \\\\
\hline
sequencing type &  %s \\\\
\hline
genome version (species) & %s \\\\
\hline
Q30 filter mapped reads & %s \\\\
\hline
custom region & %s \\\\
\hline
\end{tabularx}
\end{table}
"""%(strlatexformat(conf_dict['General']['outname']),
     strlatexformat(conf_dict['General']['inputa'].split("/")[-1]),
     strlatexformat(inputb),
     conf_dict['General']['format'].upper(),
     Seqtype,
     conf_dict['Step1_preprocess']['species'],
     q30filter,
     strlatexformat(customRegion)
     )

    ###  QC component
    QCdoc += """
\\newpage
\\newpage
\section{QC component}
we calculated three key measurements: 1) sequencing coverage, 2) AA/TT/AT dinucleotide frequency and 3) nucleosomal DNA length distribution.
\subsection{Sequencing coverage}
\\begin{quotation}
Sequencing coverage provides a direct measurement of the resolution of two features of nucleosome organization, i.e. occupancy and positioning (Struhl and Segal, 2013). Sequencing coverage is defined as: (Number of reads * 194bp)/(Effective genome size).  "Number of reads" is the number of mappable reads after MAPQ filtering (for single end data, for paired end it's the number of fragment).  "194bp" is the total length of nucleosome and linker estimated from historical data. "Effective genome size" is defined as 2.7e9 bps for humans and 1.87e9 bps for mice. Below we plotted the distribution of sequencing coverage of historical data; the sequencing coverage of input data was marked by vertical line: %s. 
\end{quotation}
\\begin{figure}[h]
        \caption{Sequencing coverage} \label{fig:profileunion}
        \setlength{\\abovecaptionskip}{0pt}
        \setlength{\\belowcaptionskip}{10pt}
        \centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\end{figure}

\\newpage
\\newpage
\subsection{AA/TT/AT di-nucleotide frequency}
\\begin{quotation}
The 10-base AA/TT/AT periodicity in nucleosomal DNA provides a measurement of nucleosome rotational positioning, which has been shown to be influenced by DNA sequence (Satchwell, et al., 1986). Mappable reads were sampled down to 10 million and were extended to 147bp in their 3'end direction. Then the aggregate AA/TT/AT di-nucleotide frequency across 4th - 143th bp of the extended reads was calculated (right). We conducted a Fourier transform on the aggregate frequency and used the energy of 10-bp periodicity (defined as rotational score) to show the extent the MNase-seq reads reflect nucleosome organization. Sample with rotational score greater than 0.08 was defined as "Pass" in this measurement, otherwise it's defined as "Fail". The cutoff 0.08 was determined from the distribution of rotational scores from all historical data (left, vertical line marked the rotational score input sample: %s [%s]).  
\end{quotation}
\\begin{figure}[h]
        \caption{AA/TT/AT di-nucleotide frequency} \label{fig:profileunion}
        \setlength{\\abovecaptionskip}{0pt}
        \setlength{\\belowcaptionskip}{10pt}
        \centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\end{figure}

\\newpage
\\newpage
\subsection{Nucleosomal DNA length distribution}
\\begin{quotation}
Nucleosomal DNA length distribution (refer to fragment length or MNase library size) is closely related and thus can reflect the degree of MNase digestion. For paired end sample, fragment length distribution from all mappable fragments was used directly to infer the nucleosomal DNA length distribution. For single end sample, we calculated a start-to-end distance to estimate the nucleosome length distribution: mappable reads were sampled down to 10 million and then we calculated the distribution of the distance from 5'end of each plus strand read to all 5'end of minus strand reads within 250bp downstream. Duplicate reads were discarded in this calculation. After the distribution of nucleosomal DNA length was generated, the length with highest frequency was defined as the estimated nucleosomal DNA length of the input sample (for both paired end and single end, left). Sample with nucleosomal DNA length within 140bp - 155bp was defined as "Pass", otherwise it's defined as "Fail". The cutoff was determined from the distribution of nucleosomal DNA length from all historical data (left, vertical line marked the nucleosomal DNA length of input sample: %s [%s]).  
\end{quotation}
\\begin{figure}[h]
        \caption{Nucleosomal DNA length distribution} \label{fig:profileunion}
        \setlength{\\abovecaptionskip}{0pt}
        \setlength{\\belowcaptionskip}{10pt}
        \centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\end{figure}
"""%(str(seq_coverage),
     conf_dict['General']['outname'] +"_seqCov.pdf",
     str(rot_score),
     str(rot_judge),
     conf_dict['General']['outname'] +"_ATfrac.pdf",
     str(nuclen),
     str(nuclen_judge),
     conf_dict['General']['outname'] +"_fraglen.pdf",
    )

    if int(conf_dict['Step2_QC']['plotcustom']) == 1:
        customTEXT = " and the custom regions: %s"%(strlatexformat(customRegion))
        siteheat_filenames = strlatexformat(conf_dict['Step2_QC']['ave_tssprofile'].split("/")[-1] + ' \& '+  conf_dict['Step2_QC']['custom_profile'].split("/")[-1])
    else:
        customTEXT = ""
        siteheat_filenames = strlatexformat(conf_dict['Step2_QC']['ave_tssprofile'].split("/")[-1])
    QCdoc += """
\\newpage
\\newpage
\subsection{Nucleosome profile on potential functional regions}
\\begin{quotation}
CAM generated the average curve and the heatmap of nucleosome organization on promoter regions%s in 10bp resolution. Signal from minus strand regions were reversed in both heatmap and aggregate curve. The signal for each regions were also outputted as matrix: %s.  
\end{quotation}
\\begin{figure}[h]
        \caption{Nucleosome profile on potential functional regions} \label{fig:profileunion}
        \setlength{\\abovecaptionskip}{0pt}
        \setlength{\\belowcaptionskip}{10pt}
        \centering
        {\includegraphics[width=0.8\\textwidth]{%s}}
\end{figure}
\\newpage
\\newpage
\subsection{Nucleosome depletion level and nucleosome fuzziness around TSS}
\\begin{quotation}
Based on nucleosome profiles on promoters, CAM generated two scores to describe the nucleosome positioning on promoters. First, nucleosome depletion level described the fold change of the MNase-seq signal of nucleosome free regions compared to the +1 nucleosome and -1 nucleosome. The higher the nucleosome depletion level is, the deeper the nucleosome free region is. Lower nucleosome depletion level associated with weak or none nucleosome free regions, which may indicate reads from open chromatins. Samples with nucleosome depletion level higher than 0.4 was defined as "Pass", otherwise it's defined as "Fail".The cutoff was determined based on the distribution of nucleosome depletion level from all historical data (left, vertical line marked the nucleosome depletion level of input sample: %s [%s]).
Next, nucleosome fuzziness downstream TSS defined whether clear nucleosome positioning pattern was observed from downstream promoters. The nucleosome fuzziness was calculated by the coefficient of variance (CV) of the linker length between the +1, +2, +3 and +4 nucleosomes. The lower the nucleosome fuzziness is, the better nucleosome positioning was observed on promoters. Samples with nucleosome fuzziness lower than 0.4 was defined as "Pass", otherwise it's defined as "Fail".The cutoff was determined based on the distribution of nucleosome fuzziness scores from all historical data (right, vertical line marked the nucleosome fuzziness of input sample: %s [%s]).
\end{quotation}
\\begin{figure}[h]
    \\begin{minipage}[t]{0.5\linewidth}
        \centering
        \includegraphics[width=1.6in]{%s}
        \caption{nucleosome depletion}
        \label{fig:side:a}
        \end{minipage}
    \\begin{minipage}[t]{0.5\linewidth}
        \centering
        \includegraphics[width=1.6in]{%s}
        \caption{nucleosome fuzziness}
        \label{fig:side:b}
    \end{minipage}
\end{figure}

"""%(customTEXT,
     siteheat_filenames,
     conf_dict['General']['outname'] +"_profile.pdf",
     str(NFRscore),
     str(NFR_judge),
     str(PSarray),
     str(PSarray_judge),
     conf_dict['General']['outname'] +"_Nucdep.pdf",
     conf_dict['General']['outname'] +"_Nucfuzz.pdf")
    QCdoc += """
\\newpage
\\newpage
\subsection{Well-positioned nucleosome arrays}
\\begin{quotation}
Regions with well-positioned nucleosome arrays are detected as previous described (Zhang, et al., 2014), and the enrichment in potential regulatory regions (downstream promoter and union DNase I hypersensitive sites (DHS sites)) is listed. Enrichment was defined as observed/expected percentage of nucleosome array on promoter ( $>$ 1 for enriched). Expected percentage was equal to the percentage of promoter length compared to the total length of effective genome. Similar approach was applied on union DHS sites.  For each region with well-positioned nucleosome array, its genomic coordinates together with nucleosome profile values were reported in the output file: %s. Nucleosome arrays with fold enrichment of DHS sites less than 2 is regarded as “Fail” in this measurement,while fold enrichment of UTR regions less than 1 is regarded as "Fail", indicating the well-positioned nucleosome arrays are more likely to be caused by random rather than the barrier model.
\end{quotation}
\\begin{table}[h]
\caption{Enrichment of well-positioned nucleosome arrays}\label{bstable}
\\begin{tabularx}{\\textwidth}{ |X|X| }    
\hline
genomic region(Category) &  enrichment \\\\
\hline
downstream promoter & %s [%s] \\\\
\hline
union DHS sites &  %s [%s]\\\\
\hline
\end{tabularx}
\end{table}
\\begin{figure}[h]   
  \\begin{minipage}[t]{0.5\linewidth}   
    \centering   
    \includegraphics[width=2in]{%s}   
    \caption{enrichment on UTR}   
    \label{fig:side:a}   
  \end{minipage} 
  \\begin{minipage}[t]{0.5\linewidth}   
    \centering   
    \includegraphics[width=2in]{%s}   
    \caption{enrichment on DHS}   
    \label{fig:side:b}   
  \end{minipage}
\end{figure}
"""%(strlatexformat(conf_dict['Step3_nucarray']['profilebw_onarray'].split("/")[-1]),
     UTR_fold,
     UTR_judge,
     DHS_fold,
     DHS_judge,
     conf_dict['General']['outname'] +"_utr.pdf",
     conf_dict['General']['outname'] +"_dhs.pdf")
    QCdoc += """
\\newpage
\\newpage
\section{Output list}
\\begin{quotation}
All output files were described in the following table
\end{quotation}
\\begin{table}[h]
\caption{output list}\label{bstable}
\\begin{tabularx}{\\textwidth}{ |l|X| }
    
\hline
filename & description  \\\\
\hline
%s/%s & mapped reads on the genome  \\\\
\hline
%s & genome-wide nucleosome profile  \\\\
\hline
%s & genome-wide nucleosome dyad profile  \\\\
\hline
%s & nucleosome signal on promoter regions  \\\\
\hline
%s & nucleosome fragment length distribution  \\\\
"""%(strlatexformat(conf_dict['Step1_preprocess']['mapped_bed'].split("/")[-1]),
     strlatexformat(conf_dict['Step1_preprocess']['mapped_bb'].split("/")[-1]),
     strlatexformat(conf_dict['Step1_preprocess']['profilebw'].split("/")[-1]),
     strlatexformat(conf_dict['Step1_preprocess']['centerbw'].split("/")[-1]),
     strlatexformat(conf_dict['Step2_QC']['ave_tssprofile'].split("/")[-1]),
     strlatexformat(conf_dict['Step2_QC']['fraglen'].split("/")[-1])
     )
    if int(conf_dict['Step2_QC']['plotcustom']) == 1:
        QCdoc += """
\hline
%s & nucleosome signal on custom regions \\\\         
"""%(strlatexformat(conf_dict['Step2_QC']['custom_profile'].split("/")[-1]))
    QCdoc += """
\hline
%s & well-positioned nucleosome arrays \\\\
\hline
%s & gene level annotation of nuc-arrays \\\\
\hline
%s & nucleosome signal on well-positioned nucleosome arrays \\\\
\hline
%s & nucleosome array score signal on the genome \\\\
\hline
%s & summary QC report \\\\
\hline

\end{tabularx}
\end{table} 
\end{document} 

"""%(
     strlatexformat(conf_dict['Step3_nucarray']['arrayselect'].split("/")[-1]),
     strlatexformat(conf_dict['Step3_nucarray']['gene_array_anno'].split("/")[-1]),
     strlatexformat(conf_dict['Step3_nucarray']['profilebw_onarray'].split("/")[-1]),
     strlatexformat(conf_dict['Step3_nucarray']['positionbw_center'].split("/")[-1]),
     strlatexformat(conf_dict['General']['outname'])+"\_summary.pdf")

    latexfile = conf_dict['General']['outname'] + '_summary.tex'
    outf = open(latexfile,'w')
    outf.write(QCdoc)
    outf.close()
    cmd = "pdflatex %s"%(latexfile)
    cmd2 = 'cp %s %s'%(conf_dict['General']['outname'] + '_summary.pdf',summarydir)
    if conf_dict['General']['latex'] == 1:
        rwlog(cmd,logfile)
        rwlog(cmd,logfile)
        rwlog(cmd2,logfile)
        for files in os.listdir(result_folder):
            if os.path.isfile(files) and files[-12:-4] == "_summary":
                if not files[-4:] in ['.tex','.pdf','.png','.txt']:
                    cmd = "rm %s"%(files)
                    rwlog(cmd,logfile)
        wlog('pdflatex was detected in default PATH, generate summary report %s'%('summary/'+conf_dict['General']['outname'] + '_summary.pdf'),logfile)
    else:
        wlog('pdflatex was not detected in default PATH, generate summary report .tex file in summary/plots folder, you can move the whole summary/plots/ folder to the environment with pdflatex installed and run cmd in the plots/ folder: "pdflatex %s"'%(conf_dict['General']['outname'] + '_summary.tex'),logfile)
   
        
    if conf_dict['clean']:
        wlog('--clean pararmeter was turned on, remove preprocess, nucarray, QC and annotation folders',logfile)
        rwlog("rm -r %s"%(conf_dict['General']['outputdirectory'] + 'preprocess/'),logfile)
        rwlog("rm -r %s"%(conf_dict['General']['outputdirectory'] + 'nucarray/'),logfile)
        rwlog("rm -r %s"%(conf_dict['General']['outputdirectory'] + 'QC/'),logfile)
        rwlog("rm -r %s"%(conf_dict['General']['outputdirectory'] + 'annotation/'),logfile)
    else:
        wlog('--clean pararmeter was turned off, remove internal file only',logfile)
        rwlog("rm -r %s"%(conf_dict['General']['outputdirectory'] + 'annotation/'),logfile)
        # rwlog("rm %s "%(conf_dict['General']['outputdirectory'] + 'preprocess/centerbed_tmp.bed'),logfile)
        rwlog("rm %s "%(conf_dict['General']['outputdirectory'] + 'preprocess/extbed_tmp.bed'),logfile)
        rwlog("rm %s "%(conf_dict['General']['outputdirectory'] + 'preprocess/centerbed_tmp.bed'),logfile)
        rwlog("rm %s "%(conf_dict['General']['outputdirectory'] + 'preprocess/rawbed_sortbychrm_tmp.bed'),logfile)
        rwlog("rm %s "%(conf_dict['General']['outputdirectory'] + 'preprocess/clipsortbdg_TMPbed2bw.bdg'),logfile)
        rwlog("rm %s "%(conf_dict['General']['outputdirectory'] + 'preprocess/freglen_tmp.txt'),logfile)
        rwlog("rm %s "%(conf_dict['General']['outputdirectory'] + 'preprocess/rawbed_sortbychrm_tmp.bed'),logfile)
        rwlog("rm %s "%(conf_dict['General']['outputdirectory'] + 'preprocess/'+conf_dict['General']['outname']+'.bed.tmp'),logfile)
        rwlog("rm %s "%(conf_dict['General']['outputdirectory'] + 'preprocess/'+conf_dict['General']['outname']+'_center.wig'),logfile)
        rwlog("rm %s "%(conf_dict['General']['outputdirectory'] + 'preprocess/'+conf_dict['General']['outname']+'_GemomeLengthTmp.genome'),logfile)
        # rwlog("rm %s "%(conf_dict['General']['outputdirectory'] + 'preprocess/'+conf_dict['General']['outname']+'_minus1bp.bed'),logfile)
        # rwlog("rm %s "%(conf_dict['General']['outputdirectory'] + 'preprocess/'+conf_dict['General']['outname']+'_PEtoSE.bed'),logfile)
        rwlog("rm %s "%(conf_dict['General']['outputdirectory'] + 'preprocess/'+conf_dict['General']['outname']+'_profile.bdg'),logfile)
        # rwlog("rm %s "%(conf_dict['General']['outputdirectory'] + 'preprocess/'+conf_dict['General']['outname']+'_sdplus1bp.bed'),logfile)
    
    rwlog("rm %s "%(conf_dict['General']['outputdirectory'] + 'progress_filelist.txt'),logfile)

    wlog('Step4 summary DONE, check %s for final outputs'%(summarydir),logfile)

    return conf_dict









