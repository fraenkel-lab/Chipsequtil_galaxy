'''
 NAME

 pieplots_macs.py

 SYNOPSIS

python pieplots_macs.py --genefile MACSoutfile_genes.txt --peakfile MACSoutfile_peaks.bed --outfile MACSdirectory/piechart.pdf


 DESCRIPTION

Peaks are assigned to the closest gene and then categorized according to their location at different genomic regions (promoter, intron, exon, or after the gene). Sites >10kb away from any gene are considered intergenic. (from Pamela)

'''

__author__='Renan Escalante'
__email__='renanec@mit.edu'

import pandas as pd
import matplotlib
matplotlib.use('pdf')
from matplotlib import pyplot as plt
matplotlib.rcParams['pdf.fonttype']=42
matplotlib.rcParams['font.size']=14
import sys
from argparse import ArgumentParser

def map_peaks(gene,peak,outfile,macsFlag):
    genefile = open(gene, 'r')
    peakfile = open(peak, 'r')

    types = {'promoter':0, 'after':0, 'intron':0, 'exon': 0}

    #read mapped gene file, store closest map for each peak
    peaks={} #{chrom:{peakStart:[dist, type]}}
    for line in genefile:
        words = line.strip().split('\t')
        #ignore first line
        if words[0] == 'knownGeneID': continue
        chrom = words[2]


        if not macsFlag:
            try:
                start = int(words[3])
                dist = abs(int(words[15]))
                maptype = words[16]
                if maptype == 'gene':
                    maptype = words[17]
            except:
                pass

        else:
            start = int(words[3])-1
            dist = abs(int(words[12]))
            maptype = words[14]
            if maptype == 'gene':
                maptype = words[15]


        if chrom not in peaks:
            #new chrom
            peaks[chrom] = {start:[dist,maptype]}
        else:
            if start in peaks[chrom].keys():
                #account for duplicate entries - choose closest gene and store type
                if dist < peaks[chrom][start][0]:
                    #closer gene
                    peaks[chrom][start] = [dist, maptype]
            else: peaks[chrom][start] = [dist, maptype]

    #count types - 1 per peak in peak file
    types = {'promoter':0, 'after':0, 'intron':0, 'exon': 0, 'inter': 0}
    totalpks = 0
    #Read peak file in bed format
    for line in peakfile:
        words = line.strip().split('\t')
        chrom = words[0]
        start = int(words[1])
        if chrom in peaks:
            if start in peaks[chrom]:
                types[peaks[chrom][start][1]] += 1
            else:
                types['inter'] += 1
        else:
            types['inter'] += 1
        totalpks += 1


    #--------------------------------------------
    #  make a square figure and axes
    #--------------------------------------------

    fig = plt.figure(figsize=(6,6))
    pie_ax = fig.add_axes((0.3,0.2,0.4,0.4))

    # The slices will be ordered and plotted counter-clockwise.
    labels = ['exon: %i'%types['exon'],'intron: %i'%types['intron'],'promoter: %i'%types['promoter'],'intergenic: %i'%types['inter'], 'after: %i'%types['after']]
    fracs = [types['exon'], types['intron'], types['promoter'], types['inter'], types['after']]

    plt.pie(fracs, labels=labels) #, autopct='%1.1f%%')

    #Export data frame with all the counts
    indexDataFrame = ['exon','intron','promoter','intergenic','after']
    df = pd.DataFrame(data=fracs, index=indexDataFrame)
    dfFileName = outfile.replace("pdf","csv")
    df.to_csv(dfFileName, sep=',')
    #plt.title('MACS peaks in %s'%(name))
    plt.figtext(.5, .1, 'Total: %i'%totalpks, ha='center')
    fig.savefig(outfile)

def main():
    usage = "usage: %prog --genefile MACSoutfile_genes.txt --peakfile MACSoutfile_peaks.bed --outfile MACSdirectory/piechart.pdf"
    parser = ArgumentParser(usage)
    parser.add_argument("--genefile", dest="genefile", help="Path to file MACS_mfold10,30_pval1e-5_genes.txt")
    parser.add_argument("--peakfile", dest="peakfile", help="Path to file MACS_mfold10,30_pval1e-5_peaks.bed")
    parser.add_argument("--outfile", dest="outfile", default="MACS_piechart.pdf", help="Path to pdf file where you want to store the piechart")
    parser.add_argument('--MACS',action='store_true',default=False,help='Set this value if you have MACS peaks')

    args=parser.parse_args()

    map_peaks(args.genefile, args.peakfile, args.outfile, args.MACS)


if __name__=='__main__':
    main()
