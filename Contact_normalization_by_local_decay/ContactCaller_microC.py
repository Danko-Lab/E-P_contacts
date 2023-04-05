import sys
import gzip
import os
import string
import datetime
import argparse
import pdb

import joblib
import numpy as np
import pandas as pd
from pandas import Series, DataFrame
import statsmodels.api as sm
import scipy
import scipy.stats as stats
import rpy2
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import IntVector as ivect
import rpy2.robjects as robjects
from joblib import Parallel, delayed
from fast_histogram import histogram1d

def parse_arguments(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument('--peakfile', type=str,
                        help='File of "bait" peaks to analyze (format: Chromosome <\t> bait-center-position')
    parser.add_argument('--preyfile' , type=str,
                        help='File of "prey" PREYs to analyze (format: Chromosome <\t> prey-center-position')
    parser.add_argument('--confile', type=str,
                        help=' Contact file in juicer short format (given from merged_nodup.txt')
    parser.add_argument('--out', type=str,
                        help='Path to output file')
    parser.add_argument('--dist', type=int,
                        help='Half size of "prey" search window around "bait" position',default=1000000)
    parser.add_argument('--cap', type=int,
                        help='Half size of "bait" and "prey" contact capture window',default=2000 )
    parser.add_argument('--index', type=int,
                        help='')
    parser.add_argument('--chrom', type=str,
                        help='')
    parser.add_argument('--position', type=int,
                        help='')

    return parser.parse_args(argv)

def optimize_lowess(contact_counts_zero, winsize, delta, dist=1000 * 1000):
    zeroModel = []
    lowess_sm = sm.nonparametric.lowess

    for i in range(dist // winsize):
        start, stop = i * winsize, (i + 1) * winsize
        pos = np.arange(1, (winsize + 1), 1)
        zeroPDF = np.zeros((stop - start), dtype=float)

        for k in range(start, stop):
            zeros = np.sum(contact_counts_zero[k - 50:k + 50])
            zeroPDF[k - start] = float(zeros / 100) if (k >= 50 and (k + 50) <= (dist - 1)) else 0.0

        zeroSM = lowess_sm(zeroPDF, pos, frac=0.01, it=3, delta=delta, return_sorted=False)
        zeroModel0 = [0.5 * (1 - i) for i in zeroSM]
        zeroModel.extend(zeroModel0)

    zeroModel = np.asarray(zeroModel)
    return zeroModel


def optimize_lowess2(x, x0, pos, winsize, delta,distance, dist=1000 * 1000 ):
    short_pos = pos[0:1000]
    pseudo_counts = x[0:1000]
    pseudo_counts = np.asarray(pseudo_counts)
    lowess_sm = sm.nonparametric.lowess
    bgModel = lowess_sm(pseudo_counts, short_pos, frac=0.05, it=3, delta=0.0, return_sorted=False)

    for i in range(dist // winsize):
        start, stop = i * winsize, (i + 1) * winsize
        short_pos = pos[start:stop + 300] if stop < dist - 300 else pos[start:stop]
        short_counts = x[start:stop + 300] if stop < dist - 300 else x[start:stop]
        zero_counts = x0[start:stop + 300] if stop < dist - 300 else x0[start:stop]
        pseudo_counts = np.add(short_counts, zero_counts)
        short_pos = np.asarray(short_pos)
        pseudo_counts = np.asarray(pseudo_counts)
        counts_sm = lowess_sm(pseudo_counts, short_pos, frac=0.01, it=3, delta=delta, return_sorted=False)

        if i < 1:
            tail = bgModel[-(300):]
            head = counts_sm[(700):(1000 + 1)]
            merge = [(i + j) / 2 for (i, j) in zip(tail, head)]
            tmpModel, counts = bgModel[:-(300)], counts_sm[(1000 + 1):]
        else:
            tail = bgModel[-(300):]
            head = counts_sm[:(300)]
            merge = [(i + j) / 2 for (i, j) in zip(tail, head)]
            tmpModel, counts = bgModel[:-(300)], counts_sm[(300):]

        bgModel = []
        bgModel.extend(tmpModel)
        bgModel.extend(merge)
        bgModel.extend(counts)

    reads = sum(val < dist for val in distance)
    print(reads)
    bgPDF = [float(i / reads) for i in bgModel]
    return bgPDF

def fisher_test(a1, a2, b1, b2, method=0):
    print(a1, a2, b1, b2)
    try:
        if method == 0:    # 0 fisher_test by R ;1 fisher_test by python
            v = robjects.FloatVector([a1, a2, b1, b2])
            robjects.r.assign("v", v)
            test = robjects.r("matrix(v, nrow = 2)")
            robjects.r.assign("matrix", test)

            robjects.r("exact <- fisher.test(matrix, alternative = 'greater')")
            result = robjects.r("exact")
            p_val = result[0][0]
            print(p_val)
        else:
            result = stats.fisher_exact([[a1,a2],[b1,b2]], alternative='greater')
            p_val = result[1]
            print(p_val)
    except:
        p_val = 999.99
    return p_val


def main(args):
    args.chrom = args.chrom.strip("chr")

    print("reading preyfile:", args.preyfile )
    plus, minus, contactProb = [], [], []
    preyFile_DF = pd.read_csv(args.preyfile, names=['args.chrom', 'pos'], sep='\t')

    # read the contacts around loci from the splited temp_contact_file
    print("reading confile:",args.confile)
    try:
        #contact_fp = gzip.open(args.confile + "_IN_" + args.peakFile + "_folder/temp_" + str(args.index) + ".gz")
        contact_fp= gzip.open(args.confile)
    except Exception as e:
        print(e)

    contacts = np.array([l.strip().split() for l in contact_fp.readlines()])

    print("contacts with locus number " + str(args.index))
    contacts = contacts[:,1:4]
    distance = contacts[:,2].astype(int) - contacts[:, 1].astype(int)
    max_distance = args.dist * 2 if max(distance) > args.dist * 2 else max(distance)

    print( "making histogram")
    contact_counts = histogram1d(distance, range=[0, max_distance], bins=(max_distance))
    a = np.asarray(contact_counts)
    contact_counts_zero = np.where(a != 0, 0, 1)
    pos = np.arange(1, (max(distance)) + 1, 1)

    print("calculating lowess")
    zeroModel = optimize_lowess(contact_counts_zero, 5000, 16, args.dist )
    bgPDF = optimize_lowess2(contact_counts, zeroModel, pos, 5000, 16, distance, args.dist)

    print("calculating fisher_test")
    baitStart, baitStop = int(args.position - args.cap), int(args.position + args.cap)
    ## Get interactions around bait
    for contact in contacts:
        if baitStart <= int(contact[1]) <= baitStop:
            plus.append(int(contact[2]))
        if baitStart <= int(contact[2]) <= baitStop:
            minus.append(int(contact[1]))

    minDist = 5000
    strand=1
    contact_in_baits=[]
    for i in preyFile_DF.index:
        if preyFile_DF['pos'][i] - args.position > minDist:
            strand=1
            contact_in_baits=plus
        if preyFile_DF['pos'][i] - args.position < (-1 * minDist):
            strand=-1
            contact_in_baits=minus

        preyPosition, preyDist = preyFile_DF['pos'][i], preyFile_DF['pos'][i] - args.position
        preyStart, preyStop = (preyPosition - (args.cap)), (preyPosition + (args.cap))
        absDist = (preyDist * strand)
        if (absDist - args.cap) >= 0:
            expStart, expStop = int(absDist - args.cap), int(absDist + args.cap)
        else:
            expStart, expStop, mirror = 0, int(absDist + args.cap), int(-1 * (absDist - args.cap))
        exp_prob = sum(bgPDF[expStart:expStop + 1])
        expected = (len(contact_in_baits) * exp_prob)
        observed = len([x for x in contact_in_baits if preyStart <= x <= preyStop])

        p_val = fisher_test(observed, expected, (len(a) - observed), (len(a) - expected), 0)
        contactProb.append([args.chrom, args.position, preyPosition, preyDist, p_val, observed, expected,
                                     (len(a) - observed), (len(a) - expected)])



    out = open(args.out, "w")
    print("output to",args.out)

    p_count = 0
    for prob in contactProb:
        # output file format - chr <\t> bait position <\t> prey position <\t> p_val <\t> observed <\t> expected <\t> FDR
        outStr = 'chr' + str(prob[0]) + "\t" + str(prob[1]) + "\t" + str(prob[2])+"\t" + str(prob[3]) + "\t" + str(prob[4]) + "\t" + str(
            prob[5]) + "\t" + str(prob[6]) + "\t" + str(prob[7]) + "\t" + str(prob[8]) +  "\n"
        out.write(outStr)
        p_count += 1


if __name__ == '__main__':
    args = parse_arguments(sys.argv[1:])
    print()
    main(args)


