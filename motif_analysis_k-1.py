#!/usr/bin/python

import re, os, sys, getopt

def calc_score(motif, a, c, g, t, glen, obs):
    p = 1
    bpp = {'A': a, 'C': c, 'G': g, 'T': t, 'N': 0, 'M': 0, 'R': 0, 'W': 0, 'S': 0, 'K': 0, 'Y': 0, 'B': 0, 'D': 0, 'H': 0, 'V': 0}
    mlen = len(motif)
    bps = list(motif)
    for i in range(mlen):
        bp = bps[i]
        pbp = bpp[bp]
        p = p * pbp
    exp = glen * p
    score = float((obs - exp))/(obs + exp)
    return score, exp

def calc_score_2(motif, n, glen):
    mlen = len(motif)
    mlen1 = mlen - 1
    mlen2 = mlen - 2
    m1a = motif[0:mlen1]
    m1b = motif[1:mlen]
    m2 = motif[1:mlen1]
    #print(m1a+" "+m1b+" "+m2)
    q = 0
    r = 0
    s = 1
    if m1a in n.keys():
         q = float(n[m1a])
    if m1b in n.keys():
        r = float(n[m1b])
    if m2 in n.keys():
        s = float(n[m2])
    exp = q * r / s # glen *
    score = float((n[motif] - exp)/(n[motif] + exp))
    return score, exp

def main(argv):
    inputfile = ''
    outputfile = ''
    motif = ''
    pos = 0
    n_ = 0
    n = dict()
    n1 = dict()
    n2 = dict()
    bplist = list()
    acgt = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0, 'M': 0, 'R': 0, 'W': 0, 'S': 0, 'K': 0, 'Y': 0, 'B': 0, 'D': 0, 'H': 0, 'V': 0}
    genome_len = 0
    nchr = 0
    spec = ''

    try:
        opts, args = getopt.getopt(argv,"hi:o:s:n:",["ifile=","ofile=","spec=","n="])
    except getopt.GetoptError:
        print ('motif_analysis.py -i <inputfile> -o <outputfile> -s <species name> -n <motif length>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('motif_analysis.py -i <inputfile> -o <outputfile> -s <species name> -n <motif length>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-n"):
            n_ = int(arg)
        elif opt in ("-s"):
            spec = arg

    fi = open(inputfile,'r')
    fo = open(outputfile,'w')

    for line in fi:
        lin = line.rstrip('\n')
        line = str(lin)
        x = re.search(">",line)
        if x :
            pos = 0; motif = ""; bplist = list(); nchr = nchr + 1
        else:
            w = list(line)
            for b in w:
                acgt[b.upper()] = acgt[b.upper()] + 1
                pos = pos+1
                genome_len = genome_len + 1

                bplist.append(b)
                if pos >= n_ + 1:
                    bplist.pop(0)

                if pos >= n_:
                    motif = ''.join(bplist).upper()
                    x1 = n_ - 1
                    x2 = n_ - 2
                    motif1 = motif[0:x1]
                    motif2 = motif[0:x2]
                    z = re.search("[NMRWSYKBDHV]",motif)
                    if not z:
                        if motif in n:
                            n[motif] = n[motif] + 1
                        else:
                            n[motif] = 1
                        if motif1 in n:
                            n[motif1] = n[motif1] + 1
                        else:
                            n[motif1] = 1
                        if motif2 in n:
                            n[motif2] = n[motif2] + 1
                        else:
                            n[motif2] = 1

    tacgt = acgt['A'] + acgt['C'] + acgt['G'] + acgt['T']
    ap = float(acgt['A'])/tacgt
    cp = float(acgt['C'])/tacgt
    gp = float(acgt['G'])/tacgt
    tp = float(acgt['T'])/tacgt

    fo.write("#Species\tNo. chr.\tGenome length\tA%\tC%\tG%\tT%\n")
    fo.write("#"+spec+"\t"+str(nchr)+"\t"+str(genome_len)+"\t"+str(ap)+"\t"+str(cp)+"\t"+str(gp)+"\t"+str(tp)+"\n")
    fo.write("#Motif\tObserved\tExpected\tScore\n")

    for motif in sorted(n):
        if len(motif) == n_:
            #score, exp = calc_score(motif, ap, cp, gp, tp, genome_len, n[motif])
            score, exp = calc_score_2(motif, n, genome_len)
            fo.write(motif+"\t"+str(n[motif])+"\t"+str(exp)+"\t"+str(score)+"\n")

    fo.close()
    fi.close()

if __name__ == "__main__":
    main(sys.argv[1:])
