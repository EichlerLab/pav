#!/usr/bin/env python3

import argparse
import re
import sys
import os
import inspect

# Append smrtsvlib to path
smrtsv_base = os.path.abspath(
    os.path.dirname(os.path.dirname(
        os.path.abspath(inspect.getfile(inspect.currentframe()))
    )),
)

sys.path.append(smrtsv_base)

from smrtsvlib import tools

# Declarations
M = 'M'
X = 'X'
E = '='
I = 'I'
D = 'D'
N = 'N'
S = 'S'
H = 'H'
P = 'P'

MATCH_SET = {M, X, E}

COORDINATE_REGEX = re.compile('.*(chr.*)\.(\d+)-(\d+).*')


# Parse arguments
ap = argparse.ArgumentParser(description="Print gaps in a SAM file.")

ap.add_argument("genome",
                help="Genome file with a .fai")
ap.add_argument("sam", nargs="+",
                help="Sam file of alignment.")
ap.add_argument("--onTarget", default=False, action='store_true',
                help="Assume the query encodes the position of the aligned sequence, and make sure at least the "
                     "chromosomes match.")
ap.add_argument("--gapFree",
                help="Print sequences without gaps.", default=None)
ap.add_argument("--minContigLength", default=0, type=int,
                help="Only parse alignments from contigs this length or greater")
ap.add_argument("--minLength", default=50, type=int,
                help="Minimum gap length (inclusive).")
ap.add_argument("--maxLength", default=None, type=int,
                help="Maximum gap length (exclusive).")
ap.add_argument("--outFile", default=None,
                help="Print output here, default= stdout")
ap.add_argument("--context", default=0, type=int,
                help="Print surrounding context")
ap.add_argument("--condense", default=0, type=int,
                help="Pack indels if the matches separating them is less than this value.")
ap.add_argument("--outsam", default=None,
                help="Write the modified condensed sam to a file.")
ap.add_argument("--minq", default=10, type=int,
                help="Minimal mapping quality to consider (10)")
ap.add_argument("--snv", default=None,
                help="Print SNVs to this file.")
ap.add_argument("--nloc", default=None,
                help="Print locations of aligned N's here.")
ap.add_argument("--contigBed", default=None,
                help="Print where contigs map.")
ap.add_argument("--status", default=False, action='store_true',
                help="Print how far along the alignments are.")
ap.add_argument("--blacklist", default=None,
                help="Exclude contigs on this list from callsets.")
ap.add_argument("--removeAdjacentIndels", default=False, action='store_true',
                help="Find instances of SNVs pushed into indels, in the format: NIXMND., and remove these operations.")
ap.add_argument("--printStrand", default=False, action='store_true',
                help="Print strand of aligned contig")
ap.add_argument('--noheader', dest='header', default=True, action='store_false',
                help='Print BED header line in output.')

args = ap.parse_args()

# Read FASTA index
fai = tools.read_fai_file(args.genome + '.fai')

# Get SAM files
sam_files = list()

for file_name in args.sam:

    file_name = file_name.strip()

    # Add FOFN
    if file_name.endswith('.fofn'):

        sam_files.extend(
            [line.strip() for line in open(file_name).readlines()]
        )

    else:
        sam_files.append(file_name)

# Open output file
if args.outFile is None:
    out_file = sys.stdout
else:
    out_file = open(args.outFile, 'w')

# Open FAI file
genome_file = open(args.genome, 'r')

# Open gap-free output file
if args.gapFree is not None:
    out_gap_free_file = open(args.gapFree, 'w')
else:
    out_gap_free_file = None

# Open contig bed
if args.contigBed is not None:
    contig_bed_file = open(args.contigBed, 'w')
else:
    contig_bed_file = None

# Read blacklist
blacklist = {}

if args.blacklist is not None:

    with open(args.blacklist, 'r') as bl:

        for line in bl:
            v = line.split()

            if v[0] not in blacklist:
                blacklist[v[0]] = []

            if len(v) > 1:
                blacklist[v[0]].append(int(v[1])+1)

# Open output SAM file
if args.outsam is not None:
    out_sam_file = open(args.outsam, 'w')
else:
    out_sam_file = None

# Open SNV output file
if args.snv is not None:
    out_snv_file = open(args.snv, 'w')
else:
    out_snv_file = None

# Open N locations file
if args.nloc is not None:
    nLocOut = open(args.nloc, 'w')
else:
    nLocOut = None

# Write header
if args.header:
    out_file.write('#CHROM\tPOS\tEND\tSVTYPE\tSVLEN\tSEQ\tQUERY_ID\tQUERY_POS\tQUERY_END\tQUERY_STRAND')

    if args.context > 0:
        out_file.write('\tCONTEXT')

    out_file.write('\n')

    if out_snv_file is not None:
        out_snv_file.write('#CHROM\tPOS\tEND\tSVTYPE\tREF\tALT\tQUERY_ID\tQUERY_POS\tQUERY_STRAND\n')

# Process each SAM file
for sam_file_name in sam_files:

    sam_file = open(sam_file_name, 'r')

    line_number = 0

    for line in sam_file:
        line_number += 1

        if line[0] == '@':

            if out_sam_file is not None:
                out_sam_file.write(line)

            continue

        if len(line) <= 1:
            continue

        aln = tools.SAMEntry(line)

        if aln.title is None:
            continue

        #
        # Use 0-based coordinate system
        #

        aln.tStart -= 1
        aln.tEnd -= 1

        if args.onTarget:

            coordReMatch = COORDINATE_REGEX.match(aln.title)

            if coordReMatch is not None:

                coordMatchGroups = coordReMatch.groups()
                srcChrom = coordMatchGroups[0]
                srcStart = int(coordMatchGroups[1])
                srcEnd = int(coordMatchGroups[2])

                if srcChrom != aln.tName:
                    sys.stderr.write("off target chromosome: " + srcChrom + " " + aln.tName + "\n")
                    continue

                if not (aln.tStart <= srcStart < aln.tEnd or
                        aln.tStart <= srcEnd < aln.tEnd or
                        (srcStart < aln.tStart and srcEnd > aln.tEnd)):

                    sys.stderr.write(
                        "no overlap " + srcChrom + " " + str(srcStart) + " " + str(srcEnd) + " alignment: " +
                        str(aln.tStart) + " " + str(aln.tEnd) + "\n"
                    )

                    continue

        if aln.mapqv < args.minq:
            sys.stderr.write("low mapqv " + str(aln.mapqv) + " , skipping " + aln.title + "\n")
            continue

        if contig_bed_file is not None:
            contig_bed_file.write("{}\t{}\t{}\t{}\n".format(aln.tName, aln.tStart, aln.tStart + aln.tlen, aln.title))

        if args.minContigLength > len(aln.seq):
            sys.stderr.write("too short, skipping " + aln.title + "\n")
            continue

        if args.blacklist is not None:
            if aln.title in blacklist:

                if len(blacklist[aln.title]) == 0:
                    sys.stderr.write("Skipping " + aln.title + " in blacklist.\n")
                    continue

                else:
                    foundPos = False

                    for p in blacklist[aln.title]:
                        if int(aln.tPos) == p:
                            foundPos = True
                            break

                    if foundPos:
                        sys.stderr.write("Skipping " + aln.title + " in blacklist.\n")
                        continue

        tPos = aln.tStart
        qPos = 0

        #
        # condense matches.
        #

        packedCigar = []
        i = 0
        i1 = 1
        niter = 0
        maxGap = 0

        foundGap = False

        if args.removeAdjacentIndels:
            for i in range(1, len(aln.lengths) - 1):

                if (aln.ops[i-1] != 'M' and
                        aln.ops[i+1] != 'M' and
                        aln.ops[i-1] != aln.ops[i + 1] and
                        aln.ops[i] == 'M' and
                        aln.lengths[i - 1] == aln.lengths[i + 1] and
                        aln.lengths[i] < 4):

                    aln.lengths[i-1] = 0
                    aln.lengths[i+1] = 0

            newLengths = []
            newOps = []

            for i in range(0, len(aln.lengths)):
                if aln.lengths[i] != 0:
                    newLengths.append(aln.lengths[i])
                    newOps.append(aln.ops[i])

            aln.lengths = newLengths
            aln.ops = newOps

        packedOps = []
        packedLengths = []

        i = 0

        if args.condense > 0:

            while i < len(aln.lengths):
                l = aln.lengths[i]
                op = aln.ops[i]
                j = i

                if op == I or op == D:
                    if l > maxGap:
                        maxGap = l

                if op == I or op == D and \
                        i < len(aln.ops) - 2 and \
                        aln.ops[i+2][0] == op:

                    matchLen = 0
                    gapLen = 0

                    while j+2 < len(aln.ops) and \
                            aln.ops[j+2][0] == op and \
                            aln.ops[j+1][0] in MATCH_SET and \
                            aln.lengths[j+1] < args.condense:

                        matchLen += aln.lengths[j+1]
                        gapLen += aln.lengths[j+2]

                        j += 2

                    if j > i:
                        newIndel = (op, l+gapLen)
                        newMatch = (M, matchLen)
                        packedOps.append(op)
                        packedLengths.append(l+gapLen)

                        packedOps.append(M)
                        packedLengths.append(matchLen)

                    else:
                        packedLengths.append(l)
                        packedOps.append(op)

                else:
                    packedLengths.append(l)
                    packedOps.append(op)

                i = j + 1
                niter += 1

                if niter > len(aln.ops):
                    sys.stderr.write("ERROR! too many iterations.\n")

        else:
            packedOps = aln.ops
            packedLengths = aln.lengths

        for i in range(len(packedOps)):

            op = packedOps[i]
            oplen = packedLengths[i]

            if op == N or op == S:
                # Inside match block (if op == M)
                qPos += oplen

            if op in MATCH_SET:

                # Inside match block (if op == M)
                if out_snv_file is not None:

                    targetSeq = tools.extract_sequence((aln.tName, tPos, tPos + oplen), genome_file, fai)

                    querySeq = aln.seq[qPos:qPos+oplen]
                    nMis = 0

                    for mp in range(0,len(targetSeq)):
                        if mp >= len(querySeq) or mp >= len(targetSeq):
                            sys.stderr.write("ERROR with seq " + aln.title + "\n")
                            continue

                        if querySeq[mp].upper() != targetSeq[mp].upper() and \
                                targetSeq[mp].upper() != 'N' and \
                                querySeq[mp].upper() != 'N':

                            nMis += 1

                            out_snv_file.write(
                                "{}\t{}\t{}\tSNV\t{}\t{}\t{}\t{}\t{}\n".format(
                                    aln.tName, tPos + mp, tPos + mp + 1,
                                    targetSeq[mp], querySeq[mp],
                                    aln.title, mp + qPos, aln.strand
                                )
                            )

                        if args.nloc is not None \
                                and (targetSeq[mp].upper() == 'N' or querySeq[mp].upper() == 'N'):
                            nLocOut.write("{}\t{}\t{}\n".format(aln.tName, tPos+mp, tPos + mp + 1))

                tPos += oplen
                qPos += oplen

            if op == I:
                if oplen >= args.minLength and (args.maxLength is None or oplen < args.maxLength):

                    foundGap = True
                    chrName = aln.tName
                    gapSeq = aln.seq[qPos:qPos+oplen]

                    context = aln.seq[qPos+oplen:min(qPos+oplen+args.context, len(aln.seq))]

                    if context == "A" * len(context) or context == "T" * len(context):
                        homopolymer = "T"
                    else:
                        homopolymer = "F"

                    if len(gapSeq) == 0:
                        sys.stderr.write('ERROR, gap seq is of zero length\n')

                    nucs = ['A', 'C', 'G', 'T']
                    fracs = [float(gapSeq.count(n))/(len(gapSeq)+1) for n in nucs]

                    doPrint = True

                    for frac in fracs:
                        if frac > 0.85:
                            doPrint = False

                    if doPrint:
                        out_file.write(
                            '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                                chrName, tPos, tPos + oplen,
                                'INS', oplen, gapSeq,
                                aln.title, qPos, qPos + oplen, aln.strand
                            )
                        )

                        if args.context > 0:
                            out_file.write("\t{}".format(homopolymer))

                        out_file.write("\n")

                qPos += oplen

            if op == D:
                if oplen >= args.minLength and (args.maxLength is None or oplen < args.maxLength):
                    foundGap = True
                    chrName = aln.tName

                    if tPos > fai[chrName][0]:
                        sys.stderr.write(
                            'ERROR! tpos is past the genome end.' + str(tPos) + ' ' + str(fai[chrName][0]) + '\n'
                        )

                    delStart = max(tPos - args.context, 0)
                    delEnd = min(tPos + args.context + oplen, fai[chrName][0])

                    if delEnd < delStart:
                        continue

                    context = aln.seq[qPos+oplen:min(qPos+oplen+args.context, len(aln.seq))]

                    if context == "A"*len(context) or context == "T"*len(context):
                        homopolymer = "T"
                    else:
                        homopolymer = "F"

                    delSeq = tools.extract_sequence([chrName, delStart, delEnd], genome_file, fai)

                    out_file.write(
                        '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                            chrName, tPos, tPos + oplen,
                            'DEL', oplen, delSeq,
                            aln.title, qPos, qPos + 1, aln.strand
                        )
                    )

                    if args.context > 0:
                        out_file.write("\t{}".format(homopolymer))

                    out_file.write("\n")

                tPos += oplen

            if op == H:
                pass

        if out_gap_free_file is not None and not foundGap:
            out_gap_free_file.write(aln.tName + "\t" + str(aln.tStart) + "\t" + str(aln.tEnd) + "\t" + aln.title + "\n")

        if out_sam_file is not None:
            packedCigar = ''.join([str(v[0]) + str(v[1]) for v in zip(packedLengths, packedOps)])
            vals = line.split()
            packedLine = '\t'.join(vals[0:5]) + "\t" + packedCigar + '\t'.join(vals[6:]) + "\n"
            out_sam_file.write(packedLine)

# Close files
out_file.close()
genome_file.close()

if out_gap_free_file is not None:
    out_gap_free_file.close()
    out_gap_free_file = None

if out_sam_file is not None:
    out_sam_file.close()
    out_sam_file = None

if contig_bed_file is not None:
    contig_bed_file.close()
    contig_bed_file = None
