#!/usr/bin/env python3

import argparse
import re
import sys

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

COORDINATE_REGEX = re.compile(r'.*(chr.*)\.(\d+)-(\d+).*')


#
# Functions
#

def GetKV(key, vals):

    lk = len(key)

    for v in vals:
        if len(v) >= len(key) and v[0:lk] == key:
            return v[lk:]

    else:
        return None


def GetStrand(value):
    if value & 16 != 0:
        return '-'
    else:
        return '+'


def CIGARToArrays(cigar):
    ops = []
    lengths = []
    i1, i2 = 0, 0
    end = len(cigar)
    opVals = re.findall(r'(\d+)([\w=])', cigar)
    lengths = [int(opVals[i][0]) for i in range(0, len(opVals))]
    ops = [opVals[i][1] for i in range(0, len(opVals))]

    return ops, lengths


class Reference:
    def __init__(self, name, length, index):
        self.name = name
        self.length = length
        self.index = index


class SAMEntry:
    def __init__(self, line):
        v = ParseSamLine(line)
        if v is None:
            self.title = None
            return None
        else:
            (self.title, self.flag, self.tName, self.tPos, self.mapqv, self.qStart, self.qEnd, self.readlen, self.seq, self.tlen) = (v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9])

        vals = line.split("\t")
        self.cigar = vals[5]
        self.ops, self.lengths = CIGARToArrays(self.cigar)
        self.strand = GetStrand(self.flag)
        self.tLen = 0

        for i in range(len(self.ops)):
            if self.ops[i] in ('M', 'D', '=', 'X'):
                self.tLen += self.lengths[i]

        prefixSoftClip = 0
        suffixSoftClip = 0

        # find the first soft clip cigar
        l = len(self.ops)
        if l > 1 and self.ops[0] == 'S':
            prefixSoftClip = self.lengths[0]

        elif l > 2 and self.ops[1] == 'S':
            prefixSoftClip = self.lengths[1]


        # The SAM alignment is in the direction of the target, so
        # the soft clipped end is the beginning of the reverse
        # strand.
        l = len(self.ops)
        if l > 2 and self.ops[-1] == 'S' and self.ops[-2] != 'S':
            suffixSoftClip = self.lengths[-1]

        elif l > 3 and self.ops[-2] == 'S' and self.ops[-3] != 'S':
            suffixSoftClip = self.lengths[-2]

        if self.strand == '-':
            tmp = prefixSoftClip
            prefixSoftClip = suffixSoftClip
            suffixSoftClip = tmp

        self.qStart += prefixSoftClip
        self.qEnd   -= suffixSoftClip
        self.tStart = self.tPos
        self.tEnd = self.tPos + self.tLen
        self.line = line
        self.fullReadLength = GetKV("XQ:i:", vals[11:])
        self.vals = vals

        if self.fullReadLength is not None:
            self.fullReadLength = int(self.fullReadLength)

    def PrintIntervals(self, out):
        out.write("{},{}\t{},{}\n".format(self.tStart, self.tEnd, self.qStart, self.qEnd))

titlei = 0
flagi = 1
tnamei = 2
tposi = 3
mapqvi = 4
qstarti = 5
qendi = 6
readleni = 7
seqi = 8
tleni = 9


def ParseSamLine(line):
    try:
        vals = line.split("\t")

        title = vals[0]
        flag = int(vals[1])
        tName = vals[2]
        tPos = int(vals[3])
        mapqv = int(vals[4])
        seq = vals[9]
        idx = 11
        readlen = len(vals[9])

        start = GetKV("XS:i:", vals[11:])

        if start is not None:
            start = int(start)
        else:
            start = 0

        end = GetKV("XE:i:", vals[11:])

        if end is not None:
            end = int(end)
        else:
            end = len(seq)

        tLen = int(vals[8])

    except:
        print("Error parsing")
        print(line)
        return None


    #       0      1     2      3     4      5      6    7,     8      9
    return (title, flag, tName, tPos, mapqv, start, end, readlen, seq, tLen)


def Overlap(a, b):

    v = (a, b)

    if a[1] < b[1]:
        i, j = 0, 1

    else:
        i, j = 1, 0

    if v[i][1] < v[j][0]:
        return 0.0

    if v[i][0] < v[j][0]:
        overlap = v[i][1] - v[j][0]
    else:
        overlap = v[i][1] - v[j][0]

    return abs(float(overlap))


def read_fai_file(fai_file_name):
    """
    Get a dictionary of an FAI file (reference index).

    :param fai_file_name: Name of the file to read.

    :return: A dictionary keyed on the sequence names with one tuple of the record values per record (all integers).
    """

    fai_record_dict = {}

    with open(fai_file_name, 'r') as fai_file:

        for line in fai_file:

            line = line.strip()

            if not line:
                continue

            fai_record = line.split()

            fai_record_dict[fai_record[0]] = [int(x) for x in fai_record[1:]]

    return fai_record_dict


def parse_region_string(region):
    """
    Separate a region, such as "chr1:1000-2000", into a tuple of "name", "pos", and "end". The coordinates will
    be stripped of commas, if present, and cast as integers.

    :param region: Region string.

    :return: Tuple of "name", "pos", and "end" or `None` if the string could not be parsed.
    """

    a = region.split(':')

    if len(a) != 2:
        return None

    b = a[1].split("-")

    if len(b) != 2:
        return None

    b[0] = b[0].replace(',', '')
    b[1] = b[1].replace(',', '')

    return a[0], int(b[0]), int(b[1])


def FormatRegion(regions):

    if type(regions) == list and len(regions) >= 3:

        chrom = regions[0]
        start = regions[1]
        end = regions[2]

        start = start.replace(',', '')
        end = end.replace(',', '')

    elif type(regions) == list and len(regions) == 1:
        val = parse_region_string(regions[0])

        if (val is None):
            return None
        else:
            (chrom, start, end) = (val[0], val[1], val[2])

    elif type(regions) == str:
        if regions.find(":") >= 0:
            val = parse_region_string(regions)
            if val is None:
                return None

            else:
                (chrom, start, end) = (val[0], int(val[1]), int(val[2]))

        else:
            vals = regions.split()
            (chrom, start, end) = (vals[0], int(vals[1]), int(vals[2]))
    else:
        return None

    start = int(start)
    end = int(end)

    return chrom, start, end


def BedToRegion(bedline):
    return str(bedline[0]) + ":" + str(bedline[1]) + "-" + str(bedline[2])


def extract_sequence(region, sequence_file, fai):
    """
    Get a sequence region from an indexed sequence file.

    :param region: A tuple of the genomic coordinates with elements: name, pos, end.
    :param sequence_file: Open file to extract sequences from.
    :param fai: Dictionary of FAI records keyed by sequence name (see `read_fai_file()`).

    :return: Sequence in `seqFile` specified by `region`.
    """

    # Check FAI for region
    if region[0] not in fai:
        raise ValueError('Sequence {} is missing in the index'.format(region[0]))

    # Get coordinates and lengths
    chr_start = int(fai[region[0]][1])
    seq_len = int(fai[region[0]][2])
    line_len = int(fai[region[0]][3])

    region_pos = int(region[1])
    region_end = int(region[2])

    # Calculate psotions
    start_line = int(region_pos / seq_len)
    start_line_pos = region_pos % seq_len

    end_line = int(region_end / seq_len)
    end_line_pos = region_end % seq_len

    start_file_pos = chr_start + start_line * line_len + start_line_pos
    end_file_pos = chr_start + end_line * line_len + end_line_pos

    # Check file positions
    if start_file_pos < 0:
        raise ValueError('Region {0}:{1}-{2} attempts to seek before 0 in the sequence file'.format(*region))

    if end_file_pos < start_file_pos:
        raise ValueError(
            'Region {0}:{1}-{2} attempts to seek before the start position of its record in the '
            'sequence file'.format(*region)
        )

    # Read sequence
    sequence_file.seek(start_file_pos)

    return sequence_file.read(end_file_pos - start_file_pos).replace('\n', '')



#
# Main
#

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
fai = read_fai_file(args.genome + '.fai')

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

        aln = SAMEntry(line)

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
        ref_len = 0  # DBGTMP

        if args.condense > 0:

            while i < len(aln.lengths):
                l = aln.lengths[i]
                op = aln.ops[i]
                j = i

                if op in {M, D, E, X}:  # DBGTMP
                    ref_len += l  # DBGTMP

                if (op == I or op == D) and \
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
                        if gapLen >= 50:  # DBGTMP
                            sys.stderr.write('Condensing: {}:{} ({}): op={}, match={}, gap={}\n'.format(aln.tName, tPos + ref_len, aln.title, op, matchLen, gapLen))  # DBGTMP

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

                    targetSeq = extract_sequence((aln.tName, tPos, tPos + oplen), genome_file, fai)

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

                    delSeq = extract_sequence([chrName, delStart, delEnd], genome_file, fai)

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
