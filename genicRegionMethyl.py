#!/nv/hp10/dsun33/data2/anaconda/bin/python
import sys
import argparse
import textwrap
import re
import math

### Usage
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
     	description=textwrap.dedent('''\
		-----------------------------------------------------------------------------------------------------------------------------------------
		This program will calculate methylation level across different genic regions (average per gene/average genome-wide) for a given pertage window size.
		It requires an output from bed2coord.py and intersectBed that contains overlapped CpGs for each gene (see WGBS onenote 2.5 for how to get the output).
		Author: Dan Sun
		Created: 1/15/2018
		Last modified: 1/17/2018
		------------------------------------------------------------------------------------------------------------------------------------------
        '''))
parser.add_argument('-f', help='a file containing methylation levels for CpGs that overlap genic regions', required=True)
parser.add_argument('--wSize', type = int, default = 10, help='percentage window (divided a region into 100/x parts); mean is calculated per window', required=False)
parser.add_argument('--outA', help='Average methylation per gene', required=True)
parser.add_argument('--outB', help='Average methylation genome-wide', required=True)
args = parser.parse_args()

### This function is for a simple region
def simpleMeth(strand,start,end,cpgCoords,cpgMeths):
	bin_num = int(100 / args.wSize) # number of bins
	bin_meths = []
	if start == None: ## if the region does not exist
		for i in range(0, bin_num):
			bin_meths.append(None)
		return bin_meths
	region_len = end - start + 1
	if region_len < 100: ## if the region is too short for calculation (shorter than 100bp)
		for i in range(0, bin_num):
			bin_meths.append(None)
		return bin_meths
	bin_len = args.wSize * region_len / 100 # length of bins

	for i in range(0, bin_num):
		bin_start = math.floor(start + bin_len * i)
		bin_end = math.floor(start + bin_len * (i + 1) - 1)
		#print(bin_len, bin_start, bin_end)
		bin_meth_sum = 0
		bin_cpg_num = 0
		## Check if a CpG is inside this bin
		for j in range(0, len(cpgCoords)):
			if cpgCoords[j] >= bin_start and cpgCoords[j] <= bin_end:
				bin_cpg_num += 1
				bin_meth_sum += cpgMeths[j]
		## Calculate average methylation levels for each bin
		if bin_cpg_num > 0:
			bin_meths.append(float(bin_meth_sum/bin_cpg_num))
		else:
			bin_meths.append(None)
	if strand == "-":
		bin_meths.reverse() # for '-' strand, reverse the list (now 5'->3' direction)
	return bin_meths

### This function is for a complex region that is not continuous (introns and exons)
def complexMeth(strand,starts,ends,cpgCoords,cpgMeths):
	bin_num = int(100 / args.wSize) # number of bins
	bin_meths = []
	if starts == None: ## if the region does not exist
		for i in range(0, bin_num):
			bin_meths.append(None)
		return bin_meths
	region_len = 0
	loc_rel = 0
	loc_map_a = {} # relative position -> absolute position
	loc_map_r = {} # absolute position -> relative position
	for idx in range(0, len(starts)):
		region_len += ends[idx] - starts[idx] + 1
		for loc_abs in range(starts[idx], ends[idx]+1):
			loc_map_a[loc_rel] = loc_abs
			loc_map_r[loc_abs] = loc_rel
			loc_rel += 1
	if region_len < 100: ## if the region is too short for calculation (shorter than 100bp)
		for i in range(0, bin_num):
			bin_meths.append(None)
		return bin_meths
	bin_len = args.wSize * region_len / 100 # length of bins
	for i in range(0, bin_num):
		bin_start_r = math.floor(bin_len * i)
		bin_end_r = math.floor(bin_len * (i + 1) - 1)
		bin_start_a = loc_map_a[bin_start_r]
		bin_end_a = loc_map_a[bin_end_r]
		bin_meth_sum = 0
		bin_cpg_num = 0
		for j in range(0, len(cpgCoords)):
			if cpgCoords[j] >= bin_start_a and cpgCoords[j] <= bin_end_a: ## Check if a CpG is within a large region (including gaps)
				if cpgCoords[j] in loc_map_r: ## If so, make sure the CpG is not within the gap
					bin_cpg_num += 1
					bin_meth_sum += cpgMeths[j]
		## Calculate average methylation levels for each bin
		if bin_cpg_num > 0:
			bin_meths.append(float(bin_meth_sum/bin_cpg_num))
		else:
			bin_meths.append(None)
	if strand == "-":
		bin_meths.reverse() # for '-' strand, reverse the list (now 5'->3' direction)
	return bin_meths

def addToSum(aList, sumList, numList):
	for i in range(0, len(aList)):
		if aList[i] is not None:
			sumList[i] += aList[i]
			numList[i] += 1
	return (sumList, numList)

### Declare a bunch of dictionaries for later use
strand = {}
upstream_start = {}
upstream_end = {}
fiveUTR_start = {}
fiveUTR_end = {}
cds_starts = {}
cds_ends = {}
intron_starts = {}
intron_ends = {}
threeUTR_start = {}
threeUTR_end = {}
downstream_start = {}
downstream_end = {}
cpg_coord = {}
cpg_methyl = {}
### Reading coordinates and methylation data into dictionaries
with open(args.f) as fHandle:
	for line in fHandle:
		line = line.strip("\n")
		cols = line.split()
		acc = cols[3]
		if not acc in upstream_start:
			strand[acc] = cols[4]
			upstream_start[acc] = int(cols[5]) if cols[5] != "NA" else None
			upstream_end[acc] = int(cols[6]) if cols[6] != "NA" else None
			fiveUTR_start[acc] = int(cols[7]) if cols[7] != "NA" else None
			fiveUTR_end[acc] = int(cols[8]) if cols[8] != "NA" else None
			cds_starts[acc] = list(map(int, cols[9].split(',')))
			cds_ends[acc] = list(map(int, cols[10].split(',')))
			intron_starts[acc] = list(map(int, cols[11].split(','))) if cols[11] != "NA" else None
			intron_ends[acc] = list(map(int, cols[12].split(','))) if cols[12] != "NA" else None
			threeUTR_start[acc] = int(cols[13]) if cols[13] != "NA" else None
			threeUTR_end[acc] = int(cols[14]) if cols[14] != "NA" else None
			downstream_start[acc] = int(cols[15]) if cols[15] != "NA" else None
			downstream_end[acc] = int(cols[16]) if cols[16] != "NA" else None
			cpg_coord[acc] = []
			cpg_methyl[acc] = []
		cpg_coord[acc].append(int(cols[19]))
		cpg_methyl[acc].append(float(cols[20]))

### For each transcript, calculate average methylation levels per window for each genic region
upstream_meths_sum = []
fiveUTR_meths_sum = []
cds_meths_sum = []
intron_meths_sum = []
threeUTR_meths_sum = []
downstream_meths_sum = []
upstream_meths_num = []
fiveUTR_meths_num = []
cds_meths_num = []
intron_meths_num = []
threeUTR_meths_num = []
downstream_meths_num = []
for i in range(0, int(100 / args.wSize)):
	upstream_meths_sum.append(0)
	fiveUTR_meths_sum.append(0)
	cds_meths_sum.append(0)
	intron_meths_sum.append(0)
	threeUTR_meths_sum.append(0)
	downstream_meths_sum.append(0)
	upstream_meths_num.append(0)
	fiveUTR_meths_num.append(0)
	cds_meths_num.append(0)
	intron_meths_num.append(0)
	threeUTR_meths_num.append(0)
	downstream_meths_num.append(0)
outAHandle = open(args.outA,'w')
for acc in upstream_start:
	## 1. Upstream of TSS
	upstream_meths = simpleMeth(strand[acc], upstream_start[acc], upstream_end[acc], cpg_coord[acc], cpg_methyl[acc])
	## 2. 5'UTR
	fiveUTR_meths = simpleMeth(strand[acc], fiveUTR_start[acc], fiveUTR_end[acc], cpg_coord[acc], cpg_methyl[acc])
	## 3. CDS
	cds_meths = complexMeth(strand[acc], cds_starts[acc], cds_ends[acc], cpg_coord[acc], cpg_methyl[acc])
	## 4. Intron
	intron_meths = complexMeth(strand[acc], intron_starts[acc], intron_ends[acc], cpg_coord[acc], cpg_methyl[acc])
	## 5. 3'UTR
	threeUTR_meths = simpleMeth(strand[acc], threeUTR_start[acc], threeUTR_end[acc], cpg_coord[acc], cpg_methyl[acc])
	## 6.Downstream of TES
	downstream_meths = simpleMeth(strand[acc], downstream_start[acc], downstream_end[acc], cpg_coord[acc], cpg_methyl[acc])
	## Per gene pattern
	upstream_meths_str = ("\t".join(str(e) for e in upstream_meths)).replace("None","NA")
	fiveUTR_meths_str = ("\t".join(str(e) for e in fiveUTR_meths)).replace("None","NA")
	cds_meths_str = ("\t".join(str(e) for e in cds_meths)).replace("None","NA")
	intron_meths_str = ("\t".join(str(e) for e in intron_meths)).replace("None","NA")
	threeUTR_meths_str = ("\t".join(str(e) for e in threeUTR_meths)).replace("None","NA")
	downstream_meths_str = ("\t".join(str(e) for e in downstream_meths)).replace("None","NA")
	## Print per gene pattern
	outAHandle.write(acc + "\t" + upstream_meths_str + "\t" + fiveUTR_meths_str + "\t" + cds_meths_str \
				 + "\t" + intron_meths_str + "\t" + threeUTR_meths_str + "\t" + downstream_meths_str + "\n")
	## Genome-wide pattern 
	for i in range(0, int(100 / args.wSize)):
		(upstream_meths_sum, upstream_meths_num) = addToSum(upstream_meths, upstream_meths_sum, upstream_meths_num)
		(fiveUTR_meths_sum, fiveUTR_meths_num) = addToSum(fiveUTR_meths, fiveUTR_meths_sum, fiveUTR_meths_num)
		(cds_meths_sum, cds_meths_num) = addToSum(cds_meths, cds_meths_sum, cds_meths_num)
		(intron_meths_sum, intron_meths_num) = addToSum(intron_meths, intron_meths_sum, intron_meths_num)
		(threeUTR_meths_sum, threeUTR_meths_num) = addToSum(threeUTR_meths, threeUTR_meths_sum, threeUTR_meths_num)
		(downstream_meths_sum, downstream_meths_num) = addToSum(downstream_meths, downstream_meths_sum, downstream_meths_num)
outAHandle.close()

### Print genome-wide pattern
upstream_meths_avg = []
fiveUTR_meths_avg = []
cds_meths_avg = []
intron_meths_avg = []
threeUTR_meths_avg = []
downstream_meths_avg = []
for i in range(0, int(100 / args.wSize)):
	upstream_meths_avg.append(str(round(upstream_meths_sum[i] / upstream_meths_num[i], 4)))
	fiveUTR_meths_avg.append(str(round(fiveUTR_meths_sum[i] / fiveUTR_meths_num[i], 4)))
	cds_meths_avg.append(str(round(cds_meths_sum[i] / cds_meths_num[i], 4)))
	intron_meths_avg.append(str(round(intron_meths_sum[i] / intron_meths_num[i], 4)))
	threeUTR_meths_avg.append(str(round(threeUTR_meths_sum[i] / threeUTR_meths_num[i], 4)))
	downstream_meths_avg.append(str(round(downstream_meths_sum[i] / downstream_meths_num[i], 4)))

outBHandle = open(args.outB,'w')
outBHandle.write("\t".join(upstream_meths_avg) + "\t" + "\t".join(fiveUTR_meths_avg) + "\t" + "\t".join(cds_meths_avg) + "\t" +
				"\t".join(intron_meths_avg) + "\t" + "\t".join(threeUTR_meths_avg) + "\t" + "\t".join(downstream_meths_avg) + "\n")
outBHandle.close()
