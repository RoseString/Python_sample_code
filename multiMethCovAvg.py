#!/nv/hp10/dsun33/data2/anaconda/envs/py37/bin/python
import sys
import argparse
import textwrap
import pandas as pd
from functools import reduce
import numpy as np
import math
import multiprocess
from multiprocess import Pool

### Usage
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
     	description=textwrap.dedent('''\
		-----------------------------------------------------------------------------------------------------------------------------------------
		This script will take in multiple Bismark coverage files and calculate the mean methylation levels per site across samples in multi-threads.
		For a site in a sample to be taken into calculation, the site must have at least >= cutoff reads in that sample.
		
		The average for a site will be calculated or not, depending on the following three modes:
		MODE 1) All samples have enough reads at that site
		MODE 2) At least half of the samples have enough reads at that site
		MODE 3) At least 3 samples have enough reads at that site

		The output will have the following format: Scaffold    Position    Average_Methylation_Level.

		Author: Dan Sun
		Created: 01/11/2019
		Last modified: 02/18/2019
		------------------------------------------------------------------------------------------------------------------------------------------
        '''))
parser.add_argument('--covs', nargs='+', help='Multiple bismark coverage files (at least two); no need to be sorted', required=True)
parser.add_argument('--mode', default = 1, type=int, help='1) All 2) At least half 3) At least N samples (specify with --N) should have enough coverage for the average of a site to be computed', required=False)
parser.add_argument('--N', default = 3, type=int, help='if mode == 3, specify the minimum number of samples for calculation (default is three samples)', required=False)
parser.add_argument('--cutoff', default = 5, type=int, help='depth cutoff', required=False)
parser.add_argument('--threads', default = 4, type=int, help='Number of threads', required=False)
parser.add_argument('--out', help='output name', required=True)
args = parser.parse_args()

def parallelize_dataframe(df, func):
    df_split = np.array_split(df, num_cores*2)
    pool = Pool(num_cores)
    df = pd.concat(pool.map(func, df_split))
    pool.close()
    pool.join()
    return df

def calcMean(row):
	nGood = 0 # Number of samples with enough coverage
	mls = []
	for c in range(1, nSamples+1): # ML: 3n+1; M: 3n+2; U: 3n+3
		depth = row[3*c+1] + row[3*c+2]
		if depth >= args.cutoff:
			nGood += 1
			mls.append(row[3*c])
	if args.mode == 1:
		if nGood == nSamples:
			row["Avg_ML"] = round(np.mean(mls) * 100, 6)
	elif args.mode == 2:
		if nGood >= math.ceil(1/2 * nSamples):
			row["Avg_ML"] = round(np.mean(mls) * 100, 6)
	elif args.mode == 3:
		if nGood >= args.N:
			row["Avg_ML"] = round(np.mean(mls) * 100, 6)
	return row

def applyCalcMean(df):
	df = df.apply(calcMean, axis=1)
	return df

if __name__ == "__main__":
	num_cores = args.threads
	print("%d threads will be used." % num_cores)
	print("Reading and merging bismark coverage files...", file=sys.stderr)
	dfs = []
	nSamples = len(args.covs)
	for cov in args.covs:
		df = pd.read_csv(cov, header=None, names = ['Scaffold', 'Start', 'End', 'ML', 'M', 'U'], sep='\t')
		dfs.append(df)

	df_final = reduce(lambda left,right: pd.merge(left,right,on=['Scaffold','Start','End']), dfs)
	df_final["Avg_ML"] = np.nan

	print("Calculating average methylation level per site across samples...", file=sys.stderr)

	df_final = parallelize_dataframe(df_final, applyCalcMean) # multi-threaded (better performance)
	#df_final = applyCalcMean(df_final) # single-threaded
	df_final = df_final[df_final.Avg_ML.notnull()] # Remove Cytosine sites without enough coverage
	df_final = df_final.loc[:,['Scaffold','Start','End','Avg_ML']]

	df_final.to_csv(args.out, index = None, header = None, sep="\t")

	print("Done!", file=sys.stderr)
