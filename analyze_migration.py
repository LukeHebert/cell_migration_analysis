'''This script takes a .csv file as input and compares DAPI counts from 
migrated cells of various genotypes. An Kruskal-Wallis test is conducted and 
a boxplot is generated from the cell group comparison.

EXAMPLE USE CASE:
python3 analyze_migration.py tidy_DAPI_counts.csv
'''

__author__ = 'Luke S Hebert'
__copywrite__ = 'Luke S Hebert, 2022-11-29'

import sys
import pandas as pd
from scipy import stats
import scikit_posthocs as sp
import seaborn as sns
import matplotlib.pyplot as plt

def main():
	#load the user-provided data
	in_path = sys.argv[1]
	df = pd.read_csv(in_path)
	
	#get a list of unique cell groups
	genotypes = list(set(df['genotype'].to_list()))
	
	#make a cell group dictionary {cellgroup:[DAPI count list]}
	#order the groups by arithmetic mean
	dapis = {g:df[df['genotype']==g]['DAPI'].to_list() for g in genotypes}
	order = sorted(
		genotypes, 
		key=lambda x: sum(dapis[x])/len(dapis[x]), 
		reverse=True)
	
	#perform Kruskal-Wallis test on the groups' DAPI counts
	#make a dictionary {cellgroup:[DAPI count Kruskal-Wallis stats]}
	kw = stats.kruskal(*dapis.values())
	dunn = sp.posthoc_dunn(list(dapis.values()), p_adjust='bonferroni')
	print('\nKey for Dunn test results interpretation:')
	for i, val in enumerate(list(dapis.keys())):
		print(f'{i+1}\t{val}')
	print(f'\nThe Kruskal-Wallis test results:'
		f'\nH stat:\t{kw.statistic}\np-val:\t{kw.pvalue}'
		f'\nThe posthoc Dunn test results:'
		f'\np-vals:\n{dunn}')
	dunn_path = in_path.replace('.csv','_dunn.csv')
	dunn.to_csv(dunn_path,index=True)
	
	#make combo boxplot + swarmplot of the groups' DAPI counts
	#group labels on X axis (minus "DAPI")
	#display key statistical take-aways on the plot
	plt.grid(axis='y', alpha=0.25)
	sns.boxplot(
		data=[dapis[x] for x in order], 
		color='powderblue')
	sns.swarmplot(
		data=[dapis[x] for x in order], 
		color='midnightblue', 
		alpha=0.75)
	plt.xlabel('PTK2 Genotype')
	plt.ylabel('Cell Counts (DAPI)')
	plt.xticks(plt.xticks()[0], order)
	plt.xticks(rotation=90)
	plt.tight_layout()
	out_path = in_path.replace('.csv','_boxplot.png')
	plt.savefig(out_path, dpi=800)

if __name__ == '__main__':
	main()