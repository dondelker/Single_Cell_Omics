import scrublet as scr
import scipy.io
import matplotlib.font_manager
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
#import h5py
import pandas as pd
import argparse

#fm = matplotlib.font_manager.json_load(os.path.expanduser("~/.cache/matplotlib/fontlist-v310.json"))
#fm.findfont("serif", rebuild_if_missing=False)

parser = argparse.ArgumentParser()
parser.add_argument('inputdir', help='"filtered_feature_bc_matrix" directory')
parser.add_argument('outdir', help='Directory path for outputs')
parser.add_argument('--threshold', type=float, nargs='?', help='Set threshold for doublet. By default, scrublet automatically sets a threshold.')
args = parser.parse_args()

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42

os.chdir(args.outdir)

input_dir = args.inputdir
#h5 = h5py.File(input_dir + '/filtered_feature_bc_matrix.h5','r')
#counts_matrix = h5['matrix']
#genes = h5['matrix/data']
counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
#genes = np.array(scr.load_genes(input_dir + '/features.tsv.gz', delimiter='\t', column=1))
genes=pd.read_csv(input_dir + '/features.tsv.gz', compression='gzip', delimiter='\t',encoding='utf-8')

print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
print('Number of genes in gene list: {}'.format(len(genes)))


scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)

doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)
#Set threshold for doublet
n = args.threshold
scrub.call_doublets(threshold=n)

df = pd.DataFrame({
    'doublet_score': scrub.doublet_scores_obs_,
    'predicted_doublet': scrub.predicted_doublets_
})
#df.to_csv('scrublet_output_table.csv', index=False)

scrub.plot_histogram();

plt.savefig('histogram.png')

#read barcode file
barcode=pd.read_csv(input_dir + '/barcodes.tsv.gz', compression='gzip', delimiter='\t',encoding='utf-8',names=['Barcode'])
print(list(barcode.columns.values)) #file header
print(barcode.tail(35)) #last N rows
barcode_df = pd.DataFrame(barcode)

#merge scores to barcode
df2 = pd.concat([df,barcode_df],axis=1)
df2.to_csv('scrublet_output_barcode_table.csv', index=False)

def getCurrentMemoryUsage():
    ''' Memory usage in kB '''

    with open('/proc/self/status') as f:
        memusage = f.read().split('VmRSS:')[1].split('\n')[0][:-3]

    return int(memusage.strip())

print('Total memory usage:{} kB'.format(getCurrentMemoryUsage()))
