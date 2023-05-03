import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from Bio import SeqIO
from Bio.SeqUtils import GC

def plot_fastq_qualities(filename, ax=None, limit=10000):
    fastq_parser = SeqIO.parse(open(filename, "rt"), "fastq")
    res = []
    c = 0
    for record in fastq_parser:
        score = record.letter_annotations["phred_quality"]
        res.append(score)
        c += 1
        if c > limit:
            break
    df = pd.DataFrame(res)
    l = len(df.T) + 1

    if ax == None:
        f, ax = plt.subplots(figsize=(12, 5))
    rect = patches.Rectangle((0, 0), l, 20, linewidth=0, facecolor='r', alpha=.4)
    ax.add_patch(rect)
    rect = patches.Rectangle((0, 20), l, 8, linewidth=0, facecolor='yellow', alpha=.4)
    ax.add_patch(rect)
    rect = patches.Rectangle((0, 28), l, 12, linewidth=0, facecolor='g', alpha=.4)
    ax.add_patch(rect)
    df.mean().plot(ax=ax, c='black')
    boxprops = dict(linestyle='-', linewidth=1, color='black')
    df.plot(kind='box', ax=ax, grid=False, showfliers=False, color=dict(boxes='black', whiskers='black'))
    ax.set_xticks(np.arange(0, l, 5))
    ax.set_xticklabels(np.arange(0, l, 5))
    ax.set_xlabel('position(bp)')
    ax.set_xlim((0, l))
    ax.set_ylim((0, 40))
    ax.set_title('per base sequence quality')
    return

def fastq_to_dataframe(filename, size=1000):
    ext = os.path.splitext(filename)[1]
    if ext == '.fastq' or ext == '.gz':
        fastq_parser = SeqIO.parse(open(filename, "rt"), "fastq")
    else:
        fastq_parser = SeqIO.parse(open(filename, "r"), "fastq")
    i = 0
    res = []
    for fastq_rec in fastq_parser:
        i += 1
        if i > size:
            break
        res.append([fastq_rec.id, str(fastq_rec.seq)])
    df = pd.DataFrame(res, columns=['id', 'seq'])
    df['length'] = df.seq.str.len()
    return df

def normpdf(x, mean, sd):
    var = float(sd) ** 2
    denom = (2 * math.pi * var) ** .5
    num = math.exp(-(float(x) - float(mean)) ** 2 / (2 * var))
    return num / denom

def plot_fastq_gc_content(filename, ax=None, limit=50000):
    if ax == None:
        f, ax = plt.subplots(figsize=(12, 5))
    df = fastq_to_dataframe(filename, size=limit)
    gc = df.seq.apply(lambda x: GC(x))
    gc.hist(ax=ax, bins=150, color='black', grid=False, histtype='step', lw=2)
    ax.set_xlim((0, 100))
    x = np.arange(1, 100)
    f = [normpdf(i, gc.mean(), gc.std()) for i in x]
    ax2 = ax.twinx()
    ax2.plot(x, f)
    ax2.set_ylim(0, max(f))
    ax.set_title('GC content', size=15)
    return

# Set the directory containing the FASTQ files
fastq_dir = 'C:/Users/00090473/PycharmProjects/pythonProject'

# Loop over each file in the directory
for filename in os.listdir(fastq_dir):
    if filename.endswith('.fastq'):
        filepath = os.path.join(fastq_dir, filename)

        # Perform QC analysis
        fig, (ax1, ax2) = plt.subplots(2, figsize=(12, 10))
        plot_fastq_qualities(filepath, ax=ax1, limit=100000)
        plot_fastq_gc_content(filepath, ax=ax2)

        # Save the QC analysis figures
        output_filename = os.path.splitext(filename)[0] + '_QC.png'
        output_filepath = os.path.join(fastq_dir, output_filename)
        plt.savefig(output_filepath)
        plt.close(fig)

        print(f"QC analysis completed for '{filename}'. Output saved as '{output_filename}'.")

