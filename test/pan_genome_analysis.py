import argparse
import numpy as np
import os

# Setup argument parsing
parser = argparse.ArgumentParser(description="Pan-genome analysis script.")
parser.add_argument('-f', '--input_file', required=True, help='Path to the input gene presence absence Rtab file')
parser.add_argument('-i', '--iterations', type=int, required=True, help='Number of iterations for analysis', default=10)
parser.add_argument('-o', '--output_dir', required=True, help='Directory path for output files')
args = parser.parse_args()

# Extract command-line arguments
input_file = args.input_file
iterations = args.iterations
output_dir = args.output_dir

# Create the output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

import random
import sys
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from scipy.optimize import curve_fit

def heaps_law(x, k, alpha):
    return k * x ** alpha

def exponential_decay(x, a, b):
    return a * np.exp(-b * x)

print("Loading data!")

# Read and process input data
r = []

for h, i in enumerate(open(input_file)):
    tmp = i.strip().split('\t')[1:]
    if h == 0:
        continue
    tmp = [int(x) for x in tmp]
    r.append(tmp)

r = np.array(r).T
num_genomes, num_genes = r.shape

print("Number of genomes: ", num_genomes, "\nNumber of genes: ", num_genes, "\n")

# Number of iterations for analysis
num_iterations = iterations

# Create a list to store results for each iteration
results = []

for iteration in range(num_iterations):
    print('Iteration: ', iteration + 1, end="\r")
    x, y = [], []
    result = np.zeros(num_genes)
    c = list(range(num_genomes))
    random.shuffle(c)
    for h, i in enumerate(c):
        x.append(h + 1)
        result += r[i]
        y.append(np.count_nonzero(result))
    results.append((x, y))

# Calculate the average parameters (k and gamma) from iterations
avg_k, avg_gamma = 0, 0
for (x, y) in results:
    pars, cov = curve_fit(f=heaps_law, xdata=x, ydata=y, p0=[0, 0], bounds=(-np.inf, np.inf))
    avg_k += pars[0]
    avg_gamma += pars[1]

avg_k /= num_iterations
avg_gamma /= num_iterations

# Determine if the pan-genome is open or closed based on gamma
pan_genome_status = "open" if avg_gamma > 0 else "closed"

# Create a single figure with two subplots
plt.figure(figsize=(10, 8))

# Subplot 1: Pan-genome and Modified Core-genome
plt.subplot(1, 1, 1)
title_font = FontProperties(weight='bold', size=14)

# Use a single color for pan-genome iterations
color = 'b'

for i, (x, y) in enumerate(results):
    if i == 0:
        plt.scatter(x, y, color=color, s=180, alpha=0.1)  # Bubble-like spheres with light color
    else:
        plt.scatter(x, y, color=color, s=180, alpha=0.1)  # Bubble-like spheres with light color

plt.plot(x, heaps_law(np.array(x), avg_k, avg_gamma), 'r-', label="Average Fit (Pan-genome)", linewidth=2)

# Initialize lists to store core-genome data
x_core, y_core = [], []

for num_selected_genomes in range(1, num_genomes + 1):  # Iterate from 1 to N
    avg_common_gene_counts = []
    if num_selected_genomes == 1:
        # Iteratively select 10 random genomes and count their genes
        for _ in range(num_iterations):
            selected_genomes = random.sample(range(num_genomes), num_selected_genomes)
            genes_count = [np.sum(r[i]) for i in selected_genomes]
            avg_genes_count = np.mean(genes_count)
            avg_common_gene_counts.append(avg_genes_count)
    else:
        for iteration in range(num_iterations):
            common_gene_counts = []
            for _ in range(num_selected_genomes):
                selected_genomes = random.sample(range(num_genomes), num_selected_genomes)
                common_genes = np.where(r[selected_genomes[0]] == 1)
                for genome_index in selected_genomes[1:]:
                    common_genes = np.intersect1d(common_genes, np.where(r[genome_index] == 1))
                common_gene_count = len(common_genes)
                common_gene_counts.append(common_gene_count)
            avg_common_gene_count = sum(common_gene_counts) / num_selected_genomes
            avg_common_gene_counts.append(avg_common_gene_count)
    x_core.extend([num_selected_genomes] * num_iterations)
    y_core.extend(avg_common_gene_counts)

# Perform average curve fitting for the modified core-genome plot
x_fit = np.array(list(range(1, num_genomes + 1)))  # Iterate from 1 to N
avg_y_fit = np.mean([y_core[i*num_iterations:(i+1)*num_iterations] for i in range(num_genomes)], axis=1)
pars, cov = curve_fit(f=heaps_law, xdata=x_fit, ydata=avg_y_fit, p0=[0, 0], bounds=(-np.inf, np.inf))

# Plot the modified core-genome with average fit curve
plt.scatter(x_core, y_core, color='orange', marker='o', s=180, alpha=0.1, edgecolor='orange')  # Bubble-like spheres with light color
plt.plot(x_fit, heaps_law(x_fit, pars[0], pars[1]), 'g-', linewidth=2)

xlabel_font = FontProperties(weight='bold', size=12)
ylabel_font = FontProperties(weight='bold', size=12)

plt.xlabel("Number of Genomes (N)", fontproperties=xlabel_font)
plt.ylabel("Number of Genes", fontproperties=ylabel_font)

plt.tight_layout()

# Save the figure as a high-resolution image (e.g., PNG or PDF) for publication
# Use os.path.join to create the full path including output_dir
png_output_path = os.path.join(output_dir, "pan_core_genome_plot.png")
pdf_output_path = os.path.join(output_dir, "pan_core_genome_plot.pdf")
plt.savefig(png_output_path, dpi=300, bbox_inches='tight')
plt.savefig(pdf_output_path, format='pdf', dpi=300, bbox_inches='tight')

# Save statistics for pan-genome to a text file
# Use os.path.join to create the full path including output_dir
stats_output_path = os.path.join(output_dir, "pan_genome_statistics.txt")
with open(stats_output_path, "w") as stats_file:
    stats_file.write(f"k = {avg_k}\n")
    stats_file.write(f"gamma = {avg_gamma}\n")
    stats_file.write(f"Pan-genome is {pan_genome_status}\n")