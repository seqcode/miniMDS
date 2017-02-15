# miniMDS

miniMDS is a software tool for inferring and plotting 3D structures from normalized Hi-C data, using an efficient approximation to multidimensional scaling (MDS). It produces a single 3D structure from a Hi-C BED file, representing an ensemble average of chromosome conformations within the population of cells. By performing structural inference at multiple resolutions, it is able to process high-resolution data quickly with limited memory requirements. Kilobase-resolution structures can be inferred within several hours on a desktop computer. Standard MDS results in inaccuracies for sparse high-resolution data, but miniMDS focuses on local substructures to achieve greater accuracy. miniMDS also supports interchromosomal structural inference. Together with Mayavi, miniMDS produces publication-quality images and gifs. 

## Installation

Tested on python 2.7

Prerequisites:
* numpy
* scikit-learn
* pymp
* mayavi (optional; for plotting)
* scipy (optional; for creating figures from paper)
* matplotlib (optional; for creating figures paper)

## Usage

### Input file format

miniMDS uses intra- or inter-chromosomal BED files as input. Data must be normalized prior to use (for example, using https://bitbucket.org/mirnylab/hiclib). 

Format:

chrA	bin1\_start	bin1\_end	chrB	bin2\_start	bin2\_end	normalized\_contact\_frequency

Example - chr22 intra-chromosomal data at 10-Kbp resolution:

chr22	16050000	16060000	chr22	16050000	16060000	12441.5189291

### Intra-chromosomal structural inference

To view help:

python minimds.py -h

By default, standard MDS (not the miniMDS algorithm) is used:

python minimds.py GM12878\_combined\_22\_10kb.bed

Structures are not saved by default. Use the -o option with the path where you want to save the structure.

python minimds.py -o GM12878\_combined\_22\_10kb_structure.tsv GM12878\_combined\_22\_10kb.bed

Structures are saved to tsv files. The header contains the name of the chromosome, the resolution, and the starting genomic coordinate. Each line in the file contains the point number followed by the 3D coordinates (with "nan" for missing data). 

Example - chr22 at 10-Kbp resolution:

chr22

10000

16050000

0	0.589878298045	0.200029092422	0.182515056542

1	0.592088232028	0.213915817254	0.186657230841

2	nan	nan	nan

...

Structures can be read into Cluster objects:

cluster = data\_tools.clusterFromFile("GM12878\_combined\_22\_100kb_structure.tsv")

To run the miniMDS algorithm, you must have a normalized BED file at a lower resolution than the BED file you want to infer. For example, to use a 100-Kbp-resolution BED file to aid in the inference of a 10-Kbp-resolution file:

python minimds.py -l GM12878\_combined\_22\_100kb.bed -o GM12878\_combined\_22\_10kb_structure.tsv GM12878\_combined\_22\_10kb.bed

The resolution you choose for the low-res file depends on your tradeoff between speed and accuracy. Lower resolutions are faster but less accurate. 
