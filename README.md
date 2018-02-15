# miniMDS

miniMDS is a tool for inferring and plotting 3D structures from normalized Hi-C data, using a novel approximation to multidimensional scaling (MDS). It produces a single 3D structure from a Hi-C BED file, representing an ensemble average of chromosome conformations within the population of cells. Using parallelization, it is able to process high-resolution data quickly with limited memory requirements. The human genome can be inferred at kilobase-resolution within several hours on a desktop computer. Standard MDS results in inaccuracies for sparse high-resolution data, but miniMDS focuses on local substructures to achieve greater accuracy. miniMDS also supports interchromosomal structural inference. Together with Mayavi, miniMDS produces publication-quality images and gifs. 

## Installation

Requirements:
* python 2.7
* Python dependencies can be installed using
``pip install -r requirements.txt``
* The following optional dependencies can be installed manually:
    * [mayavi](http://docs.enthought.com/mayavi/mayavi/) (for plotting)
    * [ImageMagick](https://www.imagemagick.org/script/index.php) (for creating gifs)

## Testing

2/15/18: Test scripts need to be updated, ignore for now. 
On Linux, please run test.sh (in the scripts directory) and report any issues. (miniMDS is compatible with Mac, but the shell scripts only run on Linux.) 

## TLDR

``python minimds.py -o [output path] [BED path]``

## Usage

### Input file format

miniMDS uses intra- or inter-chromosomal BED files as input. Data must be normalized prior to use (for example, using [HiC-Pro](http://nservant.github.io/HiC-Pro/)). 

Format:

>chrA	bin1\_start	bin1\_end	chrB	bin2\_start	bin2\_end	normalized\_contact\_frequency

Example - chr22 intra-chromosomal data at 10-Kbp resolution:

>chr22	16050000	16060000	chr22	16050000	16060000	12441.5189291
> 
>...

### Intra-chromosomal miniMDS

Intra-chromosomal analysis is performed using minimds.py.

To view help:

``python minimds.py -h``

By default, partitioned MDS is used:

``python minimds.py GM12878_combined_22_100kb.bed``

Full MDS is recommended for low-resolution high-quality (not sparse) files:

``python minimds.py --full GM12878_combined_22_100kb.bed``

Structures are not saved by default. Use the -o option with the path where you want to save the structure.

``python minimds.py -o GM12878_combined_22_100kb_structure.tsv GM12878_combined_22_100kb.bed``

Structures are saved to tsv files. The header contains the name of the chromosome, the resolution, and the starting genomic coordinate. Each line in the file contains the genomic bin number followed by the 3D coordinates (with "nan" for missing data). 

Example - chr22 at 10-Kbp resolution:

>chr22
> 
>10000
> 
>16050000
> 
>0	0.589878298045	0.200029092422	0.182515056542
> 
>1	0.592088232028	0.213915817254	0.186657230841
> 
>2	nan	nan	nan
> 
>...

0 corresponds to the bin 16050000-16060000, 1 corresponds to the bin 16060000-16070000, etc. 

#### Parameters (optional)

#### Resolution ratio

miniMDS first infers a global intrachromosomal structure at low resolution, which it uses as a scaffold for high-resolution inference. By default a resolution ratio of 10 is used. So if your input file is 100-kb resolution, a 1-Mb structure will be used for approximation. The resolution ratio can be changed with the l option. 

``python minimds.py -l 20 -o GM12878_combined_22_10kb_structure.tsv GM12878_combined_22_10kb.bed``

The resolution you choose for the low-res file depends on your tradeoff between speed and accuracy. Lower resolutions (i.e. higher ratios) are faster but less accurate. For now, the high resolution must be a factor of the low resolution. For example, a 500-Kb-resolution file can be used to infer a 100-Kb-resolution structure, but a 250-Kb-resolution file cannot. 

##### Controlling the number of partitions

The miniMDS algorithm creates partitions in the high-resolution data and performs MDS on each partition individually. A greater number of partitions can increase speed but also reduce accuracy. On the other hand, for very sparse data a greater number of partitions can actually increase accuracy. If your output appears "clumpy", increase the number of partitions.

The number of partitions cannot be set directly because partitions are created empirically to maximize clustering of the data. However, the degree of clustering of the data can be tweaked with the following parameters:

>-m: minimum partition size (as a fraction of the data). Default = 0.05
>
>-p: smoothing parameter (between 0 and 1). Default = 0.1

Make these parameters smaller to increase the number of partitions. For very high resolution data (such as 5-Kbp), m=0.01 and p=0.01 is recommended:

``python minimds.py -o GM12878_combined_22_5kb_structure.tsv -m 0.01 -p 0.01 GM12878_combined_22_5kb.bed``

You can limit the maximum RAM (in Kb) used by any given partition using -R (default = 32000):

``python minimds.py -o GM12878_combined_22_5kb_structure.tsv -R 50000 GM12878_combined_22_5kb.bed``

##### Number of threads

miniMDS uses multithreading to achieve greater speed. By default, 3 threads are requested, because this is safe for standard 4-core desktop computers. However, the number of threads used will never exceed the number of processors or the number of partitions, regardless of what is requested. You can change the number of requested threads using -n.

For example, to run miniMDS with four threads:

``python minimds.py -o GM12878_combined_22_10kb_structure.tsv -n 4 GM12878_combined_22_10kb.bed``

##### Scaling factor

The scaling factor a describes the assumed relationship between contact frequencies and physical distances: distance = contact_frequency^(-1/a). The default value is 4, based on Wang et al 2016. You can change the scaling factor using -a. 

``python minimds.py -o GM12878_combined_22_10kb_structure.tsv -a 3 GM12878_combined_22_10kb.bed``

a can be any value >1, including non-integer.

##### Classical MDS

Classical MDS (cMDS), also called principal coordinates analysis, is a variant of MDS that is faster under certain circumstances. The miniMDS tool supports cMDS but NOT with partitioned MDS. Use the --classical option. 

``python minimds.py --classical GM12878_combined_22_10kb.bed``

This mode is mainly used for testing. 

### Inter-chromosomal miniMDS

Inter-chromosomal analysis is performed using minimds_inter.py

To view help:

``python minimds_inter.py -h``

The usage of minimds_inter.py is similar to minimds.py, however inter-chromosomal files are required in addition to intra-chromosomal. To avoid entering filenames separately for each chromosome, you must name your files using a standard format.

Intra-chromosomal format:

>{prefix}\_{ChrA}\_{resolution}{kb or Mb}.bed

Example:

>GM12878_combined_22_100kb.bed

Inter-chromosomal format:

>{prefix}\_{ChrA}\_{ChrB}_{resolution}{kb or Mb}.bed

where A is before B in:

>1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X

Example:

>GM12878_combined_21_22_100kb.bed

Enter the prefix and resolution of the inter-chromosomal and intra-chromosomal files, respectively:

``python minimds_inter.py [inter-chromosomal file prefix] [intra-chromosomal file prefix] [inter-chromosomal resolution] [intra-chromosomal resolution]``

For example, if your files are stored in the directory _data_:

``python minimds_inter.py data/GM12878_combined_interchromosomal data/GM12878_combined_intrachromosomal 1000000 10000``

Because of the challenges of inter-chromosomal inference, it is recommended that a resolution no greater than 1-Mbp be used for inter-chromosomal data. 

#### Other parameters (optional)

All of the parameters from minimds.py are also available for minimds_inter.py

###### Specifying chromosomes

By default, minimds_inter.py uses all human chromosomes other than Y. You can specify chromosomes using the option -c.

To perform interchromosomal analysis on chromosomes 1 and 8:

``python minimds_inter.py -c 1 8 data/GM12878_combined_interchromosomal data/GM12878_combined_intrachromosomal 1000000 10000``

Note: it is often necessary to use this option if you are using a genome other than human, so that it won't search for chromosomes that don't exist.

### Plotting

Read a structure into a Cluster object:

``cluster = data_tools.clusterFromFile(path)``

Example:

``cluster = data_tools.clusterFromFile("GM12878_combined_22_100kb_structure.tsv")``

Create an interactive 3D plot in Mayavi. (Mayavi allows you to rotate the image and save a view.)

``plotting.plot_cluster_interactive(cluster, color=(1,0,0), radius=None)``

By default, the radius is the to-scale radius of heterochromatin. 

Multiple clusters can be plotted simultaneously:

    chroms = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, X)
    clusters = [data_tools.clusterFromFile("GM12878_combined_{}_100kb_structure.tsv".format(chrom) for chrom in chroms)]
    plotting.plot_clusters_interactive(clusters)

plotting.py has 23 built-in colors designed to be as different to the human eye as possible. By default, these colors are used when plotting multiple clusters. You can also specify a list of colors:

    chroms = (1, 2)
    clusters = [data_tools.clusterFromFile("GM12878_combined_{}_100kb_structure.tsv".format(chrom) for chrom in chroms)]
    plotting.plot_clusters_interactive(clusters, colors=[(1,0,0), (0,0,1)])

The radius can also be specified, as above. 

The option _cut_ creates a cross-section of the plot. For example, this is useful for viewing the interior of the nucleus.

    chroms = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, X)
    clusters = [data_tools.clusterFromFile("GM12878_combined_{}_100kb_structure.tsv".format(chrom) for chrom in chroms)]
    plotting.plot_clusters_interactive(clusters, cut=True)

A plot can be saved as a gif:

``plotting.plot_cluster_gif(cluster, outname, color=(1,0,0), radius=None, increment=10)``

A smaller value of _increment_ will lead to a smoother gif.

Multiple clusters can also be plotted in a single gif:

``plotting.plot_clusters_gif(clusters, outname, colors=default_colors, radius=None, increment=10)``
