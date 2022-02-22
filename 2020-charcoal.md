---
output:
  pdf_document: default
  html_document: default
---
# Charcoal: filtering contamination in metagenome-assembled genome bins
authors:
Taylor Reiter, N. Tessa Pierce-Ward, Luiz Irber, Erich Schwarz, Other People?, C. Titus Brown

keywords:

  - metagenome-assembled genome bin
  - contamination
  - FracMinHash
  - minsetcov

## Abstract

## Introduction

Metagenomic sequencing has expanded our knowledge of microbial communities and their diversity [@doi:10.1038/nmicrobiol.2016.48; @doi:10.1038/s41587-020-0718-6; @doi:10.1038/s41564-021-00928-6].
*De novo* metagenome analysis has generated thousands of draft genomes, termed metagenome-assembled genomes (MAGs), from organisms from diverse environments [@doi:10.1038/nature02340; @doi:10.1016/j.cell.2019.01.001; @doi:10.1038/s41587-020-0603-3; @doi:10.1038/s41587-020-0718-6].
Recently, large-scale re-analysis efforts have led to a rapid expansion in draft genomes in public repositories like the European Nucleotide Archive [@doi:10.1038/s41587-020-0603-3] and the Joint Genome Institute IMG/M [@doi:10.1038/s41587-020-0718-6].
Increased observation of draft genomes across the tree of life better enables researchers to contextualize new sequencing data and the roles that microorganisms play in diverse metabolic processes [@doi:10.1038/s41564-017-0012-7; @doi:10.1093/nar/gkz1002].

MAG inference relies on assembly and binning of metagenomic sequencing data.
Assembly produces long contiguous sequences by identifying overlaps between short sequencing reads, while binning groups assembled sequences into MAGs using read coverage and tetramernucleotide frequency.
Both processes are subject to biases that can reduce the completeness of or increase the contamination in a MAG: low sequencing coverage or high genomic variation causes short read assemblers to break contiguous sequences into shorter pieces (CITE), which decreases the signal for and accuracy of binning (CITE).
Commonly, the completeness and purity of MAGs is estimated through the presence and sequence composition of single-copy marker genes [@doi:10.1101/gr.186072.114] [@doi:10.1093/bioinformatics/btv351], with MAGs that reach >90% completeness and <5% contamination considered high quality [@doi:10.1101/gr.186072.114].
Single-copy marker genes are sets of genes that are present once in a genome of almost all members of a taxonomic group [@doi:10.1101/gr.186072.114; @doi:10.7717/peerj.243].
Using these genes to estimate contamination leads to two important biases.
First, given the assumption that marker genes are universally present in genomes, if a marker gene resides on a contaminant sequence but no other sequence for that marker gene is present in the genome, it will not be detected as contamination [@doi:10.1101/gr.186072.114; @doi:10.3389/fmicb.2017.02264].
When a MAG is substantially complete, this may lead to a small underestimation in contamination (~3% [@doi:10.1101/gr.186072.114]), but as completeness decreases, contamination may be substantially underestimated [@doi:10.3389/fmicb.2017.02264].
Second, contiguous sequences which do not contain marker sequences are not included among contamination estimates [@doi:10.1101/gr.186072.114] (note checkm called out plasmids and phages for this bias specifically).

Given these biases, methods that do not rely solely on marker genes may be better suited for contamination estimation.
Many tools have recently been developed that rely on different strategies for detecting contamination [@doi:10.1186/s13059-022-02619-9].
GUNC expands beyond marker genes and uses the full complement of genes in a genome to identify contiguous sequences that fall outside of the predicted lineage (CITE: gunc [@doi:10.1186/s13059-021-02393-0]).
It estimates the distribution of taxonomic assignments within and across contiguous sequences to identify contaminants even among sequences that are poorly taxonomically labelled [@doi:10.1186/s13059-021-02393-0].
Conterminator detects cross-kingdom contamination using all-vs-all alignment and is geared toward contamination detection in large collections of genomes (e.g. databases) [@doi:10.1186/s13059-020-02023-1].
Sequences are considered contaminants when at least 100 nucleotides and no more than 20 kb share a sequence identity of at least 90% in genomes that originate from different kingdoms [@doi:10.1186/s13059-020-02023-1].
In addition, conterminator finds contamination in protein sequences by identifying protein clusters that contain sequences from cross-kingdom members [@doi:10.1186/s13059-020-02023-1].
 
Removing contamination is a separate problem from estimating its extent.
Tools like RefineM [@doi:10.1038/s41564-017-0012-7], MAGPpurify [@doi:10.1038/s41586-019-1058-x], and BlobTools [@doi:10.12688/f1000research.12232.1] rely on a combination of GC content, tetramernucleotide frequency, read coverage, phylogenetic or clade markers, conspecific sequences, and known contaminants to flag contaminant contiguous sequences for removal.
The inclusion of read coverage profiles necessitates the use of the original sequencing reads to produce coverage profiles in BAM format for the genomes.
Given that sequencing data and BAM files are usually orders of magnitude larger than genome sequences, this requirement limits the utility of such approaches, especially for large-scale genome databases. 
Tools like GUNC are poised for contamination removal without relying on read coverage profiles, but as of yet such approaches are unimplemented [@doi:10.1186/s13059-021-02393-0].

Long k-mers capture relatedness between organisms, where a k=31 captures species-level similarity (CITE: METAPALETTE).
K-mers offer an alternative metric to identify contamination, especially in sequences lacking marker genes. 
Here we describe Charcoal, an automated method for filtering contaminant contiguous sequences from MAGs. 
We show...
We show...
We envisage that charcoal will complement marker gene-based approaches for contamination estimation, removing problematic sequences before they are further analyzed or propagated in public databases. 

## Materials and methods

### Overview

Charcoal identifies and removes contamination in metagenome-assembled genomes or other genomes that may contain contaminant sequences using k-mer based methods (**FIGURE 1A**).
Charcoal is a snakemake workflow developed around the tool sourmash, which enables Jaccard similarity estimation between sequence sets of different sizes using FracMinHash k-mer sketches (CITE: Snakemake, F1000, JOSS, GATHER). 
A k-mer is a nucleotide sequence of length _k_. 
When _k_ is sufficiently long, k-mers are generally specific to a taxonomic lineage (CITE: METAPELLETE, GATHER PAPER, TESSA).
Taking advantage of this property of k-mers, charcoal identifies majority and minority lineages for each contiguous sequence in a genome and removes contiguous sequences belonging to minority lineages when those lineages occur below a user-defined taxonomic threshold.
By default, charcoal filters contiguous sequences that are assigned a different order than that of the majority lineage.

To identify contamination, charcoal first creates a FracMinHash sketch for each contiguous sequence in an input genome.
Charcoal then identifies all genomes in a database ("reference genomes") that share sequence overlap with the input genome using sourmash `prefetch`.
Subsequent operations subset the original database to include only the genomes identified by `prefetch`, reducing search volumes.
Using these matches, charcoal uses sourmash `gather` to identify the minimum set of genomes that cover (or contain) the k-mers in each contiguous sequence (CITE: gather). 
Charcoal then determines the taxonomic lineage of each contiguous sequence using a lineage spreadsheet that records the taxonomy of each reference genome in the database; the taxonomic assignment occurs at the lowest common ancestor of all taxonomic assignments given to a contiguous sequence.

Charcoal compares the taxonomic lineage of each contiguous sequence against the lineage of the input genome.
If the contiguous sequence has a different lineage before or at the configured taxonomic rank (order by default) than that of the majority lineage, the contiguous sequence is considered a contaminant. 
Charcoal reports the majority lineage, the fraction of identified hashes and the fraction of identified hashes that match the majority lineage, an estimate of the number of contaminant base pairs at each level of taxonomy up to the filter taxonomic rank, the number of contiguous sequences and an estimate number of base pairs that were ignored by charcoal, and the number of contiguous sequences and an estimate number of base pairs that were no identified in the data base.

The input genome lineage can be user-provided, or it can be determined by charcoal via majority vote of all lineages assigned to all contiguous sequences.
If the lineage is determined by charcoal, by default a minimum 10% of the input genome must have been assigned a taxonomic lineage, and 20% of those sequences assigned a taxonomic lineage must match to the majority lineage.
If these specifications are not met, charcoal will not decontaminate the input genome unless the user specifies a lineage.
The user can optionally specify a lineage (e.g. `d__Eukaryota`), and charcoal will remove contiguous sequences that have a lineage different from the user-specified one.
This allows charcoal to remove contiguous sequences from an input genome when some contiguous sequences from that genome occur in the database and when the input genome is not related to anything in a database.
Charcoal reports whether the provided lineage agrees with k-mer classification at or above the genus level.

After the initial stage of contaminant detection, charcoal can perform additional tasks to verify, summarize, or remove the contaminant sequences.
Contaminant verification downloads the reference genome sequence for any genome that was detected among the input genome contiguous sequences and aligns the contiguous sequences against those genomes using mashmap (CITE: mashmap).
Contaminant removal separates "clean" from "dirty" contiguous sequences and outputs each set into a FASTA file.
Importantly, charcoal will not remove a contiguous sequence if it is unidentifiable, whether it is too short to be sketched or does not contain sequences in the reference database. 
While these contiguous sequences could still be contaminants, charcoal assumes contiguous sequences are clean for which it has no information. 
Therefore, charcoal will fail to detect contamination for very short contiguous sequences which contain no selected k-mers, as well as contiguous sequences with novel DNA content.
Lastly, charcoal has a report feature that summarizes and visualizes the taxonomic lineages detected in each input genome as well as the alignments between the input genome and reference genomes (**FIGURE**).

Identifying all lineages in a reference database present in the input genome with sourmash `prefetch` is the most compute intensive step in the decontamination process.
The RAM and CPU time needed for this step depend on the size of the database used for decontamination.

+ trick to remove exact matches
+ statement that charcoal is database dependent (probably can go in last paragraph with RAM/CPU considerations)

### Availability

Charcoal is written in python3 and can be installed via conda or pip. 
The core algorithms (contaminant detection and removal) depend on sourmash and snakemake. 
Contaminant verification depends on lxml and mashmap. 
Reporting depends on mummer, papermill, notebook, and plotly.
The source code is available at github.com/dib-lab/charcoal.
ZENODO DOI.

### Datasets and benchmarking

## Results


## Future work

+ strain heterogeneity estimation or strain insertion in newick tree using prefetch & clean results
+ completeness estimation, potentially via comparison to conspecifics via genome size
+ database with eukaryotic sequences
+ estimation via protein sequences for divergent organisms not well represented in databases

## Notes

+ fastani is an accepted way to do average nucleotide identity calculate, and it relies on k-mers.
+ from checkm paper:
> Bias in genome quality estimates: Quality estimates based on individual marker genes or collocated marker sets exhibit a bias resulting in completeness being overestimated and contamination being underestimated. This bias is the result of marker genes residing on foreign DNA that are otherwise absent in a genome being mistakenly interpreted as an indication of increased completeness as opposed to contamination. 
+ we don't deal with strain heterogeneity, as this occurs below the species-level aggregation in the LCA

## old results outline:

+ charcoal estimates low contamination in non-representative GTDB genomes
+ charcoal assigns correct taxonomy to all non-representative GTDB genomes
+ charcoal vs. checkm
  + charcoal vs. checkm: mgnify, tara
  + checkm on charcoal clean
  + prokka on charcoal dirty
+ charcoal vs. refineM and magpurify
+ verification of contamination/contam case studies
  + case studies
  + user-provided lineages
  + cDBG/spacegraphcats?
  + user-specified lineages
