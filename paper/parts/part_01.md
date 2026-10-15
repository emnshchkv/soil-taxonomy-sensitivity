# Часть 1: с. 1–3

[с. 1]
PLOS COMPUTATIONAL BIOLOGY

RESEARCH ARTICLE

RESCRIPt: Reproducible sequence taxonomy
reference database management
Michael S. Robeson, II ID1, Devon R. O’Rourke ID2, Benjamin D. Kaehler ID3,
Michal Ziemski ID4, Matthew R. Dillon ID2, Jeffrey T. Foster2, Nicholas A. Bokulich ID4*

a1111111111
a1111111111
a1111111111
a1111111111
a1111111111

OPEN ACCESS
Citation: Robeson MS, II, O’Rourke DR, Kaehler
BD, Ziemski M, Dillon MR, Foster JT, et al. (2021)
RESCRIPt: Reproducible sequence taxonomy
reference database management. PLoS Comput
Biol 17(11): e1009581. https://doi.org/10.1371/
journal.pcbi.1009581
Editor: Mihaela Pertea, Johns Hopkins University,
UNITED STATES
Received: April 22, 2021
Accepted: October 21, 2021
Published: November 8, 2021
Copyright: © 2021 Robeson, II et al. This is an
open access article distributed under the terms of
the Creative Commons Attribution License, which
permits unrestricted use, distribution, and
reproduction in any medium, provided the original
author and source are credited.
Data Availability Statement: Data reporting: All
data analysed herein, were retrieved either using
RESCRIPt (for SILVA [https://www.arb-silva.de/]
and NCBI [https://www.ncbi.nlm.nih.gov/genbank/
]) data, or by direct download of release data (for
UNITE [https://unite.ut.ee/], Greengenes [ftp://
greengenes.microbio.me/greengenes_release/gg_
13_5/], and GTDB [https://gtdb.ecogenomic.org/]
or by direct download (for BOLD [https://www.
boldsystems.org/] data; accessed July 1, 2020 and
updated August 8, 2020). Availability of data and
materials: Workflows and data from our

1 University of Arkansas for Medical Sciences, Department of Biomedical Informatics, Little Rock, Arkansas,
United States of America, 2 Pathogen and Microbiome Institute, Northern Arizona University, Flagstaff,
Arizona, United States of America, 3 School of Science, University of New South Wales, Canberra, Australia,
4 Laboratory of Food Systems Biotechnology, Institute of Food, Nutrition, and Health, ETH Zürich,
Switzerland
* nicholas.bokulich@hest.ethz.ch

Abstract
Nucleotide sequence and taxonomy reference databases are critical resources for widespread applications including marker-gene and metagenome sequencing for microbiome
analysis, diet metabarcoding, and environmental DNA (eDNA) surveys. Reproducibly generating, managing, using, and evaluating nucleotide sequence and taxonomy reference
databases creates a significant bottleneck for researchers aiming to generate custom
sequence databases. Furthermore, database composition drastically influences results,
and lack of standardization limits cross-study comparisons. To address these challenges,
we developed RESCRIPt, a Python 3 software package and QIIME 2 plugin for reproducible
generation and management of reference sequence taxonomy databases, including dedicated functions that streamline creating databases from popular sources, and functions for
evaluating, comparing, and interactively exploring qualitative and quantitative characteristics across reference databases. To highlight the breadth and capabilities of RESCRIPt, we
provide several examples for working with popular databases for microbiome profiling
(SILVA, Greengenes, NCBI-RefSeq, GTDB), eDNA and diet metabarcoding surveys
(BOLD, GenBank), as well as for genome comparison. We show that bigger is not always
better, and reference databases with standardized taxonomies and those that focus on type
strains have quantitative advantages, though may not be appropriate for all use cases. Most
databases appear to benefit from some curation (quality filtering), though sequence clustering appears detrimental to database quality. Finally, we demonstrate the breadth and extensibility of RESCRIPt for reproducible workflows with a comparison of global hepatitis
genomes. RESCRIPt provides tools to democratize the process of reference database
acquisition and management, enabling researchers to reproducibly and transparently create
reference materials for diverse research applications. RESCRIPt is released under a permissive BSD-3 license at https://github.com/bokulich-lab/RESCRIPt.

PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1009581 November 8, 2021

1 / 37

[с. 2]
PLOS COMPUTATIONAL BIOLOGY

benchmarks can be found at https://github.com/
bokulich-lab/db-benchmarks-2020 and https://
github.com/devonorourke/COIdatabases/. Code
reporting: Source code, installation and usage
instructions, and tutorials for RESCRIPt can be
found at the project page: https://github.com/
bokulich-lab/RESCRIPt.
Funding: The authors received no specific funding
for this work.
Competing interests: The authors declare that they
have no competing interests.

Reproducible reference database management

Author summary
Generating and managing sequence and taxonomy reference data presents a bottleneck to
many researchers, whether they are generating custom databases or attempting to format
existing, curated reference databases for use with standard sequence analysis tools. Evaluating database quality and choosing the “best” database can be an equally formidable challenge. We developed RESCRIPt to alleviate this bottleneck, supporting reproducible,
streamlined generation, curation, and evaluation of reference sequence databases.
RESCRIPt uses QIIME 2 artifact file formats, which store all processing steps as data provenance within each file, allowing researchers to retrace the computational steps used to
generate any given file. We used RESCRIPt to benchmark several commonly used
marker-gene sequence databases for 16S rRNA genes, ITS, and COI sequences, demonstrating both the utility of RESCRIPt to streamline use of these databases, but also to evaluate several qualitative and quantitative characteristics of each database. We show that
larger databases are not always best, and curation steps to reduce redundancy and filter
out noisy sequences may be beneficial for some applications. We anticipate that RESCRIPt
will streamline the use, management, and evaluation/selection of reference database materials for microbiomics, diet metabarcoding, eDNA, and other diverse applications.

This is a PLOS Computational Biology Software paper.

Introduction
Marker-gene amplicon and metagenome sequencing have become attractive methods for
characterizing microbial community composition and function [1,2] in human health [3–5]
and agriculture [6–8], as well as macroorganism diversity through diet metabarcoding studies
[9–11] and environmental DNA (eDNA) surveys [12–15]. Taxonomic classification is often a
primary goal in marker-gene and metagenome sequencing studies to identify the composition
of a mixed community, or to detect species of interest (e.g., pathogens or invasive species).
This is accomplished by comparing the observed sequences to a reference database consisting
of target marker-gene or genome sequences from known species. The selection of a reference
database can significantly impact both marker-gene and metagenome sequencing results
[16,17], and methods for assessing database quality and fitness for a given sample type or
hypothesis remain an undermet need.
Identification of Bacteria and Archaea is most commonly performed using the 16S rRNA
gene, due to its historical use as a phylogenetic marker [18,19] and the existence of curated reference databases [20,21]. The SILVA [20,22] rRNA gene database and Greengenes [21,23] 16S
rRNA gene database are commonly used for identifying Bacteria and Archaea, containing
curated taxonomies, sequences, and phylogenies. More recently, the Genome Taxonomy Database (GTDB) was developed with the intent to provide a standardized bacterial and archaeal
taxonomy based on genome phylogeny [24,25], and provides 16S rRNA reference sequences.
NCBI-RefSeq also provides several targeted loci sequence databases from curated records,
including Internal Transcribed Spacer (ITS), and both the small and large sub-unit (SSU &
LSU) rRNA genes [26]. Non-16S genes are also attractive targets for bacterial and archaeal species identification due to the degree of species resolution that they afford, but their application
is limited by the relative lack of curated reference materials [27–29].

PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1009581 November 8, 2021

2 / 37

[с. 3]
PLOS COMPUTATIONAL BIOLOGY

Reproducible reference database management

Fungal classification is most commonly performed using the ITS domain, the designated
fungal “barcode of life”, though the SSU and LSU rRNA genes are also common targets [30].
Both NCBI-RefSeq [26] and the UNITE database [31] provide curated ITS sequences from
fungi and other eukaryotes, as well as the RDP Warcup fungal ITS training set [32], which was
prepared from an earlier release of the UNITE+INSD database. Both SILVA [22] and RDP
[33] provide LSU databases for fungal sequence classification. NCBI RefSeq releases databases
for both fungal SSU and LSU [26].
Identification of other eukaryotes, including for diet metabarcoding, microbial eukaryote,
and eDNA surveys, is commonly accomplished using the mitochondrial cytochrome oxidase
subunit I (COI) gene for metazoa [34–36], 18S rRNA gene or other rRNA gene subunits [37],
ITS2 and chloroplast trnL (UAA) intron [38–40] for plants, 12S rRNA for fish [41,42], and a
variety of other clade-specific marker genes. For some of these marker genes, curated reference
databases exist, such as BOLD for COI [34] and PLANiTS for plant ITS2 [40], but for others
the process of generating custom reference databases poses a research bottleneck.
Taxonomic profiling studies rely on high-quality sequence taxonomy reference databases.
However, errors in public sequence databases are well documented [12,43,44] and can lead to
misclassification errors in downstream results [12]. Different reference databases can yield
widely different classification results for biological data, but standards are lacking to objectively assess the quality of individual databases [45]. Revisions to taxonomic naming [46–51]
and the rapid pace at which new sequences and genomes are added to public databases mean
that curated reference releases may lag behind [52]. Additionally, issues with amplicon length
and sequence heterogeneity can limit the ability to identify species, especially from short
marker-gene sequences or metagenome fragments [53]. Hence, many researchers choose to
perform additional curation to focus on type strains [54], quality filtering [14,55], or construct
environment-specific databases that are constrained to contain species found within a given
environment [54,56–61]. Database customization is also often performed to add new accessions that are absent in some database releases to increase database coverage [52], or to incorporate outgroups [14]. However, generating such databases can be technically challenging,
subjective, and difficult to document, leading to issues with transparency and reproducibility,
and limiting the ability of many researchers to acquire appropriate reference materials for
their studies, or leading to reliance on proprietary resources and services (limiting scientific
transparency and increasing research costs). Sequence curation is a significant hurdle in this
process, as taxonomic misannotations, sequence errors, and other errors in existing (and
inchoate) reference sequence databases reduce the accuracy of taxonomic classifiers that rely
on these data [43,44,62–64]. For example, inconsistent genus-level annotation of identical
sequences labeled as either ‘Escherichia’, ‘Shigella’, or even as the combined group ‘EscherichiaShigella’ [43] can result in queried sequences being incorrectly classified with their last common ancestor (LCA) family label “unclassified Enterobacteriaceae” instead of a more informative genus label.
The need for scientific results to be reproducible, replicable, and transparent has taken on
new urgency in the digital age [65]. On the one hand, increasing experimental and analytical
complexity pose mounting challenges to effective documentation and sharing of methodological procedures and results [65–68]. On the other hand, digital tools present opportunities to
address these challenges, and various reporting standards have been published to guide
researchers in reporting and publishing new types of data, software, and other resources [69–
71]. Following guidelines such as these is important for reporting, but also for standardization
of methods during data reuse and metaanalysis. Given the fundamental importance of reference databases to reporting results from marker-gene and metagenome experiments, it is

PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1009581 November 8, 2021

3 / 37

