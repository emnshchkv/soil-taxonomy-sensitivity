# Часть 2: с. 4–7

[с. 4]
PLOS COMPUTATIONAL BIOLOGY

Reproducible reference database management

critical that principles such as these be followed by researchers when they acquire, modify, and
use reference data.
To address the need for reproducible bioinformatics workflows to streamline database generation and curation, we developed RESCRIPt (REference Sequence annotation and CuRatIon
Pipeline; https://github.com/bokulich-lab/RESCRIPt). Below we describe the RESCRIPt software package, and demonstrate its use via a series of benchmarks to evaluate sequence and taxonomy information in several widely used reference databases.

Results
Contents
1. Software Description
2. Comparison of 16S rRNA Gene Sequence Databases
3. Effects of processing steps on the SILVA 16S rRNA gene database
4. Effect of Clustering on Sequence and Taxonomic Information: lessons from the Greengenes
16S rRNA gene database
5. Reference Curation Improves Taxonomic Classification: lessons from the UNITE Fungal
ITS database
6. Clustering and primer-region trimming effects on a BOLD COI gene database
7. Comparison of metazoan COI gene sequences in BOLD and GenBank
8. Fetching Reference Genomes for Classification

Software description
RESCRIPt (https://github.com/bokulich-lab/RESCRIPt) is a Python 3 package for retrieving,
filtering, and evaluating nucleotide sequence and taxonomic data (Fig 1). It was motivated by
the need for scientists to transparently and reproducibly generate and curate reference
sequence databases, to facilitate interoperability and comparison downstream. RESCRIPt was
implemented as a QIIME 2 [72] plugin, in order to incorporate the integrated data provenance
and multiple user interfaces of QIIME 2. Hence all processing steps used to generate a database
are recorded in provenance that is stored both in the database files as well as in all downstream
results, enhancing scientific transparency, reproducibility of results, and replication of the
database (and the processing steps used in its creation) by other researchers, following the
FAIR data principles to make data findable, accessible, interoperable and reusable [69].
RESCRIPt enables efficient and transparent construction of reference databases for any
amplicon targets for which source data exist, as well as for full genomes from NCBI. RESCRIPt
currently supports a variety of methods for acquiring, formatting, and evaluating nucleotide
sequence and taxonomy databases (Fig 1). These include methods for automated download of
sequences and taxonomy from SILVA (selecting marker gene and release version) and NCBI
GenBank (selection based on NCBI query or a list of accession numbers). Additional functions
allow traceable formatting (e.g., sequence manipulation, taxonomy parsing), filtering (e.g.,
automatic quality and length filtering), and modification of data (e.g., to merge taxonomies
and dereplicate or cluster sequence/taxonomy databases). Several functions allow evaluation
of information content in sequence and taxonomy databases, e.g., based on sequence and taxonomy label entropy and cross-validated taxonomic classification of sequences (to simulate

PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1009581 November 8, 2021

4 / 37

[с. 5]
PLOS COMPUTATIONAL BIOLOGY

Reproducible reference database management

Fig 1. Current RESCRIPt functionality for processing and curating reference sequence data. Arrows indicate suggested workflows. Dotted arrows and
edges indicate optional steps for customized workflows.
https://doi.org/10.1371/journal.pcbi.1009581.g001

classification accuracy). Future development plans are discussed below (see Future Goals), and
further details on functionality are described in the Methods section. Complete usage details
and tutorials are available on the project’s source code repository (https://github.com/
bokulich-lab/RESCRIPt).
To demonstrate the diverse applications of RESCRIPt for sequence and taxonomy analysis,
we used RESCRIPt to benchmark several commonly used reference sequence databases,
including several popular databases of bacterial 16S rRNA genes, fungal ITS sequences, and
eukaryote COI genes that are commonly used for diet metabarcoding and eDNA studies.
These include side-by-side comparisons to evaluate relative information and performance
characteristics, and the impacts of various sequence curation steps on database characteristics.
An example tutorial is provided for the acquisition and construction of a 12S rRNA marker
gene reference database (S1 File). Finally, we demonstrate the application of RESCRIPt for
reproducible and extensible genomics analysis workflows via a comparison of hepatitis
genomes from several global sources.

Comparison of 16S rRNA gene sequence databases
Researchers investigating bacterial and archaeal community compositions using 16S rRNA
gene sequences are faced with myriad options for reference sequence databases. Those using
non-16S genes will be quick to remind them that having choices is a good problem to have—
but selecting the “best” reference materials for a specific task is still indeed a problem faced by
many researchers. Although the Greengenes database [21,23] is popular among the microbiome research community, the last release was in 2013 and much has changed in the world of

PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1009581 November 8, 2021

5 / 37

[с. 6]
PLOS COMPUTATIONAL BIOLOGY

Reproducible reference database management

microbial taxonomy in the interim. The SILVA database [22] has been a popular alternative,
boasting a regular release cycle, curated taxonomy and sequences, and a large database size.
More recently, the GTDB project [24] seeks to create a standardized bacterial and archaeal taxonomy based on genome phylogeny, making it an attractive database for some researchers.
Meanwhile, many other options exist, such as NCBI-RefSeq [26] for a curated set of type
strains and high-quality reference genomes. We conducted a benchmark of these four 16S
rRNA gene databases using RESCRIPt to compare various qualitative and quantitative characteristics and performance metrics. This benchmark was performed using full-length 16S rRNA
gene sequences from each database, so the goal was to compare full-length sequence and taxonomy information, not to simulate performance for commonly used short-read sequencing
technologies.
RESCRIPt’s evaluate-seqs action was used to examine sequence length distributions (Fig
2A), the number of unique sequences (Fig 2B), and sequence and kmer entropy as measures of
both richness and evenness of unique sequences in these databases (Fig 2C). The evaluate-taxonomy action was used to examine the number of unique taxonomic labels (Fig 3A), taxonomic entropy (Fig 3B), and the number of unclassified labels at each taxonomic rank (Fig
3C). The evaluate-fit-classifier and evaluate-classifications actions were used to compare optimal classification performance for each classifier (Fig 3D).
Results illustrate varying length distributions across databases, reflecting different proportions of Bacteria and Archaea in each, as well as different methods used to identify start and
end sites in each of these databases. Notably, some outliers (as short as 200 nt and as long as
3983 nt) were initially observed in SILVA, NCBI, and GTDB, presumably representing partial
and untrimmed 16S rRNA gene sequences. These were removed prior to downstream evaluation to avoid biasing performance metrics (see Methods section). Researchers should be aware
of length aberrations in these and other reference databases, and can use the evaluate-seqs
action in RESCRIPt to check length distributions in their own databases before proceeding.

Proportion

Sequence Length Distribution

Length (nt)

N Unique Sequences

300,000

silva_99

12

200,000
100,000
0

Entropy

13

H

Counts

400,000

greengenes_99
GTDB
NCBI-RefSeq

11
10

silva_99

greengenes_99

GTDB

NCBI-RefSeq

Sequence

32mer

16mer

8mer

Fig 2. Comparison of sequence information from SILVA, Greengenes, GTDB, and NCBI-RefSeq 16S rRNA gene databases. A, Sequence length distributions (after
removing outliers, see materials and methods). B, Number of unique sequences in each database. C, Entropy of full-length sequences and different kmer lengths in each
database.
https://doi.org/10.1371/journal.pcbi.1009581.g002

PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1009581 November 8, 2021

6 / 37

[с. 7]
PLOS COMPUTATIONAL BIOLOGY

Reproducible reference database management

Fig 3. Comparison of taxonomic information and simulated classification accuracy from SILVA, Greengenes, GTDB, and NCBI-RefSeq 16S rRNA gene
databases. A, Number of unique taxonomic labels; B, Taxonomic entropy; C, proportion of unclassified taxa at each rank; D, optimal classification accuracy (as
F-Measure) without cross-validation (simulating best possible classification accuracy when the true label is known but classification accuracy may be confounded by
other similar hits in the database). Cross-validation was not used because two of the databases (GTDB and NCBI-RefSeq) lack replicate species. Rank labels on x-axis:
D = domain, P = phylum, C = class, O = order, F = family, G = genus, S = species.
https://doi.org/10.1371/journal.pcbi.1009581.g003

The SILVA database exhibited the highest number of unique sequences (Fig 2B) and species
labels among the databases compared here (Fig 3A). However, ~72% of species labels present
in SILVA consist of unidentified, uncultured, or unknown organisms, and ~2.5% (excluding
chloroplast and mitochondrial sequences) do not match the genus label, leaving only ~25% of
sequences with meaningful species labels. Notably, this is because SILVA only curates the taxonomy to genus level but provides the “organism name” given to the sequence in the NCBI
GenBank source data, and hence genus–species mismatches can occur. Furthermore, GTDB
and Greengenes display similar levels of kmer entropy (Fig 2C), suggesting that these databases
cover a similar sequence space and taxonomic diversity to SILVA. GTDB actually has more
species-level annotations than SILVA (once unidentified and mismatched species labels are
discounted), and higher species label entropy (Fig 3B), indicating less redundancy. The lack of
species-rank curation in SILVA leads to poor optimal classification performance at the species
level (Fig 3D), yielding a species-level F-measure of 0.73, far below the other 16S rRNA gene
databases. By comparison, classification accuracy at the genus level is much higher for SILVA,
consistent with the level of curation performed (Fig 3D).
NCBI-RefSeq and GTDB share the highest species-level taxonomic entropy, and the most
unique species, after discounting the unknown/unmatched labels in SILVA (Fig 3A and 3B).
However, NCBI-RefSeq has fewer unique sequences and lower sequence entropy (Fig 2B and
2C), indicating that a lower amount of sequence space is covered, most likely because of the
stringent quality control and assessment process employed. Thus, NCBI-RefSeq exhibits high
quality reference sequences, but may have limited coverage for characterization of microbial
communities in some environments, e.g., where a large number of unknown species may be
encountered. In well characterized environments, this database is likely to offer competitive
advantages in terms of its curated taxonomy, size, and extensive use of type strains. NCBI-RefSeq exhibited the highest optimal classification accuracy (F = 0.94, Fig 3D), though this is likely
aided by the smaller database size and use of many genomes sequenced from type material,
reducing taxonomic ambiguities that are likely to occur in nature and are reflected in the other
databases. GTDB exhibited a slightly lower optimal classification accuracy (F = 0.92), indicating very high optimal accuracy in spite of its size, suggesting that the curation efforts and taxonomic re-classification strategies employed by the GTDB curators lead to a very well resolved
(if currently not officially recognized) taxonomic labeling scheme that closely aligns with the
16S rRNA gene sequence space.

PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1009581 November 8, 2021

7 / 37

