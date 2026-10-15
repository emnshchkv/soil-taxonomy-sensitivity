# Часть 8: с. 27–29

[с. 27]
PLOS COMPUTATIONAL BIOLOGY

Reproducible reference database management

1. Demonstrate various aspects of RESCRIPt’s current functionality for retrieving and curating reference sequences.
2. Compare the information content and classification performance of four commonly used
sequence databases for 16S rRNA gene classification of Bacteria and Archaea, retrieved and
formatted using RESCRIPt: SILVA [22], Greengenes [21, 23], NCBI-RefSeq [26,120,121],
and GTDB [25].
3. Evaluate the effects of sequence filtering on classification accuracy and information content
of the SILVA [22] rRNA gene sequence database.
4. Evaluate the effects of sequence clustering on sequence and taxonomic information content, using the Greengenes [23] 16S rRNA gene sequence database.
5. Evaluate the effects of sequence filtering on classification accuracy and information content
of the UNITE [31] ITS sequence database.
6. Evaluate the effects of sequence filtering on classification accuracy and information content
of the BOLD [34] and NCBI GenBank [64,121] COI gene sequence databases.
7. Evaluate the effects of sequence clustering and primer-coordinate trimming on classification accuracy and information content of the BOLD gene sequence database.
Data were retrieved either using RESCRIPt (for SILVA and NCBI data) or by direct download of release data (for UNITE, Greengenes, and GTDB) or by direct download (for BOLD
data; accessed July 1, 2020 and updated August 8, 2020).
SILVA data were filtered to remove sequences containing homopolymer lengths > 8 and/
or > 5 ambiguous characters, using the RESCRIPt action cull-seqs with default settings.
SILVA, GTDB, and NCBI data were all found to contain unusually long and short 16S rRNA
gene sequences (using the RESCRIPt action evaluate-seqs), thus Archaea sequences < 900 nt
[50] and Bacteria < 1200 nt [25], as performed for the SILVA releases (e.g. https://www.arbsilva.de/documentation/release-138/) were filtered out using the RESCRIPt action filter-seqslength-by-taxon.
The raw NCBI COI data were obtained with the ‘get-ncbi-data‘action in RESCRIPt (see
tutorials at https://github.com/bokulich-lab/RESCRIPt). Sequences containing the “BARCODE” keyword were considered cross-referenced NCBI data from BOLD, and labeled as
“NCBIob”, while those sequences without this keyword were labeled as “NCIBnb”. These data
were processed only in one fashion: initially dereplicated and primer trimmed to the ANML
primer coordinates, then dereplicated once again. This produced a pair of NCBI COI databases with (“ncbiOB”) or without (“ncbiNB”) the BOLD cross-referenced label. These data
were also combined, and again dereplicated, to represent the full NCBI set of COI sequences
(“ncbiAll”). In all cases, every database was filtered prior to running benchmark tests to retain
only those taxonomy labels and sequences associated with the phylum “Chordata” or
“Arthropoda”.
The raw BOLD COI data were obtained using a custom R script (https://osf.io/m5cgs/). All
data sets were filtered with similar RESCRIPt actions to remove sequences containing homopolymer lengths > 12 and/or > 5 ambiguous characters (via ‘cull-seqs’), and to retain
sequences between 250 to 1600 bp (‘filter-seqs-length’). Sequences were then either dereplicated (‘dereplicate’) or clustered using one of three percent identities (97, 98, 99). These “boldFull” datasets contained sequences of variable lengths of COI sequence, and were further
trimmed to the boundaries that map to ANML [122] sequences—a primer pair commonly
used in animal diet metabarcoding experiments. Trimming was performed using MAFFT

PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1009581 November 8, 2021

27 / 37

[с. 28]
PLOS COMPUTATIONAL BIOLOGY

Reproducible reference database management

[123] in three stages: first, a subset of reference sequences created a high quality alignment
with ‘mafft—auto’; second, primers were aligned to this small alignment file with ‘mafft—multipair—addfragments‘; third, the remaining reference sequences were aligned with ‘mafft—
auto—addfull’. These trimmed data were dereplicated once more to produce the “boldANML”
datasets.
In several benchmarks, reference sequences and taxonomy were dereplicated and/or clustered using the RESCRIPt action ‘dereplicate’. This action uses VSEARCH to dereplicate
sequences and optionally cluster them at a specified % similarity to form operational taxonomic units (OTUs), then RESCRIPt finds the most appropriate taxonomic label for each
sequence cluster using one of several available modes of operation to find the last common
ancestor (LCA) for the cluster, or to preserve identical sequences with unique taxonomic
labels.
Taxonomic information in each database was evaluated using the RESCRIPt action evaluate-taxonomy. This action measures the number of unique labels, label entropy, and the number of unknown/unclassified labels at each taxonomic rank. Shannon’s entropy (H) [79] is
defined as:
n
X

HðXÞ ¼

Pðxi Þ loge Pðxi Þ
i¼1

Thus, entropy relates to both the evenness and richness of information content: e.g., the
number and evenness of taxonomic label frequencies or sequence/kmer frequencies.
Sequence information in each database was evaluated using the RESCRIPt action ‘evaluateseqs’. This action measures the number of unique sequences, sequence entropy, sequence
length distribution, and kmer entropy.
Taxonomic classification accuracy was simulated for each database using the RESCRIPt
actions ‘evaluate-cross-validate’ and ‘evaluate-fit-classifier’, followed by accuracy evaluation
with evaluate-classifications (which measures precision, recall, and F-measure in taxonomic
classification results [115] and visualizes these metrics at each taxonomic rank). RESCRIPt
implements two different classification evaluation methods to simulate different classification
conditions, including both “unrealistic” (easy) and “realistic” classification tasks, i.e., that
many of the sequences detected in a true environmental survey, e.g., of microbial or eDNA
sequences in most sample types, will not have exact matches in any reference database, either
because they represent novel strains, species, or higher-order taxonomic clades. The ‘evaluatecross-validate’ action uses k-fold cross-validation (implemented in scikit-learn [114]) to perform a pseudo-realistic classification task whereby a set of query sequences may not have an
exact match in the reference database, but other similar taxonomic groups may be present, as
implemented and described previously [56,115]. This action splits a database into K test sets
(such that each sequence appears in a test set exactly once) and classification is performed in
each fold with the remaining sequences as the training set. Splitting is stratified by taxonomic
groups, so that taxonomic groups are evenly stratified across training and test sets. An
“expected” taxonomy is generated by this method, in which taxonomic singletons (i.e.,
sequence queries that do not have a representative from the same taxonomic group in the
training set) are truncated so that the expected taxonomic classification is the LCA between
that taxon and its nearest taxon in the training set, as this is the “correct” answer when the true
taxonomic label for that query is not “known” (i.e., absent from the training set). The ‘evaluate-fit-classifier’ action trains and tests classification on the full dataset, without cross-validation, to report the best-case performance, i.e., when each query sequence has an exact match in
the reference database (and hence the correct taxonomic label is known, but other matches

PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1009581 November 8, 2021

28 / 37

[с. 29]
PLOS COMPUTATIONAL BIOLOGY

Reproducible reference database management

may also be present). In this case, data leakage (where information is shared between the test
and training sets) is intentional, in order to estimate the upper bound of classification accuracy
for a given database by simulating an unrealistically easy classification task (using the definition of “realistic” described above). Both of these classification evaluation methods can be
adapted in RESCRIPt to simulate the expected level of challenge in a given ecosystem—e.g.,
for well characterized sample types and clades, ‘evaluate-fit-classifier’ method may actually
represent a realistic classification scenario, and users can set different levels of K for ‘evaluatecross-validate’ (to adjust the number of splits performed) to adjust the degree of “challenge”
(i.e., lower levels of K result in larger splits and more uncertainty).
In addition to dereplicating and clustering sequences, the RESCRIPt ‘dereplicate’ action is
useful as a quick assessment of taxonomic resolution, versus using the classification simulation
methods described above. By dereplicating the sequences and using the ‘lca’ mode to find the
LCA for each cluster, followed by using the ‘evaluate-taxonomy’ visualizer (described above)
to examine the number of unclassified labels at each rank, this action allows a quick assessment
of how well the different taxonomic groups contained in the database can be resolved based on
sequence information alone, as we demonstrate in some of our benchmarks. The classification
simulation/evaluation methods are better suited for simulating realistic classification tasks,
and for estimating actual taxonomy classifier performance, but are computationally intensive
and time-consuming.

Reproducible genomics workflows
To demonstrate the ability of RESCRIPt to process and compare genome data, we used the
scalable MinHash (MASH) approach [108] through q2-sourmash (https://github.com/dib-lab/
q2-sourmash) [88]. In brief, MASH generates compressed sketch representations of large
genome sequence sets, making large genome comparisons possible through dimensionality
reduction. We queried the NCBI Virus portal [124], for the Hepatitis E Virus (HEV) on September 23rd, 2020 and downloaded accessions within the viral lineage "Hepatitis E virus
taxid:12461", that were annotated as having a nucleotide completeness status of "complete".
Only HEV data with sufficient regional representation, and associated with humans, were
retained. These accessions were downloaded using RESCRIPt’s “get-ncbi-data” function, and
processed using a modified version of q2-sourmash (https://github.com/mikerobeson/
q2-sourmash/tree/use-fasta). The q2-sourmash functions “compute-fasta” and “compare”
were used to generate MASH signatures for each genome and perform pairwise genome comparisons respectively. PCoA was performed through the QIIME 2 q2-diversity plugin and visualized with EMPeror [90]. Finally, q2-sample-classifier [89], was used to determine if MASH
signatures are predictive of geographic source, based on k-nearest-neighbors classification
with leave-one-out cross-validation.

Supporting information
S1 Text. Earth Microbiome Project reference database comparison summary. Methods and
results summary on how database preparation affects taxonomic classification of real-world
biological data from the Earth Microbiome Project.
(DOCX)
S1 Fig. Comparison of reference database classification of Earth Microbiome Project
sequences. Comparison of taxonomic classification of Earth Microbiome Project sequences
with SILVA, Greengenes, GTDB, and NCBI-RefSeq 16S rRNA gene databases. A, Number of
unique taxonomic labels; B, Taxonomic entropy; and C, proportion of unclassified taxa at each

PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1009581 November 8, 2021

29 / 37

