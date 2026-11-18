# Часть 6: с. 21–23

[с. 21]
PLOS COMPUTATIONAL BIOLOGY

Reproducible reference database management

The reproducibility problem
The transparent reporting, replication, and reproduction of scientific discoveries is not a new
problem, but it is one that has become complicated in the digital age, as experiments and
computational methods become increasingly sophisticated and datasets have become both
larger and more likely to be re-used by others [65,67,68,91]. Reference database selection, and
any subsequent curation, critically impact the findings from marker-gene and metagenome
experiments [44,77,80,92], and hence must be carefully documented to allow others to interpret scientific findings. For example, inappropriate use of environment-specific databases have
been shown to yield alarming false-positive rates in metagenomics datasets [80]. As reference
data are circulated and re-used by other researchers, it is critical that database curation steps
be transparently documented and transmitted so that the impact of these decisions can be evaluated downstream both by researchers re-using those datasets, as well as by the wider scientific
community when interpreting results. To address this issue, RESCRIPt utilizes QIIME 2’s integrated data provenance tracking system [72] to record and store processing steps inside each
individual file generated as part of a workflow. Hence, provenance can be retrieved from any
terminal result file to document and replicate the entire processing chain, from data acquisition to filtering to downstream use (e.g., for taxonomic classification or sequence analysis).
We believe that embedding provenance within reference database files should become a standard in the field whenever reference data are modified by a researcher or destined for re-use
by others, and hope that others utilize RESCRIPt to facilitate greater transparency and reproducibility within the microbiome, eDNA, and metagenomics communities.

The curation problem
Our results have shown that formatting and correcting taxonomy and other metadata are critical components of reference database generation and management prior to applications for
sequence classification, e.g., to standardize taxonomic ranks across entries [93]. We have
already implemented functions in RESCRIPt to format the popular SILVA rRNA gene and
NCBI GenBank databases, and are planning future support for parsing and editing other taxonomy formats, as well as mapping between these formats [73].
There are four codes of nomenclature as reviewed in [94], the International Code of
Nomenclature for algae, fungi and plants (ICNafp; [95]) International Code of Nomenclature
of Prokaryotes (ICNP; [96]), International Code of Zoological Nomenclature (ICZN; [97]),
International Code of Virus Classification and Nomenclature (ICTV, [98]). Each of these have
their own rules for taxonomic curation within their respective areas. Combining and curating
taxonomic information across multiple databases can be an onerous task as new issues can
arise when attempting to merge information across the respective authorities on nomenclature
[94]. For example, not all taxonomic ranks are recognized, available, or even treated equally
across the various databases. Valid rank fields in one database may not be recognized, or even
useful, in other databases, e.g. INCafp formally recognizes the below-species ranks varieties
(varietas) and forms (forma), which are not recognised within the other codes of nomenclature. Inconsistencies in taxonomic labeling, updates, and rank suffixes can cause additional
incompatibilities between databases. Current efforts seek to standardize higher level suffixes of
microbial nomenclature [99], streamlining bioinformatic extraction and inference of rank
information.
In recent years, the explosion of high-throughput sequencing technologies has allowed
researchers to generate genomic data on many as yet uncultured microbial taxa. In fact, the
rate at which novel genomic data can be acquired [86], and rapidly placed within a phylogenetic context [24], has surpassed our ability to appropriately resolve any conflicts with

PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1009581 November 8, 2021

21 / 37

[с. 22]
PLOS COMPUTATIONAL BIOLOGY

Reproducible reference database management

traditional Linnaean taxonomy. This has resulted in some proposals on how to taxonomically
organize genomic data from uncultured microbes, and with greater emphasis on phylogenetic
systematics [24,100]. For example, commonly used labels may have no officially recognized
rank (e.g. Opisthokonta), but are quite informative, as they may refer to monophyletic grouping of taxa.
Sequences with missing or incomplete taxonomy/metadata and low-quality sequences are
common in some databases. RESCRIPt introduces easy-to-use and user-customizable tools for
detection and removing low-quality entries based on sequence filtering criteria (e.g., presence
of homopolymer or ambiguous bases) or based on taxonomic information. Although some of
these are trivial tasks for experienced bioinformaticists, using RESCRIPt to perform these
functions results in those processing decisions being recorded automatically in the provenance
stored in the downstream results files, so that this information can be recovered at any stage of
downstream processing (provided that the data are maintained in a QIIME 2 archive format).
The aim of RESCRIPt is to democratize the tools for database acquisition, formatting, and
curation, but RESCRIPt is not in itself a tool to automatically curate data. It is the responsibility of users to check the validity of their source data and to use those reference data and
RESCRIPt appropriately. Even popular sequence taxonomy databases are prone to error [43],
and issues with taxonomic naming, polyphyly, and inconsistent degrees of sequence similarity
at different taxonomic ranks can complicate the use, accuracy, and contemporaneity of reference databases [101]. Defining definitions of taxonomic boundaries has historically been a
challenge and can vary based on the characters used to define them, e.g. biochemical, metabolomic, ecological phenotypes, with more recent definitions, and circumscription of taxa, based
upon phylogenetic relatedness and average nucleotide identity. However, these new phylogenetic approaches have resulted in either the lumping or splitting of taxa creating ever more
inconsistencies between taxonomy and phylogeny [16,25,102,103].
We hope that by decreasing the number of technical hurdles involved with the generation
and curation of custom databases, we will ease the point of entry, and create an interest in the
taxonomic sciences among the research community. In the future, RESCRIPt could help facilitate data curation in large-scale citizen science projects in which the greater research community and general public can contribute to the growth and curation of sequence and taxonomy
databases.

The evaluation problem
After constructing a custom database comes the critical question: is the database actually useful
and informative? How does it compare to other databases? User-friendly methods for
sequence reference database evaluation are not currently available, making database evaluation
and benchmarking a formidable challenge to the research community. RESCRIPt has implemented multiple methods for database evaluation, which generate interactive visualizations to
allow users to explore and better interpret database quality characteristics (see example gallery
at https://github.com/bokulich-lab/RESCRIPt). These involve both qualitative metrics, for
evaluating sequence and taxonomy information within and between databases, as well as
quantitative metrics for evaluating taxonomic classification performance of marker-gene
sequences with different cross-validation schemes. Furthermore, we provide reproducible
examples in the online tutorials (https://github.com/bokulich-lab/RESCRIPt) to guide users in
the use and interpretation of these methods, to make these methods widely available and
usable by the research community.
Using these evaluation methods, we have benchmarked performance characteristics of
some of the most popular reference databases for bacterial/archaeal 16S rRNA gene, eukaryote

PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1009581 November 8, 2021

22 / 37

[с. 23]
PLOS COMPUTATIONAL BIOLOGY

Reproducible reference database management

ITS region, and animal COI gene sequences. This evaluation informs several conclusions.
First, all of these databases may require some additional curation by end-users to improve suitability for certain research applications. This includes filtering low-quality sequences and
annotations to improve database quality and classification accuracy, such as abnormally short
sequences in the GTDB, SILVA, and NCBI-RefSeq databases. Second, we compared several of
these databases side-by-side to measure relative performance metrics. In the case of 16S rRNA
gene analysis specifically, we conclude that the size and taxonomic comprehensiveness of
SILVA are major assets, though GTDB and NCBI-RefSeq may be more suitable for various
applications that respectively require greater taxonomic and phylogenetic rigor. The use of
genomes sequenced from type material provides these two databases with a robust taxonomic
and phylogenetic backbone that enables users to link natural history and experimental science
[94,104].
NCBI-RefSeq’s species records are extracted from data submissions to the International
Nucleotide Sequence Database Collaboration (INSDC), i.e., NCBI-GenBank, the European
Nucleotide Archive (ENA), and the DNA Data Bank of Japan (DDBJ). NCBI-RefSeq continually curates these species records from the primary data, often by collaborating with other
groups, e.g. authorities in sequence data curation, taxonomic nomenclature, phylogenetic systematics, et cetera. Furthermore, NCBI-Taxonomy continually runs taxonomic consistency
checks on assembled genomes with average nucleotide identity (ANI) [105]. These curational
efforts result in a well integrated suite of biological information that can be interrogated
through a variety of means and data types [26,94]. The GTDB extracts and curates data from
both NCBI-RefSeq and NCBI-GenBank to generate a phylogeny of Archaea and Bacteria from
roughly 120 ubiquitous single-copy proteins [24,93]. This phylogeny is used to inform microbial taxonomy, especially in cases where a given taxonomy is observed to be polyphyletic. In
this case, a conservative approach is used to remove polyphyletic groups and normalize taxonomic ranks according to their relative evolutionary divergence. Any remaining polyphyletic
groups are then flagged as “regions of instability” in the hopes that future in-depth analyses
will result in a stable set of classifiable taxa [24,93]. The efforts by NCBI and GTDB enable
researchers to not only more accurately classify uncultured microbes, but also place them into
ecological and evolutionary context based on their nearest phylogenetic neighbors.
The design, curation decisions, and ultimate quality of these databases must be considered
carefully when applying them for particular purposes. For example, the SILVA database is not
curated at species level, though the “organism name” provided in the source NCBI data are
provided. RESCRIPt’s “get-silva-data” method for acquiring and formatting SILVA data can
be configured to either report these organism names as the species labels in the output hierarchical taxonomy annotations, or to only report the SILVA taxonomy, which is curated from
the domain to genus rank. In our evaluation, we found that although SILVA species labels can
be informative, 72% consist of unidentified, uncultured, or unknown organisms, and 2.5% do
not match the genus. Downstream users should be aware of such caveats when using SILVA,
or any reference data, and understand the limitations that these impose when interpreting
results.

What RESCRIPt does not do, and other limitations
RESCRIPt is designed to give researchers access to tools for reproducible nucleotide sequence
and taxonomy database generation, curation, and evaluation. RESCRIPt is not in itself a data
source nor an authority on taxonomy, systematics, or data quality, and the qualitative and
quantitative metrics that RESCRIPt can generate are not infallible indicators of quality or accuracy. As with any bioinformatics methods, the quality of RESCRIPt’s outputs is dependent on

PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1009581 November 8, 2021

23 / 37

