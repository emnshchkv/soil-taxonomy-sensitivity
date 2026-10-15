# Часть 7: с. 24–26

[с. 24]
PLOS COMPUTATIONAL BIOLOGY

Reproducible reference database management

the quality of its inputs and the processing decisions made by the user. In general, users should
use multiple metrics to guide their interpretations of RESCRIPt’s results, but also need to be
aware of the composition of input data before making conclusions about database quality.
For example, classification accuracy metrics output by RESCRIPt can be artificially high if
the input database is of low quality. This is clearly seen in the OTU clustering benchmark performed using the Greengenes database (Fig 8D and 8E): the “evaluate-fit-classifier” method,
which classifies sequences without cross-validation, reports perfect and near-perfect specieslevel classification accuracy on the highly clustered sequences (e.g., 64 and 79% OTUs) (Fig
8D). This is, however, because these sequences are clustered to the extent that the remaining
taxonomic coverage becomes relatively poor and sparse as near-neighbors become clustered
into fewer OTUs (Fig 8A). This sparsity becomes reflected in the poor classification accuracy
when cross-validation is used, as the lack of near neighbors leads to high misclassification rates
for the 64% and 79% OTUs (Fig 8E). Using multiple classification evaluation methods and
both qualitative and quantitative metrics (e.g., comparing classification accuracy to taxonomic
and sequence entropy and coverage information) will help guide more robust conclusions
about database quality.

Future goals
RESCRIPt currently contains a range of tools for sequence reference database acquisition,
management, and evaluation. Curation tools remain manual and qualitative. In the future, we
plan to explore and incorporate methods for quantitative curation of sequences and taxonomies. For example, the methods used by GTDB [24] to inform species clusters based on average nucleotide identity (ANI) [106] could be incorporated for similar curation of species
clusters in RESCRIPt. Other methods to detect and re-annotate mis-annotated and unannotated sequences would be valuable for guiding sequence curation efforts. However, methods
such as these can be difficult to apply generally, and while ANI is useful for defining species
clusters based on whole-genome sequences, it may not scale appropriately to incomplete
genomes or marker-gene sequences [107].
Although the scope and benchmarks included in this study focus on marker-gene sequencing applications, genome and metagenome databases are already compatible with RESCRIPt.
For example, the “get-ncbi-data” method could be used to automatically download reference
genomes from GenBank, and the filtering and taxonomy functions are general purpose. More
genome- and metagenome-focused functionality is planned for future releases of RESCRIPt,
such as ANI [106] and MASH [108] for (meta)genome distance estimation, and methods for
estimating the taxonomic classification accuracy of (meta)genome databases.
RESCRIPt’s developers remain committed to working with researchers to provide access to
leading reference materials and reproducible and transparent reference sequence processing
workflows. We plan to add more methods for sequence and taxonomy acquisition from public
online databases that are commonly used by the marker-gene and metagenome research community, and welcome collaboration with database curators who want to better integrate their
databases with RESCRIPt. The most up-to-date information related to feature requests, usage,
and troubleshooting can be found on the project’s GitHub page (https://github.com/bokulichlab/RESCRIPt).

Methods
Implementation
RESCRIPt is implemented as a free, open-source QIIME 2 [72] plugin, in order to leverage
QIIME 2’s data provenance tracking system to ensure that users can trace the steps used to

PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1009581 November 8, 2021

24 / 37

[с. 25]
PLOS COMPUTATIONAL BIOLOGY

Reproducible reference database management

make their custom reference databases. QIIME 2 “Artifacts” (results files) consist of zip
archives containing the result file (in typical, interoperable formats, e.g., FASTA for nucleotide
sequence data) as well as data provenance information and other file metadata. Users unfamiliar or unwilling to work with QIIME 2 Artifacts can extract the data using either QIIME 2 or
standard methods (e.g., the UnZip utility (http://infozip.sourceforge.net/) that is included in
most Linux and Unix distributions), making QIIME 2 Artifacts a fully interoperable, portable
solution for storing reference databases with integrated provenance information.
RESCRIPt is written primarily in Python 3 and depends on pandas [109,110] for dataframe
operations; VSEARCH [111] and scikit-bio (scikit-bio.org) for parsing nucleotide sequences;
numpy [112], scipy [113], and scikit-learn [114] for numerical and statistical operations;
xmltodict (https://github.com/martinblech/xmltodict) for parsing XML; urllib and requests
(https://github.com/psf/requests) for HTTP requests; q2-feature-classifier [115] for taxonomic
classification; and matplotlib [116], seaborn [117], VEGA [118], and q2-longitudinal [119] for
plotting and data visualization, including interactive visualizations.
The current release of RESCRIPt (2021.8) contains a variety of functions for retrieving,
managing, and evaluating sequence and taxonomy reference databases (Fig 1). Details on specific functions, usage, and tutorials can be found at the project website (https://github.com/
bokulich-lab/RESCRIPt), and are described in the sections below.
SILVA data retrieval and taxonomy formatting. RESCIPt supports retrieval of SSU and
LSU marker-gene data from SILVA via an automated method, “get-silva-data”, or manual
import of the necessary sequence and taxonomy files (Fig 1). The “get-silva-data” pipeline
allows selection of (a) which version of the database to download, (b) whether to download
LSU, SSU sequences, or the SSU NR99 sequences, and (c) which taxonomic ranks to use and
other options for taxonomy parsing (see software documentation for more details). These
options are all stored in the data provenance of the output files, for later retrieval and reproducibility. RESCRIPt parses the SILVA taxonomy, using three files as input:
• taxonomy rank (taxrank) file, containing both the taxonomic rank and taxonomy for each
numeric taxonomy identifier (taxid);
• taxonomy mapping (taxmap) file, which maps each sequence accession to a taxid and the
“organism name” provided by NCBI;
• taxonomy tree (taxtree) file, which contains the hierarchical taxonomy in Newick format,
with the taxids used as node labels, in which the daughter nodes contain the taxids of lowerlevel taxonomies.
The “parse-silva-taxonomy” method utilizes the taxrank, taxmap, and taxtree files to
generate a consistent user-defined rank-associated taxonomy. Although the set of ranks
can be configured by the user, the following ranks are extracted by default: domain (d_),
phylum (p_), class (c_), order (o_), family (f_), and genus (g_). Any ranks not associated
with taxonomy have their upper-level taxonomic lineage propagated downward (i.e. the
values are forward filled with the last observed taxonomic value) towards lower-level ranks.
This ensures general compatibility with downstream taxonomy classification tools, many
of which may require non-empty fields at each rank. Rank propagation can be optionally
disabled.
Finally, the user can choose to append the organism name (from the taxmap file) for use as
the species (s_) rank taxonomy. We generally warn against this due to the myriad of inconsistent information found within the organism name field (based on our benchmarking results
described herein), but it can occasionally be useful. If the user does decide to leverage the
organism name, we currently only return the first two words, to remove subspecies-level

PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1009581 November 8, 2021

25 / 37

[с. 26]
PLOS COMPUTATIONAL BIOLOGY

Reproducible reference database management

information that is often included in the given organism name and which can degrade classification accuracy (e.g., because the extra information causes that species to be interpreted as a
unique label).
Rank propagation is provided to allow users to extract more taxonomic information, rather
than explicitly pulling down only the ranks of interest. For example, if a user opted to download sequence data along with only the six standard taxonomic ranks (see above), they may
obtain the following taxonomic output when rank propagation is not used:
Z27393.1.1722 d__Eukaryota; k__Fungi; p__Ascomycota; c__; o__; f__; g__
AB671439.1.2071 d__Eukaryota; k__Fungi; p__Ascomycota; c__; o__; f__; g__
The user might assume that query sequences that “hit” either of these reference sequences
would be unable to classify beyond the phylum level. However, applying rank propagation will
yield the following for these same accessions:
Z27393.1.1722 d__Eukaryota; k__Fungi; p__Ascomycota; c__Taphrinomycotina;
o__Taphrinomycotina; f__Taphrinomycotina; g__Taphrinomycotina
AB671439.1.2071 d__Eukaryota; k__Fungi; p__Ascomycota; c__Pezizomycotina; o__Pezizomycotina; f__Pezizomycotina; g__Pezizomycotina
This is because intermediate ranks not selected by the user (e.g., sub-phyla Taphrinomycotina and Pezizomycotina) were propagated downward and used to fill in the unannotated
ranks. Hence, forward filling allows users to disambiguate incompletely annotated reference
sequences. The drawback is the conflation of taxonomy by mixing ranks from other levels.
The RESCRIPt project page (https://github.com/bokulich-lab/RESCRIPt) lists several tutorials describing how to use various RESCRIPt functions, including methods to import and
parse SILVA data.
NCBI GenBank data retrieval and taxonomy formatting. RESCRIPt supports automated retrieval of sequence taxonomy databases from the NCBI Nucleotide and Taxonomy
databases [85] using the “get-ncbi-data” method. Sequences can be selected using a standard
NCBI query, by specifying a list of sequence accession ids, or as a combination of the two.
Downloads of large sequence databases can be made faster using parallel connections and
batch downloads.
The NCBI download method retrieves the requested sequences from the NCBI Nucleotide
database, cross-references their taxids with the NCBI Taxonomy database to obtain their taxonomic classifications, then standardizes the taxonomies to adhere to a fixed set of ranks. The
set of ranks can be configured by the user but are kingdom (k__), phylum (p_), class (c_),
order (o_), family (f_), genus (g_), and species (s__) by default. If a given rank is not present in
the NCBI Taxonomy it is propagated down from the nearest higher rank, as described above.
Rank propagation can be optionally disabled.
The RESCRIPt project page (https://github.com/bokulich-lab/RESCRIPt) lists several tutorials describing how to use various RESCRIPt functions, including methods to download and
save NCBI data.

Reference database benchmarks
To demonstrate some examples of how users can use RESCRIPt to process and evaluate custom reference databases from popular source data, we used RESCRIPt to conduct several
benchmarks. Tutorials demonstrating this functionality, based on some of the following
benchmarks, can be found on the project website (https://github.com/bokulich-lab/
RESCRIPt). Workflows and data from our benchmarks can be found at https://github.com/
bokulich-lab/db-benchmarks-2020 and https://github.com/devonorourke/COIdatabases.
Benchmarks were designed and executed with the following aims:

PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1009581 November 8, 2021

26 / 37

