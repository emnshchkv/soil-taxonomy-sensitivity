# Часть 5: с. 16–20

[с. 16]
PLOS COMPUTATIONAL BIOLOGY

Reproducible reference database management

Fig 10. Comparison of sequence information from BOLD COI gene database for available arthropod and chordate sequences. Differences in
datasets reflect whether sequences were trimmed to a particular primer region (boldANML) or not (boldFull), and whether sequences were dereplicated
(100) or clustered at a particular percent identity (97, 98, 99). A, Number of unique sequences. B, Entropy of sequences and different kmer lengths.
https://doi.org/10.1371/journal.pcbi.1009581.g010

reduction in the number of terminal labels (Fig 11C). Although both arthropod and chordate
references contain the most sequences with species-rank labels in BOLD, arthropods uniquely
contain more sequences ending with family-level information than genus-level, which is perhaps an indication of the greater challenge in classifying many arthropod specimens. The
entropy results, paired with the number of terminal labels at a given rank, indicate that
sequence clustering and primer trimming are both more consequential with respect to reducing the total number of sequences (richness) than the redundancy of labels (evenness).
While trimming and clustering followed similar trends with respect to taxonomic information, the results of the evaluate-fit-classifier, our “best-case” classification accuracy (Fig 11D),
indicate opposing outcomes with respect to these two processes: primer trimming reduces true
taxonomic classification accuracy (due to reduced sequence information and thus lowered
ability to distinguish taxa), while clustering artificially increases it (by reducing genetic complexity and clustering genetically similar clades). Thus, the “best-case” accuracy for these COI
sequences was obtained when references were untrimmed and clustered at 97%. Notably, the
magnitude of these effects varied by taxonomic group: chordate sequences were more sensitive
to clustering than arthropods, and primer trimming was more impactful for arthropods than

PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1009581 November 8, 2021

16 / 37

[с. 17]
PLOS COMPUTATIONAL BIOLOGY

Reproducible reference database management

Fig 11. Comparison of taxonomic information and simulated classification accuracy from BOLD COI gene database for available arthropod and chordate
sequences. Differences in datasets reflect whether sequences were trimmed to a particular primer region (boldANML) or not (boldFull), and whether sequences were
dereplicated (_100) or clustered at a particular percent identity (_97, _98, _99). A, Number of unique taxonomic labels; B, Taxonomic entropy; C, proportion of
unclassified taxa at each rank; D, optimal classification accuracy (as F-Measure) without cross-validation (simulating best possible classification accuracy when the true
label is known but classification accuracy may be confounded by other similar hits in the database). E, Classification accuracy with cross-validation. Rank labels on xaxis: K = kingdom, P = phylum, C = class, O = order, F = family, G = genus, S = species.
https://doi.org/10.1371/journal.pcbi.1009581.g011

chordates. Cross-validated classification suggests a more unified pattern with respect to classification accuracy, such that trimming and clustering both reduce accuracy (Fig 11E). As mentioned previously with regards to clustering the Greengenes SSU OTUs, the differences
between “best-case” and cross-validated accuracy are likely driven by a loss of information
with increasing OTU clustering. Collectively, our data suggest that OTU clustering is detrimental for COI gene classification (Fig 11).

Comparison of metazoan COI gene sequences in BOLD and GenBank
Next, we compared dereplicated and primer-trimmed metazoan COI reference sequences
obtained from either BOLD or NCBI GenBank. Because some COI reference sequences
have been deposited in both BOLD and NCBI, the NCBI data were partitioned into
sequences cross-referenced to BOLD (“ncbiOB”) or not (“ncbiNB”), as well as represented
in its totality (“ncbiAll”). In addition, we evaluated the number of distinct taxonomic labels
shared among these databases to illustrate the degree of taxonomic information shared
between groups.
Dereplicated and primer-trimmed sequences obtained from BOLD contained a slightly
larger number of unique arthropod and chordate references compared to those obtained

PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1009581 November 8, 2021

17 / 37

[с. 18]
PLOS COMPUTATIONAL BIOLOGY

Reproducible reference database management

Dataset

Fig 12. Comparison of sequence information from BOLD and NCBI GenBank COI gene databases for available arthropod and chordate
sequences. All sequences were dereplicated and trimmed to a common primer region. NCBI references either contained a cross-reference term to
BOLD (“ncbiOB”) or not (“ncbiNB”) or were combined together (“ncbiAll”). A, Number of unique sequences (note difference in scales between
Arthropoda and Chordata). B, Entropy of sequences and different kmer lengths.
https://doi.org/10.1371/journal.pcbi.1009581.g012

through NCBI GenBank (Fig 12A). Although many thousands of sequence accessions are
cross-listed in both NCBI GenBank and BOLD, these data indicate that tens of thousands of
COI sequences publicly available via BOLD are not cross-listed in NCBI GenBank. Sequence
entropy was similar among all combined NCBI datasets and BOLD (Fig 12B) for both chordate
and arthropod references.
Despite having fewer overall unique sequences, data obtained from NCBI contained
more unique genus and species labels than BOLD arthropod references, and a similar number of BOLD chordate references (Fig 13A), leading to higher taxonomic entropy at the
genus and species levels for NCBI (Fig 13B). Similarly, the combined NCBI database
(“ncbiAll”) contains many more arthropod and slightly more chordate sequences that are
assigned species labels after truncation and dereplication (Fig 13C). For arthropod references, we find that BOLD data is the least accurate database among the “best case scenario”
classification (Fig 13D), but performs the best when subject to cross validation (Fig 13E).
However, chordate references are consistent between both measures, indicating that NCBI
references provide improved accuracy relative to BOLD references from family through
species-levels (Fig 13D and 13E).

PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1009581 November 8, 2021

18 / 37

[с. 19]
PLOS COMPUTATIONAL BIOLOGY

Reproducible reference database management

Fig 13. Comparison of taxonomic information and simulated classification accuracy from BOLD and NCBI GenBank COI gene databases for available
arthropod and chordate sequences. All sequences were dereplicated and trimmed to a common primer region. NCBI references either contained a cross-reference
term to BOLD (“ncbiOB”) or not (“ncbiNB”) or were combined together (“ncbiAll”). A, Number of unique taxonomic labels; B, Taxonomic entropy; C, proportion of
unclassified taxa at each rank; D, optimal classification accuracy (as F-Measure) without cross-validation (simulating best possible classification accuracy when the true
label is known but classification accuracy may be confounded by other similar hits in the database). E, Classification accuracy with cross-validation. Rank labels on xaxis: K = kingdom, P = phylum, C = class, O = order, F = family, G = genus, S = species.
https://doi.org/10.1371/journal.pcbi.1009581.g013

Fetching reference genomes for classification
The generation and use of genomic data is increasing [86], and allows better phylogenetic and
taxonomic resolution of microorganisms [87]. Thus the need to easily and reproducibly
retrieve and use this information for either direct analyses or generating custom metagenomics
reference databases has become imperative. To address this need, RESCRIPt supports the
automated retrieval of reference genomes from NCBI GenBank, enabling extensible and
reproducible genomics workflows via interaction with other QIIME 2 plugins, and registering
this pipeline in data provenance for greater transparency and reproducibility downstream (Fig
14A). To highlight this functionality, we used RESCRIPt along with several other plugins,
q2-sourmash (https://github.com/dib-lab/q2-sourmash) [88], q2-sample-classifier [89],
q2-diversity, and EMPeror [90] to generate a reproducible workflow to acquire and process a
set Hepatitis E virus (HEV) genomes from NCBI-GenBank (Fig 14A), generate MASH signatures, perform pairwise genome comparisons, and visualize genome similarity via PCoA (Fig
14B). Finally, we show that HEV genotype (MASH signature) is predictive of geographic
source (Accuracy = 86.1%) using k-nearest-neighbors classification with leave-one-out crossvalidation (Fig 14C). Note that this analysis was performed using a subset of HEV genomes

PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1009581 November 8, 2021

19 / 37

[с. 20]
PLOS COMPUTATIONAL BIOLOGY

Reproducible reference database management

Fig 14. An example of using RESCRIPt for reproducible genomics workflows. HEV genomes were downloaded from NCBI-GenBank and used to make a reference
genome classifier based on the following geographic locations: Bangladesh (BD), China (CN), France (FR), India (IN), and the United Kingdom (UK). The
interoperability of RESCRIPt with other QIIME 2 plugins enables users to chain together a variety of functions into fully reproducible workflows that record processing
decisions in data provenance. A, a simplified data provenance graph highlighting our workflow leveraging RESCRIPt, q2-sourmash, q2-diversity, q2-sample-classifier,
and EMPeror. B, PCoA plot of individual HEV genomes based on MASH signature comparison results. C, k-nearest-neighbor classification accuracy based on MASH
signature dissimilarities and geographic location.
https://doi.org/10.1371/journal.pcbi.1009581.g014

merely as a demonstration of RESCRIPt’s broad functionality, and does not represent a fully
structured test from which biological conclusions should be drawn.

Discussion
The acquisition problem
Curated reference materials are publicly available for commonly used marker-genes such as
rRNA genes [21,22,24,25,33] and the ITS region [31] for various domains of life. However,
public curated reference databases are currently lacking for many other marker genes and for
particular clades, creating a major bottleneck in scientific research. Even when such databases
do exist, acquiring and formatting these data for use with standard methods for sequence analysis and taxonomy classification can present a steep learning curve for scientists and clinicians
who lack the bioinformatics expertise required to generate and manage custom sequence and
taxonomy databases. RESCRIPt resolves many of these issues, providing automatic tools for
generating and formatting custom sequence and taxonomy databases (from either marker
genes or genomes) from NCBI GenBank and from SILVA. Methods to provide similar support
for other commonly used databases are planned for future releases of RESCRIPt. We hope that
these methods will democratize the process of generating custom reference databases, supporting research efforts across the microbiome, eDNA, and metagenomics communities.

PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1009581 November 8, 2021

20 / 37

