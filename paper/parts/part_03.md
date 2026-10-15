# Часть 3: с. 8–11

[с. 8]
PLOS COMPUTATIONAL BIOLOGY

Reproducible reference database management

Fig 4. Comparison of taxonomic coverage among SILVA, Greengenes, GTDB, and NCBI-RefSeq 16S rRNA gene databases. Each panel displays the proportion of
taxa represented in one reference database (as indicated in the panel title) that are shared with each other database at each taxonomic rank. The legends indicate which
groups are being compared: the reference alone (always 1.0, shown for clarity), pairs consisting of the reference and one other database, trios consisting of the reference
and two other databases, and the proportion of the reference’s labels that are shared by all four databases. Rank labels on x-axis: D = domain, P = phylum, C = class,
O = order, F = family, G = genus, S = species.
https://doi.org/10.1371/journal.pcbi.1009581.g004

Greengenes 13_8 hosts a large number of unique sequences (Fig 2B) and similar sequence
entropy to SILVA (Fig 2C), but yields many sequences that are unannotated at the genus
(54%) and species (90%) levels. This indicates that a large number of sequences in this database
are genetically similar (� 98%) but taxonomically distinct, yielding ambiguous labels. This
highlights practical disadvantages with using this database, as 99% OTU clustering (the highest
% similarity provided in the 13_8 release) limits sequence information that could be used to
differentiate groups (e.g., dereplication to 100% OTUs would preserve this sequence diversity).
In practice, the use of non-sequence information (e.g., ecological distribution) can be leveraged to guide taxonomic classification [56] and hence preserving this information can be
advantageous in some use cases. RESCRIPt provides users with a variety of taxonomy dereplication options to put such decisions in the hands of individual researchers, and to make data
processing pipelines transparent and reproducible so that others can reconstitute and adjust
processing decisions as desired.
To provide context to database size comparisons, we evaluated taxonomic overlap among
these databases, extending an earlier comparison of SILVA, Greengenes, NCBI (not RefSeq),
and other taxonomies by Balvočiūtė and Huson [73]. Fig 4 illustrates the proportion of labels
shared at each taxonomic rank, between each pair, trio, and across all databases, relative to the
total number of taxonomic labels in each database. Labels that were a prefix of another label
were collapsed into that label to avoid undercounting the number of shared labels between
databases, and to account for subclade labels used by GTDB (for example, “Lactobacillus_A”
would be considered the same label as “Lactobacillus”).

PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1009581 November 8, 2021

8 / 37

[с. 9]
PLOS COMPUTATIONAL BIOLOGY

Reproducible reference database management

We found that very low proportions of taxonomic labels were shared between and among
databases at all ranks below domain. In general, select pairs shared more labels than trios, and
~30–50% (of the total number of labels in each of those databases) of species labels were shared
by SILVA, GTDB, and NCBI-RefSeq (Fig 4). Proportions increased at the genus rank, and for
some groups at the species rank, reflecting differences in taxonomic labels (often due to taxonomic reclassifications related to database age) at the intermediate ranks. Taxonomic reclassifications have rendered many of the older taxonomies obsolete, most notably (and
unsurprisingly) the Greengenes 2013 release taxonomy, as reflected in the low proportions
shared with all other databases. Proposed taxonomic reclassifications in GTDB are unique to
that database, leading to reduced sharing with all other databases. SILVA exhibited relatively
low proportions of genus and species labels shared with other databases, but this is unsurprising given that SILVA does not curate species labels (as discussed above). Considering these
limitations, these findings indicate that a reasonably high proportion of genus and species
labels are shared across databases, and instances of non-sharing reflect taxonomic reclassifications more often than lack of coverage. Greengenes exhibited notably poorer taxonomic coverage, taking into account both the paucity of unique labels (Fig 3A), low taxonomic entropy
(Fig 3B), the markedly lower level of sharing with all other databases (compared to other pairs
of databases), and the high level of coverage of Greengenes taxonomies by the other databases
(Fig 4).
Using NCBI-RefSeq as the standard for coverage of official taxonomic names (as the only
taxonomy in this comparison that is comprised mostly of genomes sequenced from type material), SILVA exhibits the best coverage at class through species rank (Fig 4), but GTDB exhibits
only slightly lower coverage at those ranks (most likely because proposed reclassifications
reduce the degree of taxonomic sharing), hence both exhibit similar levels of coverage. By
comparison, only a small proportion of GTDB and SILVA species labels are shared with other
databases, reflecting both taxonomic inconsistencies (as described above) as well as the inclusion of uncultured and proposed taxa. Taken together, these findings reinforce the suggestion
that NCBI-RefSeq may contain the best coverage of accepted type strains, making it best suited
to some applications, although the greater inclusion of non-type and uncultured species in
SILVA and GTDB may make these databases more suitable for environmental survey applications and other studies containing many uncultured organisms.
To evaluate practical implications of reference database selection on 16S rRNA gene
sequence classification of biological data, we compared classification of the Earth Microbiome
Project (EMP) data [74] using SILVA, Greengenes, GTDB, and NCBI-RefSeq reference databases (S1 Text). Results demonstrate that classification with SILVA yielded the highest number
(S1A Fig) and entropy (S1B Fig) of taxonomic labels at genus and species ranks, and the lowest
proportion of unclassified sequences at order through species ranks (S1C Fig). Viewing taxonomic classification results at family level indicates that the different databases do not appreciably alter predicted abundances of different groups, but differences in taxonomic labeling
and taxonomic coverage between databases lead to some notable differences, and lower depth
of classification with Greengenes (Fig 5). Classification of real biological data is important for
contextualizing methodological benchmarking results with “real world” performance, but the
true result often cannot be known a priori [75]. Hence, we use the EMP data here to test performance for classification of biological sequences, but cannot use these results to compare
accuracy between the different reference databases. For example, the higher proportion of species-level classifications with SILVA is encouraging, but does not necessarily indicate that this
is a better result—indeed, the simulated classification results (Fig 3D) suggest that the specieslevel classifications achieved with SILVA have lower accuracy than the other databases.
Accordingly, classification of true biological data should be interpreted with caution.

PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1009581 November 8, 2021

9 / 37

[с. 10]
PLOS COMPUTATIONAL BIOLOGY

Reproducible reference database management

Fig 5. Average family-level taxonomic composition of EMP empo 3 types. Family-level classification as predicted by SILVA, Greengenes, GTDB, or NCBI-RefSeqs
classifiers. Samples were grouped by EMPO 3 type to look at average family-level taxonomic composition of each sample type. Only taxa detected at a minimum of 10%
relative frequency in at least one group are shown.
https://doi.org/10.1371/journal.pcbi.1009581.g005

Effects of processing steps on the SILVA 16S rRNA gene database
RESCRIPt contains multiple functions that were designed specifically for handling data from
SILVA, due to the popularity of SILVA as a reference for LSU and SSU rRNA gene sequences
(Fig 1). We tested the impact of several of these steps on the SILVA 138 release 16S rRNA
sequences to inform best practices for processing these data with RESCRIPt. Removing abnormally short sequences (fragments) (Fig 6A) and sequences with excessive ambiguity and
homopolymer content had the beneficial effect of reducing database size (Fig 6B) without substantial loss of sequence entropy (Fig 6C). Ambiguously labeled taxa (e.g., those unidentified at
genus or species ranks) may create “taxonomic noise” but clearly represent unique genotypes

PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1009581 November 8, 2021

10 / 37

[с. 11]
PLOS COMPUTATIONAL BIOLOGY

Reproducible reference database management

Fig 6. Comparison of sequence information across each successive sequence quality filtering step as applied to the SILVA 16S rRNA gene database. A, Sequence
length distributions. B, Number of unique sequences. C, Entropy of full-length sequences and different kmer lengths. Note: The subsequent sequence length filtering did
not have any effect on the data as the NR99 reference database is already pre-trimmed as specified above. Base: the complete NR99 SILVA database, Culled: after
sequences with either 8 or more homopolymers and/or 5 ambiguous bases removed, LengFiltByTax: sequence length filtering of the data based on taxonomy, i.e.
removal of archaeal and bacterial sequences less than 900 and 1200 bp in length, respectively. DereplicateUniq: Taxonomy and Sequence dereplication using “uniq”
mode (i.e. any identical sequences with differing taxonomy will not be merged), NoAmbigLabels: any sequence data associated with ambiguous labels (typically at lower
taxonomic ranks) are removed from the data set.
https://doi.org/10.1371/journal.pcbi.1009581.g006

and should not be removed (Fig 6B). Whereas the strict removal of any sequence containing
ambiguous taxonomic annotations (typically at the lower ranks) resulted in the removal of
305,636 sequences. This had a noticeable effect on sequence entropy.
The evaluate-taxonomy action was also used to examine the number of unique taxonomic
labels (Fig 7A), taxonomic entropy (Fig 7B). Optimal classification performance for each classifier without cross-validation (Fig 7C) as was performed earlier (Fig 7D), and with cross-validation (Fig 7D). Aside from the strict removal of ambiguous taxonomic annotations, quality
filtering also had minimal impact on classifier performance and entropy. Quality filtering has
a similarly subtle effect on taxonomic classification of biological data (S2 Fig).

Fig 7. Comparison of taxonomic information and simulated classification accuracy across several successive steps of quality filtering of the NR99 16S rRNA gene
databases. A, Number of unique taxonomic labels; B, Taxonomic entropy; C, optimal classification accuracy from the evaluate-fit-classifier action (as F-Measure)
without cross-validation (simulating best possible classification accuracy when the true label is known but classification accuracy may be confounded by other similar
hits in the database); D, optimal classification accuracy from the evaluate-cross-validate action (as F-Measure), which simulates pseudo-realistic classification task
whereby a set of query sequences may not have an exact match in the reference database. See Fig 6 Legend for label descriptions. Rank labels on x-axis: D = domain,
P = phylum, C = class, O = order, F = family, G = genus, S = species.
https://doi.org/10.1371/journal.pcbi.1009581.g007

PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1009581 November 8, 2021

11 / 37

