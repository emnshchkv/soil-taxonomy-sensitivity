# Часть 4: с. 12–15

[с. 12]
PLOS COMPUTATIONAL BIOLOGY

Reproducible reference database management

Effect of clustering on sequence and taxonomic information: Lessons from
the Greengenes 16S rRNA gene database
Clustering sequences into OTUs has long been practiced for dereplicating and reducing errors
in marker-gene sequencing experiments [76]. In times of yore, clustering was also often
applied to reference sequence databases for marker-gene sequencing, to reduce complexity
and thus computational requirements. We have previously shown that clustering COI diet
metabarcoding databases at 97% vs. 99% is detrimental to database information [77], but the
effects of clustering on database quality more generally are lacking. To benchmark general
effects of OTU clustering on database quality, and for 16S rRNA gene sequencing specifically,
we used RESCRIPt to evaluate multiple database quality characteristics of the Greengenes
database (13_8 release) [21,23] clustered at multiple OTU % similarity thresholds. The Greengenes public release data contain pre-clustered sequence and taxonomy files (with LCA consensus taxonomies assigned to OTU clusters), which were evaluated in this benchmark using
the RESCRIPt actions evaluate-taxonomy, evaluate-fit-classifier, evaluate-cross-validate, and
evaluate-classifications.
The loss of information as a result of clustering sequences has been highlighted by others
for bacterial SSU [44,78] and metazoan COI [77] query sequences. We build on their work by
demonstrating similar issues with the use of OTU clustering for reducing complexity in
marker-gene sequence reference databases. Decreasing % similarity threshold rapidly leads to
information loss; at the genus and species ranks the number of unique taxonomic labels rapidly declines (Fig 8A), as sequences (and genera and species) are collapsed into larger OTUs

Fig 8. Taxonomic information (A-C) and classification accuracy (D-E) of Greengenes 16S rRNA gene database clustered at different similarity thresholds. Subpanels
show taxonomic/classification characteristics at each taxonomic level: A, Number of unique taxonomic labels; B, Taxonomic entropy; C, number of taxa that terminate
at that level; D, optimal classification accuracy (as F-Measure) without cross-validation (simulating best possible classification accuracy when the true label is known but
classification accuracy may be confounded by other similar hits in the database); E, classification accuracy (F-Measure) with cross-validation (simulating realistic
classification tasks when the correct label is unknown). Rank labels on x-axis: D = domain, P = phylum, C = class, O = order, F = family, G = genus, S = species.
https://doi.org/10.1371/journal.pcbi.1009581.g008

PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1009581 November 8, 2021

12 / 37

[с. 13]
PLOS COMPUTATIONAL BIOLOGY

Reproducible reference database management

with less taxonomic resolution. Taxonomic entropy (Shannon’s entropy [79] applied to vectors
of taxonomic label counts), which measures both the richness and evenness of taxonomic
labels, registers a gradual decline as OTU % similarity is decreased from 99% to 88% at both
genus and species levels, and rapidly declines thereafter (Fig 8B). This indicates that, although
unique genus and species labels are being collapsed into larger family-level OTUs, OTU clustering is also initially reducing label redundancy, leading to increased evenness. The proportion of terminal labels (i.e., the rank at which taxonomic annotation terminates) illustrates
how the rank assignment landscape changes: a higher proportion of genus- and species-level
annotations in the 99% OTU sequences gives way to a higher proportion of class-, order-, and
family-level terminal annotations as the % similarity threshold is decreased (Fig 8C). Contrary
to this trend, “best-case” classification accuracy (with evaluate-fit-classifier, which trains and
tests a naive Bayes taxonomy classifier on the same input data without cross-validation) is seen
to increase from F = 0.88 to F = 1.0 as databases are clustered (Fig 8D), but this phenomenon
reflects the loss of information with increased OTU clustering, suggesting that the higher
apparent accuracy is actually an indicator of lost database coverage. Cross-validated classification accuracy (with evaluate-cross-validate) gives a more realistic demonstration of performance, demonstrating that classification accuracy actually declines as OTU clustering %
increases, as database coverage decreases making the classifier less effective (Fig 8E). Taken
together, these results suggest that even very modest OTU % clustering thresholds are likely to
negatively affect database information content as well as classification accuracy. Although a
small degree of dereplication and clustering may be beneficial for reducing sequence and taxonomic redundancy, we recommend against using OTU clustering < 99% similarity for any
marker-gene sequence databases.

Reference curation improves taxonomic classification: Lessons from the
UNITE Fungal ITS database
Next, we used RESCRIPt to benchmark the effects of various sequence processing steps on
fungal ITS sequences, and in particular the impacts of database clustering and representation
on taxonomic classification. Sequence reference databases are often subsetted by investigators
to focus on particular clades of interest or to perform additional curation of public datasets.
Some researchers have generated environment-specific databases, founded in the belief that
such databases increase taxonomic classification accuracy by removing sequences that are
genetically related but ecologically distinct from species found in a specific environment
[54,57–61], although this can elevate the risk of false-positive errors [80]. RESCRIPt contains
several methods to support and evaluate such filtering decisions, which then become embedded in provenance to facilitate transparent and reproducible use of these databases. To demonstrate this filtering capacity, and evaluate its effect on classification accuracy, we performed a
benchmark of the UNITE [31] ITS database. The most recent releases of UNITE contain different versions that we benchmark here: (1) different “species hypothesis” OTU clustering
thresholds (including a dynamically defined clustering threshold defined by manual curation
[31]); (2) release versions containing ITS sequences for all eukaryotes vs. only fungi; and (3)
the fungal database filtered (by RESCRIPt) to contain only sequences that are annotated at the
order level or below.
OTU clustering (97% vs. 99% vs. dynamic clustering) exhibited minimal impact on results,
though the 99% OTUs yielded the highest taxonomic information and classification accuracy
(Fig 9), corresponding to the more comprehensive clustering benchmark performed above
(Fig 8). As expected, the “all eukaryotes” version of UNITE contains more than twice as many
sequences as the fungi-only database (Fig 9A), though taxonomic entropy is slightly lower

PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1009581 November 8, 2021

13 / 37

[с. 14]
PLOS COMPUTATIONAL BIOLOGY

Reproducible reference database management

Fig 9. Taxonomic information (A–C) and classification accuracy (D–E) of UNITE ITS domain database with different filtering and clustering settings. Filtered versions
include the "all Eukaryotes" 2020.04.02 release version containing all Eukaryotes ("All Euks"), filtered to contain only Fungi, and filtered to contain only Fungi with at
least order-level taxonomic annotation ("Fungi Order"). Cluster levels indicate which UNITE release version was used: sequences clustered at 97% similarity, 99%
similarity, or the UNITE "dynamic" species hypothesis threshold. Subpanels show taxonomic/classification characteristics at each taxonomic level: A, Number of unique
taxonomic labels; B, Taxonomic entropy; C, proportion of taxa that terminate at that level; D, optimal classification accuracy (as F-Measure) without cross-validation
(simulating best possible classification accuracy when the true label is known but classification accuracy may be confounded by other similar hits in the database); E,
classification accuracy (F-Measure) with cross-validation (simulating realistic classification tasks when the correct label is unknown). Rank labels on x-axis:
K = kingdom, P = phylum, C = class, O = order, F = family, G = genus, S = species. See Materials and Methods for more details on how these databases were created and
processed.
https://doi.org/10.1371/journal.pcbi.1009581.g009

(Fig 9B), reflecting a greater proportion of non-fungal sequences that are not annotated at the
family, genus, and species ranks (Fig 9C). Taxonomic classification accuracy is lower in the
“all eukaryote” release version of the UNITE database, compared to fungal sequences alone
(Fig 9D and 9E). This indicates that removing non-target sequences from the database
improves taxonomic resolution; however, such practices (whether focusing on particular
clades or environment-specific species) is fraught with risks and should be used with caution.
If the non-target sequences can be detected (e.g., amplified by the same primer, or introduced
by cross-contamination or rare migration events), filtered databases may lead to misclassification (e.g. a sequence may be classified as a fungus, when it is actually a metazoan, or vice versa
[80]). We recommend careful consideration of these risks when selecting primers, databases,
and filtering decisions for marker-gene and metagenome sequencing studies. An advantage of
using RESCRIPt for database construction and curation is that these processing steps are
embedded in provenance, allowing the appropriateness of these steps to be re-evaluated at
later stages (e.g., in documenting results, peer review, and re-use in future studies and by other
researchers).
Filtering out fungal sequences that were not annotated at least to order level removed only
a small fraction of unique labels (Fig 9A), boosting entropy by a narrow margin (Fig 9B) and
classification accuracy by a wide margin (Fig 9D and 9E). These results indicate that, even in

PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1009581 November 8, 2021

14 / 37

[с. 15]
PLOS COMPUTATIONAL BIOLOGY

Reproducible reference database management

curated release versions of some public databases, some additional curation is beneficial to
remove sequences with missing or uninformative taxonomic labels. Filtering with RESCRIPt
enables researchers to automatically record these filtering decisions in provenance, making it
clear when, where, and why their reference materials diverge from public release versions of
these databases.

Clustering and primer-region trimming effects on a BOLD COI gene
database
The COI gene is a common target for taxonomic identification of metazoa, both for diet metabarcoding and eDNA studies [35]. The earliest versioned COI databases were available
through a single resource: the Barcode of Life Data Systems (BOLD) [34]. However, a growing
number of researchers have recently contributed updated COI databases using either BOLD
or NCBI GenBank (or both) as source data [14,81–84]. We conducted a series of benchmarks
that evaluated the compositional and performance effects of COI databases constructed by
varying the following: first, clustering reference sequences and/or using full length sequences
versus trimming references to a particular primer region; second, the reference source itself
(BOLD vs. NCBI GenBank) [34,85]. The BOLD vs. NCBI GenBank comparison is the subject
of the next section. For all COI benchmarks, the results are reported separately for the two
largest groups of (animal) COI sequence data available in BOLD: arthropods and chordates.
Previous reports have separately described the reduction of taxonomic information when
clustering COI sequences, as well as the effects of classifier performance based on the reference
source [77,83,84]. We build on those previous works to demonstrate the effects of clustering in
combination with trimming these sequences to a particular region within the COI sequence,
focusing only on reference sequences obtained from BOLD. Clustering sequences dramatically
reduces the number of unique sequences in both the untrimmed and primer-trimmed COI
datasets (Fig 10A). Similarly, trimming COI sequences to a particular primer-defined region
results in a decrease in the number of unique chordate and arthropod sequences. Nevertheless,
it is worth noting that despite trimming to a region containing sequences approximately 180
bp in length, a large amount of sequence diversity remains for both arthropod (N = 611,166)
and chordate (N = 69,924) references. Sequence entropy was similarly reduced when references were trimmed or clustered (Fig 10B).
We found that sequence clustering and primer trimming both also contribute to a decrease
in taxonomic information, most pronounced at the species rank (Fig 11A). Untrimmed chordate sequences clustered at 97% identity (“Full_97”) contained just 77% as many unique labels
as the untrimmed sequences that were only dereplicated (“Full_100”). Trimming these chordate sequences to a particular primer region also reduced taxonomic information, with dereplicated (“ANML_100”) and clustered (“ANML_97”) sequences containing 80% and 60% as
many unique labels as the untrimmed dereplicated Chordate reference set, respectively. A similar effect was observed for arthropod sequences, with 97% identity clustered untrimmed
sequences containing 83% as many unique labels compared to the dereplicated and
untrimmed arthropod references. Dereplicated and primer-trimmed arthropod sequences,
and 97% clustered primer-trimmed sequences contained 78% and 62% as many labels relative
to the untrimmed, dereplicated arthropod sequences, respectively.
Unlike the observations with Greengenes SSU data, we found that taxonomic entropy
increased marginally when a larger clustering radius (decreasing percent identity) is applied
(Fig 11B), indicating that clustering reduces redundancy of taxonomic labels. In addition, the
same effect is observed when a reference sequence is reduced to a relatively shorter subsequence. Likewise, we observed that both clustering and primer trimming led to a significant

PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1009581 November 8, 2021

15 / 37

