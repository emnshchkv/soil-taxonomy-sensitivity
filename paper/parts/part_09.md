# Часть 9: с. 30–32

[с. 30]
PLOS COMPUTATIONAL BIOLOGY

Reproducible reference database management

taxonomic rank. Rank labels on x-axis: D = Domain, P = phylum, C = class, O = order,
F = family, G = genus, S = species.
(EPS)
S2 Fig. Comparison of SILVA database quality filtering on taxonomic classification of
Earth Microbiome Project sequences. A, Number of unique taxonomic labels; B, Taxonomic
entropy; and C, proportion of unclassified taxa at each taxonomic rank. Rank labels on x-axis:
D = Domain, P = phylum, C = class, O = order, F = family, G = genus, S = species. Raw V4: the
complete NR99 SILVA database, with V4 sequences extracted, followed by dereplication,
Default V4: same as Raw, but filtered to remove sequences with either 8 or more homopolymers and/or 5 ambiguous bases removed, and to remove Archaeal and Bacterial sequences less
than 900 and 1200 bp in length, respectively. Strict V4: same as Default, but removing
sequences with any ambiguous nucleotides, NoAmbigLabels V4: same as Default, but removing any sequence associated with ambiguous labels (typically empty annotations at genus or
species ranks).
(EPS)
S1 File. Example procedure for constructing a 12S rRNA marker gene reference database.
A jupyter notebook version is available at https://github.com/bokulich-lab/db-benchmarks2020.
(PDF)

Acknowledgments
We would like to thank the QIIME 2 developers and users for providing feedback on
RESCRIPt.

Author Contributions
Conceptualization: Michael S. Robeson, II, Benjamin D. Kaehler, Nicholas A. Bokulich.
Data curation: Michael S. Robeson, II, Devon R. O’Rourke, Matthew R. Dillon, Nicholas A.
Bokulich.
Formal analysis: Michael S. Robeson, II, Devon R. O’Rourke, Benjamin D. Kaehler, Nicholas
A. Bokulich.
Investigation: Michael S. Robeson, II, Devon R. O’Rourke, Nicholas A. Bokulich.
Methodology: Michael S. Robeson, II, Devon R. O’Rourke, Benjamin D. Kaehler, Michal
Ziemski, Matthew R. Dillon, Nicholas A. Bokulich.
Project administration: Nicholas A. Bokulich.
Resources: Matthew R. Dillon.
Software: Michael S. Robeson, II, Benjamin D. Kaehler, Michal Ziemski, Matthew R. Dillon,
Nicholas A. Bokulich.
Supervision: Jeffrey T. Foster, Nicholas A. Bokulich.
Validation: Michael S. Robeson, II, Benjamin D. Kaehler, Michal Ziemski, Matthew R. Dillon,
Nicholas A. Bokulich.
Visualization: Michael S. Robeson, II, Devon R. O’Rourke, Benjamin D. Kaehler, Nicholas A.
Bokulich.

PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1009581 November 8, 2021

30 / 37

[с. 31]
PLOS COMPUTATIONAL BIOLOGY

Reproducible reference database management

Writing – original draft: Michael S. Robeson, II, Devon R. O’Rourke, Benjamin D. Kaehler,
Michal Ziemski, Nicholas A. Bokulich.
Writing – review & editing: Michael S. Robeson, II, Devon R. O’Rourke, Benjamin D. Kaehler, Michal Ziemski, Matthew R. Dillon, Jeffrey T. Foster, Nicholas A. Bokulich.

References
1.

Caporaso JG, Lauber CL, Walters WA, Berg-Lyons D, Lozupone C, Turnbaugh PJ, et al. Global patterns of 16S rRNA diversity at a depth of millions of sequences per sample. Proc Natl Acad Sci U S A.
2011; 108 Suppl 1: 4516–4522.

2.

Tedersoo L, Bahram M, Põlme S, Kõljalg U, Yorou NS, Wijesundera R, et al. Fungal biogeography.
Global diversity and geography of soil fungi. Science. 2014; 346: 1256688–1256688. https://doi.org/
10.1126/science.1256688 PMID: 25430773

3.

Consortium THMP, The Human Microbiome Project Consortium. Structure, function and diversity of
the healthy human microbiome. Nature. 2012. pp. 207–214. https://doi.org/10.1038/nature11234
PMID: 22699609

4.

Bokulich NA, Chung J, Battaglia T, Henderson N, Jay M, Li H, et al. Antibiotics, birth mode, and diet
shape microbiome maturation during early life. Sci Transl Med. 2016; 8: 343ra82. https://doi.org/10.
1126/scitranslmed.aad7121 PMID: 27306664

5.

Nayfach S, Shi ZJ, Seshadri R, Pollard KS, Kyrpides NC. New insights from uncultivated genomes of
the global human gut microbiome. Nature. 2019; 568: 505–510. https://doi.org/10.1038/s41586-0191058-x PMID: 30867587

6.

Vorholt JA, Vogel C, Carlström CI, Müller DB. Establishing Causality: Opportunities of Synthetic Communities for Plant Microbiome Research. Cell Host & Microbe. 2017. pp. 142–155. https://doi.org/10.
1016/j.chom.2017.07.004 PMID: 28799900

7.

Bokulich NA, Thorngate JH, Richardson PM, Mills DA. Microbial biogeography of wine grapes is conditioned by cultivar, vintage, and climate. Proc Natl Acad Sci U S A. 2014; 111: E139–48. https://doi.org/
10.1073/pnas.1317377110 PMID: 24277822

8.

Wagg C, Schlaeppi K, Banerjee S, Kuramae EE, van der Heijden MGA. Fungal-bacterial diversity and
microbiome complexity predict ecosystem functioning. Nat Commun. 2019; 10: 4841. https://doi.org/
10.1038/s41467-019-12798-y PMID: 31649246

9.

Robeson MS 2nd, Khanipov K, Golovko G, Wisely SM, Bodenchuck M, et al. Assessing the utility of
metabarcoding for diet analyses of the omnivorous wild pig (Sus scrofa). Ecol Evol. 2018; 8: 185–196.
https://doi.org/10.1002/ece3.3638 PMID: 29321862

10.

Bergmann GT, Craine JM, Robeson MS 2nd, Fierer N. Seasonal Shifts in Diet and Gut Microbiota of
the American Bison (Bison bison). Maldonado JE, editor. PLoS One. 2015; 10: e0142409. https://doi.
org/10.1371/journal.pone.0142409 PMID: 26562019

11.

Kartzinel TR, Chen PA, Coverdale TC, Erickson DL, Kress WJ, Kuzmina ML, et al. DNA metabarcoding illuminates dietary niche partitioning by African large herbivores. Proc Natl Acad Sci U S A. 2015;
112: 8019–8024. https://doi.org/10.1073/pnas.1503283112 PMID: 26034267

12.

Deiner K, Bik HM, Mächler E, Seymour M, Lacoursière-Roussel A, Altermatt F, et al. Environmental
DNA metabarcoding: Transforming how we survey animal and plant communities. Mol Ecol. 2017; 26:
5872–5895. https://doi.org/10.1111/mec.14350 PMID: 28921802

13.

Creer S, Deiner K, Frey S, Porazinska D, Taberlet P, Thomas WK, et al. The ecologist’s field guide to
sequence-based identification of biodiversity. Methods Ecol Evol. 2016; 7: 1008–1018.

14.

Porter TM, Hajibabaei M. Automated high throughput animal CO1 metabarcode classification. Sci
Rep. 2018; 8: 4226. https://doi.org/10.1038/s41598-018-22505-4 PMID: 29523803

15.

Taberlet P, Coissac E, Hajibabaei M, Rieseberg LH. Environmental DNA. Molecular Ecology.
2012. pp. 1789–1793. https://doi.org/10.1111/j.1365-294X.2012.05542.x PMID: 22486819

16.

Méric G, Wick RR, Watts SC, Holt KE, Inouye M. Correcting index databases improves metagenomic
studies. bioRxiv. 2019; 2: e000075.

17.

Almeida A, Mitchell AL, Tarkowska A, Finn RD. Benchmarking taxonomic assignments based on 16S
rRNA gene profiling of the microbiota from commonly sampled environments. Gigascience. 2018; 7.
https://doi.org/10.1093/gigascience/giy054 PMID: 29762668

18.

Woese CR, Fox GE. Phylogenetic structure of the prokaryotic domain: the primary kingdoms. Proc
Natl Acad Sci U S A. 1977; 74: 5088–5090. https://doi.org/10.1073/pnas.74.11.5088 PMID: 270744

PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1009581 November 8, 2021

31 / 37

[с. 32]
PLOS COMPUTATIONAL BIOLOGY

Reproducible reference database management

19.

Woese CR, Fox GE, Pechman KR. Comparative Cataloging of 16S Ribosomal Ribonucleic Acid:
Molecular Approach to Procaryotic Systematics. Int J Syst Evol Microbiol. 1977; 27: 44–57.

20.

Quast C, Pruesse E, Yilmaz P, Gerken J, Schweer T, Yarza P, et al. The SILVA ribosomal RNA gene
database project: improved data processing and web-based tools. Nucleic Acids Res. 2013; 41:
D590–6. https://doi.org/10.1093/nar/gks1219 PMID: 23193283

21.

DeSantis TZ, Hugenholtz P, Larsen N, Rojas M, Brodie EL, Keller K, et al. Greengenes, a chimerachecked 16S rRNA gene database and workbench compatible with ARB. Appl Environ Microbiol.
2006; 72: 5069–5072. https://doi.org/10.1128/AEM.03006-05 PMID: 16820507

22.

Pruesse E, Quast C, Knittel K, Fuchs BM, Ludwig W, Peplies J, et al. SILVA: a comprehensive online
resource for quality checked and aligned ribosomal RNA sequence data compatible with ARB. Nucleic
Acids Res. 2007; 35: 7188–7196. https://doi.org/10.1093/nar/gkm864 PMID: 17947321

23.

McDonald D, Price MN, Goodrich J, Nawrocki EP, DeSantis TZ, Probst A, et al. An improved Greengenes taxonomy with explicit ranks for ecological and evolutionary analyses of bacteria and archaea.
ISME J. 2012; 6: 610–618. https://doi.org/10.1038/ismej.2011.139 PMID: 22134646

24.

Parks DH, Chuvochina M, Waite DW, Rinke C, Skarshewski A, Chaumeil P-A, et al. A standardized
bacterial taxonomy based on genome phylogeny substantially revises the tree of life. Nat Biotechnol.
2018; 36: 996–1004. https://doi.org/10.1038/nbt.4229 PMID: 30148503

25.

Parks DH, Chuvochina M, Chaumeil P-A, Rinke C, Mussig AJ, Hugenholtz P. Selection of representative genomes for 24,706 bacterial and archaeal species clusters provide a complete genome-based
taxonomy. Microbiology. bioRxiv; 2019. p. 820.

26.

O’Leary NA, Wright MW, Brister JR, Ciufo S, Haddad D, McVeigh R, et al. Reference sequence
(RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation. Nucleic
Acids Res. 2016; 44: D733–45. https://doi.org/10.1093/nar/gkv1189 PMID: 26553804

27.

Roux S, Enault F, Bronner G, Debroas D. Comparison of 16S rRNA and protein-coding genes as
molecular markers for assessing microbial diversity (Bacteria and Archaea) in ecosystems. FEMS
Microbiol Ecol. 2011; 78: 617–628. https://doi.org/10.1111/j.1574-6941.2011.01190.x PMID:
22066608

28.

Stefan CP, Koehler JW, Minogue TD. Targeted next-generation sequencing for the detection of ciprofloxacin resistance markers using molecular inversion probes. Sci Rep. 2016; 6: 25904. https://doi.
org/10.1038/srep25904 PMID: 27174456

29.

Dahllöf I, Baillie H, Kjelleberg S. rpoB-based microbial community analysis avoids limitations inherent
in 16S rRNA gene intraspecies heterogeneity. Appl Environ Microbiol. 2000; 66: 3376–3380. https://
doi.org/10.1128/AEM.66.8.3376-3380.2000 PMID: 10919794

30.

Schoch CL, Seifert KA, Huhndorf S, Robert V, Spouge JL, Levesque CA, et al. Nuclear ribosomal
internal transcribed spacer (ITS) region as a universal DNA barcode marker for Fungi. Proc Natl Acad
Sci U S A. 2012; 109: 6241–6246. https://doi.org/10.1073/pnas.1117018109 PMID: 22454494

31.

Kõljalg U, Larsson K-H, Abarenkov K, Nilsson RH, Alexander IJ, Eberhardt U, et al. UNITE: a database providing web-based methods for the molecular identification of ectomycorrhizal fungi. New Phytol. 2005; 166: 1063–1068. https://doi.org/10.1111/j.1469-8137.2005.01376.x PMID: 15869663

32.

Deshpande V, Wang Q, Greenfield P, Charleston M, Porras-Alfaro A, Kuske CR, et al. Fungal identification using a Bayesian classifier and the Warcup training set of internal transcribed spacer
sequences. Mycologia. 2016; 108: 1–5. https://doi.org/10.3852/14-293 PMID: 26553774

33.

Cole JR, Wang Q, Cardenas E, Fish J, Chai B, Farris RJ, et al. The Ribosomal Database Project:
improved alignments and new tools for rRNA analysis. Nucleic Acids Res. 2009; 37: D141–5. https://
doi.org/10.1093/nar/gkn879 PMID: 19004872

34.

Ratnasingham S, Hebert PDN. BOLD: The Barcode of Life Data System (www.barcodinglife.org). Mol
Ecol Notes. 2007; 7: 355–364. https://doi.org/10.1111/j.1471-8286.2007.01678.x PMID: 18784790

35.

Hebert PDN, Cywinska A, Ball SL, deWaard JR. Biological identifications through DNA barcodes.
Proc Biol Sci. 2003; 270: 313–321. https://doi.org/10.1098/rspb.2002.2218 PMID: 12614582

36.

Folmer O, Black M, Hoeh W, Lutz R, Vrijenhoek R. DNA primers for amplification of mitochondrial
cytochrome c oxidase subunit I from diverse metazoan invertebrates. Mol Mar Biol Biotechnol. 1994;
3: 294–299. PMID: 7881515

37.

Waraniak JM, Marsh TL, Scribner KT. 18S rRNA metabarcoding diet analysis of a predatory fish community across seasonal changes in prey availability. Ecol Evol. 2019; 9: 1410–1430. https://doi.org/10.
1002/ece3.4857 PMID: 30805170

38.

James D, Schmidt A-M. Use of an intron region of a chloroplast tRNA gene (trnL) as a target for PCR
identification of specific food crops including sources of potential allergens. Food Res Int. 2004; 37:
395–402.

PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1009581 November 8, 2021

32 / 37

