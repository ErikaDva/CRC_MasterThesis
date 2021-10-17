# CRC_MasterThesis_21

### Title

<i> “Investigation of early-stage detection of colorectal cancer using machine learning models based on functional profiling of the human gut microbiome” </i> 

### Background
CRC prediction models based on functional profiling of gut microbiome.

### Results

_TBU_

### Conclusion

_TBU_

## Repository structure

This repository contains the code and brief description of the workflow used for the CRC microbiome analysis in Master's thesis that was undertaken at the [Albertsen lab](https://albertsenlab.org/) (AAU).

#### Files

This folder contains generated files in this analysis:
* e.g. `predictions_lasso.tsv`

#### Data

This folder contains data & other relevant information:
* Feature tables*
* Metadata
* Library sizes
* Run accessions numbers

#### Figures

This folder contains the generated figures.

#### Models

This folder contains models based on each feature category. Models are saved as R objects and can be loaded to R environment without re-running the script.

#### Scripts

This folder contains the source code used in this project.

>\* Files too large to be uploaded in this repository

## Workflow

### 1. Data collection & processing

#### 1.1 Data availability

The publicly available raw sequencing data from CRC studies were used in this analysis and are available on the European Nucleotide Archive (ENA) at EMBL-EBI.

| Study | Country | Accession number(s) |
| --- | --- | --- |
| Feng _et al.,_ 2015 | Austria | [ERP008729](https://www.ebi.ac.uk/ena/browser/view/ERP008729) |
| Gupta _et al.,_ 2019 | India | [PRJNA531273](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA531273) & [PRJNA39711](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA397112) |
| Yachida _et al.,_ 2019 | Japan | [DRA006684](https://www.ebi.ac.uk/ena/browser/view/DRA006684) & [DRA008156](https://www.ebi.ac.uk/ena/browser/view/DRA008156) |
| Thomas _et al.,_ 2019 | Italy | [SRP136711](https://www.ebi.ac.uk/ena/browser/view/SRP136711) |
| Vogtmann _et al.,_ 2016 | USA | [PRJEB12449](https://www.ebi.ac.uk/ena/browser/view/PRJEB12449) |
| Wirbel _et al.,_ 2019 | Germany | [PRJEB27928](https://www.ebi.ac.uk/ena/browser/view/PRJEB27928) |
| Yu _et al.,_ 2015 | China | [PRJEB12449](https://www.ebi.ac.uk/ena/browser/view/PRJEB12449) |
| Zeller _et al.,_ 2014 | France | [ERP005534](https://www.ebi.ac.uk/ena/browser/view/ERP005534) |

>Note: The list of all studies including metadata was retrieved and modified from [Wirbel _et al.,_ 2021](https://doi.org/10.1186/s13059-021-02306-1). The metadata from Gupta _et al.,_ 2019 and Yamada _et al.,_ 2019 were seperately obtained according to the information provided by the researchers in their papers.

#### 1.2 Raw data pre-processing

Command line was used to pre-process raw data and to perform taxonomic and functional profiling of high-quality reads.

[TrimGalore (v.0.6.5) wrapper package](https://github.com/FelixKrueger/TrimGalore) with filtering parameters:\
`--stringency 5 -- length 45 --quality 20 --max_n 2 --trim-n --paired`  

Bowtie2 (v.2.3.4.1) was used to align reads to the human genome (hg19) and discard contaminant reads from the host.

#### 1.3 Functional profiling
Functional profiles of high-quality metagenomic shotgun sequences were determined using [HUMAnN 3.0](https://github.com/biobakery/humann) [(Francesco _et al.,_ 2020)](https://elifesciences.org/articles/65088).

```bash
bash humann.sh
```

#### 1.4 Taxonomic profiling

Taxonomic profiles with MetaPhlAn 3.0

```bash
bash mpa.sh
```

### 2. Metadata and feature table preparation

The metadata and feature table preparation was carried out in RStudio using R.

#### 2.1 Metadata

The metadata was retrieved and modified from [Wirbel _et al.,_ 2019](https://doi.org/10.1038/s41591-019-0406-6). The metadata from Gupta _et al.,_ 2019 and Yamada _et al.,_ 2019 were seperately obtained according to the information provided by the researchers in their papers.

The final metadata used in this project can be found in [CRC_MasterThesis/data/meta](https://github.com/ErikaDva/CRC_MasterThesis/blob/main/scripts/test) folder as `"meta.crc.tsv"`.

**or**  

can be generated by running the script:  
[`2.1_prepare_metadata.R`](https://github.com/ErikaDva/CRC_MasterThesis/blob/main/scripts/test)

#### 2.2 Feature tables

The feature tables produced by [HUMAnN 3.0](https://github.com/biobakery/humann) were subjected to post-processing in R Studio.

Firstly, multi-sequenced samples were merged together taking into the account the library sizes using these scripts:  
[`2.2_prepare_functional_data.R`](https://github.com/ErikaDva/CRC_MasterThesis/blob/main/scripts/test)\
[`2.2_prepare_pathway_data.R`](https://github.com/ErikaDva/CRC_MasterThesis/blob/main/scripts/test)\
[`2.2_prepare_gene_data.R`](https://github.com/ErikaDva/CRC_MasterThesis/blob/main/scripts/test)\
[`2.2_prepare_taxonomic_data.R`](https://github.com/ErikaDva/CRC_MasterThesis/blob/main/scripts/test)

>Note: exception with CN-CRC study as the number of samples matches the metadata entries

Secondly, the feature tables were cleaned and filtered to remove low-abundant features:  
[`2.2_clean_functional_data.R`](https://github.com/ErikaDva/CRC_MasterThesis/blob/main/scripts/test)\
[`2.2_clean_pathway_data.R`](https://github.com/ErikaDva/CRC_MasterThesis/blob/main/scripts/test)\
[`2.2_clean_taxonomic_data.R`](https://github.com/ErikaDva/CRC_MasterThesis/blob/main/scripts/test)

### 3. Explorative analysis

#### Overview of different profilers

#### Ordination with ampvis2

R package, [ampvis2](https://madsalbertsen.github.io/ampvis2/index.html), was utilised for explorative analysis of functional and taxonomic feature tables. The package was originally developed for visualing amplicon data, however, it is capable of dealing with shotgun metagenomics data.  
[`3_explorative_analysis_ampvis2.R`](https://github.com/ErikaDva/CRC_MasterThesis/blob/main/scripts/3_explorative_analysis_ampvis2.R)

### 4. Machine learning models

Machine learning models were built using [SIAMCAT](https://siamcat.embl.de/) pipeline for associations between gut microbiome and host phenotype [(Wirbel _et al.,_ 2021)](https://doi.org/10.1186/s13059-021-02306-1).

#### LASSO

The machine learning scripts were run in the following order:

4.1_train_models_species.R
4.1_train_models_functions.R

4.2_model_predictions.R

4.3_ml_figures.R

## Contact
E-mail: [Erika Dvarionaite](mailto:erika.dvarionaite@outlook.com)


Extra

### Set-up & Data preparation

1. import_feat_suffix()
2. remove_name_and_export()
3. 

## Summary of Custom functions

| Function | Description |
| --- | --- |
| `load_kegg()` | Import feature table with KEGG annotations (`"kegg"` can be substituted with `"pfam"`, `"eggnog"`, `"level4ec"` & `"go"` as required) |
| `load_meta()` | Import metadata file |
| `import_feat_suffix(filename, suffix)` | Import feature table & attach suffix (e.g. `import_feat_suffix("kegg_relab.tsv", kegg)` |
| `load_meta()` | Import metadata file |

### 4. Association testing

#### 4.1 SIAMCAT/Wirbel

#### 4.2 MaAsLin2

R package, [MaAslin2](https://huttenhower.sph.harvard.edu/maaslin/), was used to determine multivariable associations between metadata and microbial features [(Mallick _et al.,_ 2021)](https://doi.org/10.1186/s13059-021-02306-1). It is based on general linear models.
