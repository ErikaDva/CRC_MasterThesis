# CRC_MasterThesis

### Title

<i> “Evaluation of early-stage detection of colorectal cancer using machine learning models based on functional profiling of human gut microbiome” </i> 

### Background
CRC prediction models based on functional profiling of gut microbiome

### Results

### Conclusion

## Repository structure

This repository contains the code and brief description of the workflow used for the CRC microbiome analysis in Master's thesis that was undertaken at the [Albertsen lab](https://albertsenlab.org/) (AAU).

#### Files

#### Data

#### Figures

#### Scripts

## Workflow

### 1. Data collection & processing

#### 1.1 Data availability

The publicly available raw sequencing data from CRC studies were used in this analysis and are available from the European Nucleotide Archive (ENA). The metadata for each study was obtained according to the information provided by the researchers in their papers.

| Study | Country | ENA |
| --- | --- | --- |
| Feng _et al.,_ 2015 | Austria | ERP008729 |
| Gupta _et al.,_ 2019 | India | PRJNA53127 |
| Yamada _et al.,_ 2019 | Japan | DRA006684 & DRA008156 |
| Thomas _et al.,_ 2019 | Italy | SRP136711 |
| Vogtmann _et al.,_ 2016 | USA | PRJEB12449 |
| Wirbel _et al.,_ 2019 | Germany | PRJEB27928 |
| Yu _et al.,_ 2015 | China | PRJEB12449 |
| Zeller _et al.,_ 2014 | France | ERP005534 |

#### 1.2 Raw data pre-processing


#### 1.3 Functional profiling
Functional profiles of high-quality metagenomic shotgun sequences were determined using [HUMAnN 3.0](https://github.com/biobakery/humann) [(Francesco _et al.,_ 2020)](https://elifesciences.org/articles/65088).

```bash
bash humann.sh
```

#### 1.4 Taxonomic profiling

Taxonomic profiles


### 2. Metadata and feature table preparation

#### 2.1 Metadata

#### 2.2 Feature tables


Multi-sequenced samples were merged together taking into the account the library sizes using these scripts:
[`2.1_prepare_functional_data.R`](https://github.com/ErikaDva/CRC_MasterThesis/blob/main/scripts/test)\
[`2.2_prepare_pathway_data.R`](https://github.com/ErikaDva/CRC_MasterThesis/blob/main/scripts/test)\
[`2.3_prepare_taxonomic_data.R`](https://github.com/ErikaDva/CRC_MasterThesis/blob/main/scripts/test)

### 3. Explorative analysis (ampvis2)

R package, [ampvis2](https://madsalbertsen.github.io/ampvis2/index.html), was utilised for explorative analysis of functional and taxonomic feature tables. The package was originally developed for visualing amplicon data, however, it is capable of dealing with shotgun metagenomics data.

### 4. Association testing

#### 4.1 SIAMCAT/Wirbel

#### 4.2 MaAsLin2

R package, [MaAslin2](https://huttenhower.sph.harvard.edu/maaslin/), was used to determine multivariable associations between metadata and microbial features [(Mallick _et al.,_ 2021)](https://doi.org/10.1186/s13059-021-02306-1). It is based on general linear models.

### 5. Machine learning models

Machine learning models were built using [SIAMCAT](https://siamcat.embl.de/) pipeline for associations between gut microbiome and host phenotype [(Wirbel _et al.,_ 2021)](https://doi.org/10.1186/s13059-021-02306-1).

#### 5.1 LASSO

#### 5.2 Random Forest

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

## Contact
