# CRC_MasterThesis

<i> “Evaluation of early-stage detection of colorectal cancer using machine learning models based on functional profiling of human gut microbiome” </i> 

### Background
CRC prediction models based on functional profiling of gut microbiome

### Results

### Conclusion

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

#### 1.2 


#### 1.3 Functional profiling
Functional profiles of high-quality metagenomic shotgun sequences were determined using [HUMAnN 3.0](https://github.com/biobakery/humann) [(Francesco _et al.,_ 2020)](https://elifesciences.org/articles/65088).

#### 1.4 Taxonomic profiling

Taxonomic profiles

```bash
humann.sh
```

### 2. Metadata and feature table preparation

### 3. Explorative analysis (ampvis2)

[ampvis2](https://madsalbertsen.github.io/ampvis2/index.html)

### 4. Association testing

### 5. Machine learning models

#

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
