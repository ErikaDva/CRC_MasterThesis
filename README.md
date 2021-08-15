# CRC_MasterThesis

<i> “Evaluation of early-stage detection of colorectal cancer using machine learning models based on functional profiling of human gut microbiome” </i> 

### Background
CRC prediction models based on functional profiling of gut microbiome

### Results

### Conclusion

## Workflow

### 1. Data collection & processing

#### 1.1 Data availability

| Study | ENA |
| --- | --- |
| Yu _et al.,_ 2015 | PRJEB12449 |
| Feng _et al.,_ 2015 | ERP008729 |
|  |  |
|  |  |
| Gupta _et al.,_ 2019 | PRJNA53127 |
| Wirbel _et al.,_ 2019 | PRJEB27928 |
| Zeller _et al.,_ 2014 | ERP005534 |

Functional profiles of high-quality metagenomic shotgun sequences were determined using [HUMAnN 3.0](https://github.com/biobakery/humann) [(Francesco _et al.,_ 2020)](https://elifesciences.org/articles/65088).

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
