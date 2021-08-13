# CRC_MasterThesis
CRC prediction models based on functional profiling of gut microbiome

# Workflow

## 1.

## 2. Data preparation

## 3. Explorative analysis (ampvis2)

## 4. Association testing

## 5. Machine learning models

## Set-up & Data preparation

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
