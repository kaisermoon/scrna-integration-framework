"""IO module: multi-source scRNA-seq data readers + CellxGene-compatible schema.

Planned readers:
- cellranger (filtered_feature_bc_matrix / raw_feature_bc_matrix)
- h5ad (AnnData native)
- RData (via anndata2ri or rds2py)
- 10x .h5
- mtx + tsv triplet

Planned schema utility:
- normalize_obs_to_cellxgene(adata, project_meta) -> AnnData
"""
