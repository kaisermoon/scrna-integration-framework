# YAML manifest for dataset facts, despite "YAML is hard to configure"

Per-source-dataset configuration lives in a `data/{source_dataset}/manifest.yaml` file read by `read_with_manifest`, **not** as a Python dict in the stage-1 notebook PARAMS cell. This appears to contradict the original 项目构思 ("YAML 比较难懂且难以配置；配置尽量放 Jupyter"), so the reasoning is recorded here.

The apparent contradiction dissolves once two kinds of configuration are separated:

- **Write-once dataset facts** — how a given GEO dataset's obs columns are named, how its clinical table joins, which author-annotation columns exist, whether the author already removed doublets, the species/ontology constants. Set once when a dataset is first ingested, then essentially frozen. These go in the **YAML manifest**.
- **Tunable analysis parameters** — QC thresholds, HVG count, clustering resolution, which embedding to promote. Adjusted repeatedly across re-runs. These stay in the **notebook `# === PARAMS ===` cell**, exactly as the 项目构思 wanted.

PI confirmed this split directly: "数据集自己本身的配置，配好了基本上不需要改；但是参数那些要改的。" The 项目构思's aversion to YAML was about the *tunable* knobs — putting those in YAML would be genuinely painful. Dataset facts are not knobs.

## Considered Options

- **Manifest as a Python dict in the stage-1 notebook** (option B in the grilling): rejected. A notebook processes one dataset at a time, so 5+ source datasets would need 5 notebook copies or a dict-list that re-introduces the complexity it tried to avoid. The declaration also stops travelling with the data — reproducing on another machine or handing a dataset to a collaborator would require carrying the notebook too. Cross-disease multi-source integration (the project's whole point) needs the per-dataset declaration to live with the dataset and be version-controlled.
- **Manifest as YAML, full 95-line schema always written**: rejected as the default presentation. The full schema is a reference, not a per-dataset requirement; presenting it as the norm would make every dataset look like it needs 95 lines and would re-create the "YAML is hard" problem the 项目构思 warned about.

## Consequences

- Manifest is YAML, kept minimal: six required fields (~8 lines) cover a clean dataset; the seven optional blocks (`obs_mapping` / `value_mapping` / `clinical_metadata` / `ontology` / `project_specific` / `original_annotations` beyond `[]` / `qc_overrides`) appear **only when that dataset actually has the heterogeneity they describe**.
- SPEC.md presents a minimal example first, then the full-featured example labelled as a reference superset.
- Students do touch YAML, but only for stable dataset facts they set once — accepted by PI as reasonable. The reproducibility win (declaration version-controlled alongside data) outweighs the learning cost.
- The write-once vs tunable boundary is the standing rule for deciding whether any new piece of configuration belongs in the manifest or in PARAMS.
