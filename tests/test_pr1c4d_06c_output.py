# ruff: noqa: E401,E501,E701,E702,I001
import ast, json, re
from pathlib import Path
from types import SimpleNamespace
import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scrna_integration.run_contract import atomic_write_json, prepare_run, promote_run, resume_run, sha256_file, snapshot_effective_parameters, validate_artifacts, validate_checkpoint

NB=Path(__file__).parents[1]/"notebooks/06c_subset.ipynb"
def code():
    return "\n\n".join("".join(c.get("source",[])) for c in json.loads(NB.read_text())["cells"] if c["cell_type"]=="code")
def fn(name,ns):
    tree=ast.parse(code()); node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name==name); exec(compile(ast.Module([node],[]),str(NB),"exec"),ns); return ns[name]
def objects():
    obs=pd.DataFrame({"cell_type_final_v1":["A","A","B","B"]},index=["c1","c2","c3","c4"])
    main=ad.AnnData(sparse.csr_matrix(np.eye(4,dtype=np.float32)),obs=obs.copy()); subset=main[["c1","c2"]].copy()
    main.uns.update({"stage":"06_annotated","status":"SUCCESS","run_id":"stage06"})
    subset.obs["cell_type_final_subset_v1"]=["A1","A2"]; main.obs["cell_type_final_subset_v1"]=["A1","A2",None,None]; main.obs["cell_type_unified_v2"]=["A1","A2","B","B"]
    return subset,main,obs["cell_type_final_v1"].copy()
def namespace(tmp,versions=("v1","v1","v2"),run_id="06c-output"):
    subset,main,coarse=objects(); manifest=tmp/"upstream.json"; checkpoint=tmp/"upstream.bin"; manifest.write_text("{}"); checkpoint.write_bytes(b"upstream")
    ns=dict(Path=Path,re=re,RUN_ROOT=tmp/"runs",RUN_ID=run_id,UPSTREAM_RUN_ID="stage06",UPSTREAM_LABEL_VERSION=versions[0],SUBSET_OUTPUT_VERSION=versions[1],MAIN_OUTPUT_VERSION=versions[2],UPSTREAM_LABEL_COL=f"cell_type_final_{versions[0]}",SUBSET_LABEL_COL=f"cell_type_final_subset_{versions[1]}",MAIN_LABEL_COL=f"cell_type_unified_{versions[2]}",SUBSET_LABELS=["A"],_verified_upstream_index=main.obs_names.copy(),_verified_upstream_coarse=coarse,_verified_upstream_obs_columns=("cell_type_final_v1",),_root=tmp,_upstream_manifest_path=manifest,_upstream_manifest_sha256=sha256_file(manifest),_upstream_checkpoint=checkpoint,_upstream_checkpoint_sha256=sha256_file(checkpoint),_OUTPUT_VERSION_RE=re.compile(r"v[1-9][0-9]*\Z"),_H5AD_BASENAME_RE=re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*\.h5ad\Z"),prepare_run=prepare_run,atomic_write_json=atomic_write_json,sha256_file=sha256_file,snapshot_effective_parameters=snapshot_effective_parameters,collect_runtime_provenance=lambda *_:{},validate_checkpoint=validate_checkpoint,validate_artifacts=validate_artifacts)
    for name in ("_series_equal","_output_contract","_manifest_base","_record_output_failure","_complete_output_write","_write_06c_outputs"): ns[name]=fn(name,ns)
    return ns,subset,main

@pytest.mark.parametrize("versions",[("v1","bad","v2"),("v1","v1","v1"),("v1","v01","v2"),("v1","v1","V2")])
def test_unsafe_or_conflicting_versions_write_failed_without_h5ad(tmp_path,versions):
    ns,subset,main=namespace(tmp_path,versions); writes=[]; subset.write_h5ad=lambda *a,**k:writes.append(a); main.write_h5ad=lambda *a,**k:writes.append(a)
    with pytest.raises(ValueError): ns["_write_06c_outputs"](subset,main)
    p=SimpleNamespace(manifest_path=tmp_path/"runs"/"06c-output"/"draft"/"manifest.json"); m=json.loads(p.manifest_path.read_text())
    assert not writes and m["stage_status"]=="FAILED" and "checkpoint" not in m and "artifacts" not in m

@pytest.mark.parametrize("failure",["snapshot","provenance"])
def test_manifest_base_failure_does_not_claim_run(tmp_path,failure):
    ns,subset,main=namespace(tmp_path); ns["snapshot_effective_parameters" if failure=="snapshot" else "collect_runtime_provenance"]=lambda *a,**k: (_ for _ in ()).throw(RuntimeError(failure))
    with pytest.raises(RuntimeError,match=failure): ns["_write_06c_outputs"](subset,main)
    assert not (tmp_path/"runs"/"06c-output").exists()

def test_wrong_unified_values_fail_without_retained_outputs(tmp_path):
    ns,subset,main=namespace(tmp_path); main.obs.loc["c1","cell_type_unified_v2"]="WRONG"
    with pytest.raises(ValueError): ns["_write_06c_outputs"](subset,main)
    draft=tmp_path/"runs"/"06c-output"/"draft"; manifest=json.loads((draft/"manifest.json").read_text())
    assert not list(draft.glob("*.h5ad")) and manifest["stage_status"]=="FAILED" and "checkpoint" not in manifest and "artifacts" not in manifest

def test_literal_na_label_never_equals_missing_or_passes_gate(tmp_path):
    ns,subset,main=namespace(tmp_path); left=pd.Series(["<NA>"],index=["c1"]); right=pd.Series([pd.NA],index=["c1"],dtype="string")
    assert not ns["_series_equal"](left,right)
    subset.obs.loc["c1","cell_type_final_subset_v1"]="<NA>"; main.obs.loc["c1","cell_type_final_subset_v1"]="<NA>"; main.obs.loc["c1","cell_type_unified_v2"]=pd.NA
    with pytest.raises(ValueError): ns["_write_06c_outputs"](subset,main)
    draft=tmp_path/"runs"/"06c-output"/"draft"; manifest=json.loads((draft/"manifest.json").read_text()); assert not list(draft.glob("*.h5ad")) and manifest["stage_status"]=="FAILED" and "checkpoint" not in manifest and "artifacts" not in manifest

def test_second_write_failure_removes_both_outputs_and_keeps_original_error(tmp_path):
    ns,subset,main=namespace(tmp_path); subset.write_h5ad=lambda path,**_:Path(path).write_bytes(b"subset")
    def fail(path,**_): Path(path).write_bytes(b"partial"); raise OSError("second write failed")
    main.write_h5ad=fail
    with pytest.raises(OSError,match="second write failed"): ns["_write_06c_outputs"](subset,main)
    draft=tmp_path/"runs"/"06c-output"/"draft"; m=json.loads((draft/"manifest.json").read_text())
    assert not list(draft.glob("*.h5ad")) and m["stage_status"]=="FAILED" and "checkpoint" not in m and "artifacts" not in m

def test_real_two_anndata_roundtrip_manifest_and_review_gate(tmp_path):
    ns,subset,main=namespace(tmp_path); upstream_top=dict(main.uns); paths,subset_path,main_path=ns["_write_06c_outputs"](subset,main)
    artifacts=validate_artifacts(paths.manifest_path); manifest=json.loads(paths.manifest_path.read_text()); subset_rt=ad.read_h5ad(subset_path); main_rt=ad.read_h5ad(main_path)
    assert validate_checkpoint(paths.manifest_path)==subset_path and set(artifacts)=={"subset","main_reflow"} and manifest["stage_status"]=="NEEDS_REVIEW"
    assert subset_rt.obs_names.equals(pd.Index(["c1","c2"])) and main_rt.obs["cell_type_final_v1"].tolist()==["A","A","B","B"]
    assert main_rt.obs["cell_type_final_subset_v1"].iloc[:2].tolist()==["A1","A2"] and main_rt.obs["cell_type_final_subset_v1"].iloc[2:].isna().all()
    assert main_rt.obs["cell_type_unified_v2"].tolist()==["A1","A2","B","B"] and set(main_rt.obs.columns)=={"cell_type_final_v1","cell_type_final_subset_v1","cell_type_unified_v2"}
    assert all(main_rt.uns[k]==v for k,v in upstream_top.items()) and subset_rt.uns["stage"]=="06c_subset" and subset_rt.uns["status"]=="NEEDS_REVIEW"
    with pytest.raises(ValueError,match="cannot be promoted"): promote_run(paths)

@pytest.mark.parametrize("role",["subset","main_reflow"])
def test_tamper_either_artifact_blocks_validation_and_resume(tmp_path,role):
    ns,subset,main=namespace(tmp_path); paths,_,_=ns["_write_06c_outputs"](subset,main); target=validate_artifacts(paths.manifest_path)[role]; target.write_bytes(b"tampered")
    with pytest.raises(ValueError): validate_checkpoint(paths.manifest_path)
    with pytest.raises(ValueError): resume_run(tmp_path/"runs",paths.run_id)

def test_collision_refuses_overwrite_and_no_legacy_flat_outputs(tmp_path):
    ns,subset,main=namespace(tmp_path); paths,subset_path,main_path=ns["_write_06c_outputs"](subset,main); before=(sha256_file(subset_path),sha256_file(main_path))
    ns2,subset2,main2=namespace(tmp_path)
    with pytest.raises(FileExistsError): ns2["_write_06c_outputs"](subset2,main2)
    assert before==(sha256_file(subset_path),sha256_file(main_path)) and not list(tmp_path.glob("*.h5ad"))
    source=code(); assert "OUTPUT_PATH" not in source and "MAIN_OUTPUT_PATH" not in source and "promote_run(" not in source
