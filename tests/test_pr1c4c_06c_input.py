# ruff: noqa: E401,E501,E701,E702,E731,I001
import ast, json, re
from pathlib import Path
from types import SimpleNamespace
import pandas as pd
import pytest
from scrna_integration.run_contract import atomic_write_json, prepare_run, promote_run, resume_run, sha256_file, snapshot_effective_parameters, validate_artifacts, validate_checkpoint

NB = Path(__file__).parents[1] / "notebooks/06c_subset.ipynb"
def cells(code=False):
    return ["".join(c.get("source", [])) for c in json.loads(NB.read_text())["cells"] if not code or c["cell_type"] == "code"]
def fn(name, ns):
    tree=ast.parse("\n\n".join(cells(True))); node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name==name)
    exec(compile(ast.Module([node],[]),str(NB),"exec"),ns); return ns[name]
def upstream(tmp, stage="06_annotated", status="SUCCESS", warning=True, artifact=True):
    p=prepare_run(tmp/"up","stage06"); cp=p.draft_dir/"checkpoint.h5ad"; cp.write_bytes(b"checkpoint")
    m={"run_id":p.run_id,"stage":stage,"stage_status":status,"checkpoint":{"path":cp.name,"sha256":sha256_file(cp)}}
    if status=="SUCCESS_WITH_WARNINGS" and warning: m["warning_acceptance"]={"accepted_by":"PI","accepted_at":"2026-07-14"}
    if artifact:
        art=p.draft_dir/"report.txt"; art.write_bytes(b"report"); m["artifacts"]=[{"role":"report","path":art.name,"sha256":sha256_file(art)}]
    atomic_write_json(p.manifest_path,m)
    if status=="SUCCESS" or status=="SUCCESS_WITH_WARNINGS" and warning: promote_run(p)
    else: p.draft_dir.replace(p.promoted_dir)
    return p,cp.name
def loader(p, reader):
    ns=dict(UPSTREAM_RUN_ROOT=p.run_dir.parent,UPSTREAM_RUN_ID=p.run_id,json=json,sc=SimpleNamespace(read_h5ad=reader),resume_run=resume_run,validate_checkpoint=validate_checkpoint,validate_artifacts=validate_artifacts,sha256_file=sha256_file)
    return fn("_load_verified_upstream",ns)

@pytest.mark.parametrize("status",["SUCCESS","SUCCESS_WITH_WARNINGS"])
def test_valid_promoted_stage06_records_selected_checkpoint(tmp_path,status):
    p,name=upstream(tmp_path,status=status); marker=object(); selected=p.promoted_dir/name
    assert loader(p,lambda path: marker if path==selected else None)()[:2]==(marker,selected)
@pytest.mark.parametrize("bad",["stage","status","warning","checkpoint","artifact","read"])
def test_bad_upstream_or_read_never_claims_current_run(tmp_path,bad):
    p,name=upstream(tmp_path,stage="wrong" if bad=="stage" else "06_annotated",status="FAILED" if bad=="status" else ("SUCCESS_WITH_WARNINGS" if bad=="warning" else "SUCCESS"),warning=bad!="warning")
    if bad=="checkpoint": (p.promoted_dir/name).write_bytes(b"bad")
    if bad=="artifact": (p.promoted_dir/"report.txt").write_bytes(b"bad")
    with pytest.raises((ValueError,OSError)): loader(p,lambda _: (_ for _ in ()).throw(OSError("read")) if bad=="read" else object())()
    assert not (tmp_path/"runs"/"06c").exists()

def preflight(frame,labels=("A",),minimum=1,versions=("v1","v1","v2")):
    current=SimpleNamespace(obs=frame,obs_names=frame.index,n_obs=len(frame)); upstream_v,subset_v,main_v=versions
    ns=dict(pd=pd,_VERSION_RE=re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*\Z"),UPSTREAM_LABEL_VERSION=upstream_v,SUBSET_OUTPUT_VERSION=subset_v,MAIN_OUTPUT_VERSION=main_v,UPSTREAM_LABEL_COL=f"cell_type_final_{upstream_v}",SUBSET_LABELS=labels,MIN_SUBSET_CELLS=minimum)
    return fn("_subset_preflight",ns)(current)
@pytest.mark.parametrize("kind",["missing_col","na","empty","empty_labels","labels_str","labels_none","duplicate","unknown","small","min_bool","min_str","min_negative","bad_version","all"])
def test_invalid_local_preflight_cases(kind):
    frame=pd.DataFrame({"cell_type_final_v1":["A","A","B"]},index=["c1","c2","c3"]); labels=("A",); minimum=1
    if kind=="missing_col": frame.columns=["other"]
    if kind=="na": frame.iloc[0,0]=None
    if kind=="empty": frame.iloc[0,0]=""
    if kind=="empty_labels": labels=()
    if kind=="labels_str": labels="A"
    if kind=="labels_none": labels=None
    if kind=="duplicate": labels=("A","A")
    if kind=="unknown": labels=("X",)
    if kind=="small": minimum=3
    if kind=="min_bool": minimum=True
    if kind=="min_str": minimum="1"
    if kind=="min_negative": minimum=-1
    versions=(None,"v1","v2") if kind=="bad_version" else ("v1","v1","v2")
    if kind=="all": labels=("A","B")
    errors,mask=preflight(frame,labels,minimum,versions); assert errors and mask.dtype==bool and mask.index.equals(frame.index) and (kind not in {"labels_str","labels_none","min_bool","min_str","min_negative","bad_version"} or not mask.any())

def test_valid_mask_versions_do_not_claim_run(tmp_path):
    errors,mask=preflight(pd.DataFrame({"cell_type_final_v3":["A","B","B"]},index=["x","y","z"]),("A",),1,("v3","s2","m4"))
    assert not errors and mask.tolist()==[True,False,False] and mask.index.tolist()==["x","y","z"] and not (tmp_path/"runs").exists()
@pytest.mark.parametrize("versions,labels",[((None,"v1","v2"),("A",)),(("v1","v1","v2"),"A"),(("v1","v1","v2"),None)])
def test_raw_invalid_parameters_write_only_failed_audit(tmp_path,versions,labels):
    errors,mask=preflight(pd.DataFrame({"cell_type_final_v1":["A","B"]}),labels,1,versions); assert errors and not mask.any()
    manifest=tmp_path/"upstream.json"; checkpoint=tmp_path/"upstream.h5ad"; manifest.write_text("{}"); checkpoint.write_bytes(b"x")
    ns=dict(RUN_ROOT=tmp_path/"runs",RUN_ID="06c",_root=tmp_path,_upstream_manifest_path=manifest,_upstream_manifest_sha256=sha256_file(manifest),_upstream_checkpoint=checkpoint,_upstream_checkpoint_sha256=sha256_file(checkpoint),prepare_run=prepare_run,atomic_write_json=atomic_write_json,snapshot_effective_parameters=snapshot_effective_parameters,collect_runtime_provenance=lambda *_:{})
    p=fn("_write_preflight_failure",ns)(["bad labels"]); audit=json.loads(p.manifest_path.read_text())
    assert audit["stage_status"]=="FAILED" and len(audit["inputs"])==2 and "checkpoint" not in audit and "artifacts" not in audit and not (p.draft_dir/"reports").exists()

def test_params_no_eval_or_fixed_path_and_science_preserved():
    code="\n".join(cells(True)); joined="\n".join(cells()); params=cells(True)[0]
    for name in ("UPSTREAM_RUN_ROOT","UPSTREAM_RUN_ID","RUN_ROOT","RUN_ID","UPSTREAM_LABEL_VERSION","SUBSET_OUTPUT_VERSION","MAIN_OUTPUT_VERSION","SUBSET_LABELS","MIN_SUBSET_CELLS"): assert name in params
    assert "UPSTREAM_PATH" not in code and "SUBSET_FILTER" not in code and "eval(" not in code and "\nOUTPUT_VERSION =" not in code
    assert all(stale not in joined for stale in ("cell_type_final_v1","cell_type_final_subset_v1","cell_type_unified_v1","OUTPUT_DIR","normalize_v1","embedding_v1","clustering_v1"))
    assert all(name in joined for name in ("UPSTREAM_LABEL_COL","SUBSET_LABEL_COL","MAIN_LABEL_COL","MODE","UPSTREAM_RUN_ID","RUN_ROOT","REPORT_DIRNAME"))
    for token in ("highly_variable_genes","run_harmony","SCVI","leiden","mLLMCelltype","pi_subset_decisions","标签回流"): assert token in code
