# ruff: noqa: E401,E501,E701,E702,E731,I001
import ast, json, re
from pathlib import Path
from types import SimpleNamespace
import pytest
from scrna_integration.run_contract import atomic_write_json, determine_stage_status, prepare_run, promote_run, resume_run, sha256_file, validate_artifacts, validate_checkpoint

NB = Path(__file__).parents[1] / "notebooks/06b_per_cluster.ipynb"
def src(code=False):
    return ["".join(c.get("source", [])) for c in json.loads(NB.read_text())["cells"] if not code or c["cell_type"] == "code"]
def fn(name, ns):
    tree = ast.parse("\n\n".join(src(True))); node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == name)
    exec(compile(ast.Module([node], []), str(NB), "exec"), ns); return ns[name]
def upstream(tmp, stage, status="SUCCESS", subset=False):
    p = prepare_run(tmp / "up", "input"); cp = p.draft_dir / "checkpoint.h5ad"; cp.write_bytes(b"cp")
    m = {"run_id": p.run_id, "stage": stage, "stage_status": status, "checkpoint": {"path": cp.name, "sha256": sha256_file(cp)}}
    if status == "SUCCESS_WITH_WARNINGS": m["warning_acceptance"] = {"accepted_by": "PI", "accepted_at": "2026-07-14"}
    if subset:
        art = p.draft_dir / "subset.h5ad"; art.write_bytes(b"sub"); m["artifacts"] = [{"role": "subset", "path": art.name, "sha256": sha256_file(art)}]
    atomic_write_json(p.manifest_path, m); promote_run(p) if status != "FAILED" else p.draft_dir.replace(p.promoted_dir)
    return p, p.promoted_dir / ("subset.h5ad" if subset else cp.name)
def loader(mode, p, reader):
    ns = dict(MODE=mode, UPSTREAM_RUN_ROOT=p.run_dir.parent, UPSTREAM_RUN_ID=p.run_id, json=json, sc=SimpleNamespace(read_h5ad=reader), resume_run=resume_run, validate_checkpoint=validate_checkpoint, validate_artifacts=validate_artifacts)
    return fn("_load_verified_upstream", ns)

@pytest.mark.parametrize("mode,stage,status,subset", [("global","06_annotated","SUCCESS",False),("subset","06c_subset","SUCCESS_WITH_WARNINGS",True)])
def test_valid_modes(tmp_path, mode, stage, status, subset):
    p, selected = upstream(tmp_path, stage, status, subset); marker = object()
    assert loader(mode, p, lambda path: marker if path == selected else None)()[:2] == (marker, selected)
@pytest.mark.parametrize("bad", ["stage","status","missing","tamper","read"])
def test_bad_upstream_never_claims_run(tmp_path, bad):
    subset = bad in {"missing","tamper"}; p, selected = upstream(tmp_path, "wrong" if bad == "stage" else ("06c_subset" if subset else "06_annotated"), "FAILED" if bad == "status" else "SUCCESS", subset and bad != "missing")
    if bad == "tamper": selected.write_bytes(b"bad")
    with pytest.raises((ValueError,OSError)): loader("subset" if subset else "global", p, lambda _: (_ for _ in ()).throw(OSError("read")) if bad == "read" else object())()
    assert not (tmp_path / "runs" / "06b").exists()

def test_reserved_names_and_real_write_failures(tmp_path):
    safe = fn("_safe_report_dirname", {"_SAFE_DIRNAME_RE": re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*\Z"), "_RESERVED_REPORT_DIRNAMES": {"manifest.json","draft","promoted"}})
    assert safe("reports") and all(not safe(v) for v in ("MANIFEST.JSON","Draft","PROMOTED","a/b",".."))
    errors=[]; ns={"_report_errors":errors}; ns["_report_error"]=fn("_report_error",ns); write=fn("_write_report",ns)
    assert not write(tmp_path,"x") and not write(tmp_path,"x","a")
    fig=SimpleNamespace(savefig=lambda *a,**k: (_ for _ in ()).throw(OSError("savefig"))); ns.update(plt=SimpleNamespace(close=lambda _:None)); save=fn("_save_report_figure",ns); links=[]
    assert not save(fig,tmp_path/"x.png",links,"![x](x.png)") and not links and {e["operation"] for e in errors} == {"text_x","text_a","savefig"}

def finalizer(tmp, ref="ok.png", warnings=()):
    p=prepare_run(tmp/"runs","06b"); report=p.draft_dir/"reports"; report.mkdir(); index=report/"index.md"; cluster=report/"cluster-a.md"; image=report/"ok.png"; image.write_bytes(b"png")
    index.write_text("index"); cluster.write_text(f"![x]({ref})" if ref else "no optional figure")
    ns=dict(Path=Path,re=re,shutil=__import__("shutil"),atomic_write_json=atomic_write_json,sha256_file=sha256_file,validate_checkpoint=validate_checkpoint,validate_artifacts=validate_artifacts,determine_stage_status=determine_stage_status,run_paths=p,REPORT_DIR=report,RUN_ID=p.run_id,MODE="global",_upstream_manifest_path=tmp/"m.json",_upstream_manifest_sha256="0"*64,effective_parameters={},runtime_provenance={},_warnings=list(warnings),_method_status={"per_cluster_reports":"success"},_report_errors=[])
    ns["_safe_cluster_slug"]=fn("_safe_cluster_slug",{"hashlib":__import__("hashlib"),"re":re}); ns["_report_error"]=fn("_report_error",ns)
    return fn("_finalize_run",ns),p,report,index,{"a":{"out_md":str(cluster)}}

@pytest.mark.parametrize("ref,warnings,status", [("ok.png",(),"SUCCESS"),(None,("optional figure unavailable",),"SUCCESS_WITH_WARNINGS"),("missing.png",(),"FAILED"),("../escape.png",(),"FAILED")])
def test_final_status_and_image_completeness(tmp_path, ref, warnings, status):
    finalize,p,report,index,summaries=finalizer(tmp_path,ref,warnings)
    if status == "FAILED":
        with pytest.raises(RuntimeError): finalize(index,summaries)
        m=json.loads(p.manifest_path.read_text()); assert m["method_status"]["per_cluster_reports"] == "failed" and "artifacts" not in m and not report.exists()
    else:
        assert finalize(index,summaries).value == status and validate_checkpoint(p.manifest_path) == index and p.draft_dir.is_dir() and not p.promoted_dir.exists()

def test_science_and_preflight_contract_retained():
    s="\n".join(src()); required=("rank_genes_groups","邻居簇对比 DEG","sc.pl.umap","sc.pl.dotplot","基因集评分","跨疾病丰度","call_llm_for_annotation","Priority Review",'"auto_promote": False')
    assert all(x in s for x in required) and s.index("_load_verified_upstream()") < s.index("run_paths = prepare_run(RUN_ROOT, RUN_ID)") and "promote_run(" not in s
    branch=s[s.index("if _preflight_errors:"):s.index("# 合法输入在任何报告写入前")]; assert '"stage_status": "FAILED"' in branch and "REPORT_DIR.mkdir" not in branch and "checkpoint" not in branch and "artifacts" not in branch
