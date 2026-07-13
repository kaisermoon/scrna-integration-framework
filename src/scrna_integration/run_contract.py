"""各阶段共用的运行状态、manifest 与 checkpoint 技术管道。"""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import platform
import re
import subprocess
import sys
import tempfile
import types
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any

_REDACTED = "<redacted>"
_SENSITIVE_PARAMETER_TOKENS = {"TOKEN", "SECRET", "PASSWORD", "CREDENTIAL", "CREDENTIALS"}
_RESEARCH_KEY_ALLOWLIST = {
    "BATCH_KEY",
    "CELL_TYPE_KEY",
    "CLUSTER_KEY",
    "DISEASE_KEY",
    "DONOR_KEY",
    "GROUP_KEY",
    "HVG_KEY",
    "LABEL_KEY",
    "SAMPLE_KEY",
}


def _is_sensitive_parameter_name(name: str) -> bool:
    if name in _RESEARCH_KEY_ALLOWLIST:
        return False
    tokens = name.upper().split("_")
    if any(token in _SENSITIVE_PARAMETER_TOKENS for token in tokens):
        return True
    sensitive_pairs = {("ACCESS", "KEY"), ("API", "KEY"), ("AUTH", "HEADER"), ("PRIVATE", "KEY")}
    return name.endswith("_KEY") or any(
        (left, right) in sensitive_pairs for left, right in zip(tokens, tokens[1:], strict=False)
    )


def _safe_parameter_path(value: str | Path, *, path_root: Path | None) -> str:
    path = Path(value)
    if not path.is_absolute():
        return path.as_posix() if isinstance(value, Path) else value
    absolute = path.resolve()
    if path_root is not None and absolute.is_relative_to(path_root):
        return absolute.relative_to(path_root).as_posix()
    return f"<external>/{absolute.name}"


def _json_parameter_value(
    value: Any, *, parameter_name: str, path_root: Path | None
) -> Any:
    if isinstance(value, Enum):
        return _json_parameter_value(value.value, parameter_name=parameter_name, path_root=path_root)
    if isinstance(value, str):
        return _safe_parameter_path(value, path_root=path_root)
    if value is None or isinstance(value, (int, bool)):
        return value
    if isinstance(value, float):
        try:
            json.dumps(value, allow_nan=False)
        except ValueError as error:
            raise TypeError(f"parameter {parameter_name} is not JSON serializable") from error
        return value
    if isinstance(value, Path):
        return _safe_parameter_path(value, path_root=path_root)
    if isinstance(value, Mapping):
        if any(not isinstance(key, str) for key in value):
            raise TypeError(f"parameter {parameter_name} has a non-string mapping key")
        return {
            key: _json_parameter_value(
                value[key], parameter_name=parameter_name, path_root=path_root
            )
            for key in sorted(value)
        }
    if isinstance(value, (tuple, list)):
        return [
            _json_parameter_value(item, parameter_name=parameter_name, path_root=path_root)
            for item in value
        ]
    item_method = getattr(value, "item", None)
    if callable(item_method):
        try:
            scalar = item_method()
        except Exception as error:
            raise TypeError(
                f"parameter {parameter_name} has an unusable scalar value"
            ) from error
        if scalar is not value:
            return _json_parameter_value(scalar, parameter_name=parameter_name, path_root=path_root)
    raise TypeError(
        f"parameter {parameter_name} has unsupported value type {type(value).__name__}"
    )


def snapshot_effective_parameters(
    namespace: Mapping[str, Any], exclude: Sequence[str] = (), *, path_root: str | Path | None = None
) -> dict[str, Any]:
    """捕获 notebook 中公开的大写参数，并生成稳定、脱敏的 JSON 数据。"""

    excluded = set(exclude)
    root = Path(path_root).resolve() if path_root is not None else None
    snapshot: dict[str, Any] = {}
    for name in sorted(namespace):
        if not isinstance(name, str) or name in excluded or name.startswith("_") or not name.isupper():
            continue
        if _is_sensitive_parameter_name(name):
            snapshot[name] = _REDACTED
            continue
        value = namespace[name]
        if isinstance(value, types.ModuleType) or isinstance(value, type) or callable(value):
            continue
        snapshot[name] = _json_parameter_value(value, parameter_name=name, path_root=root)
    return snapshot


def collect_runtime_provenance(
    project_root: str | Path, package_names: Sequence[str] = ()
) -> dict[str, Any]:
    """采集可复现运行所需的解释器、平台、Git 与包版本事实。"""

    root = Path(project_root).resolve()
    try:
        git_commit = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
            timeout=10,
        ).stdout.strip()
        git_diff = subprocess.run(
            ["git", "-C", str(root), "diff", "--binary", "HEAD"],
            check=True,
            capture_output=True,
            timeout=10,
        ).stdout
        git_available = True
        git_dirty: bool | None = bool(git_diff)
        git_diff_sha256 = hashlib.sha256(git_diff).hexdigest() if git_diff else None
    except (FileNotFoundError, subprocess.SubprocessError):
        git_available = False
        git_commit = "unavailable"
        git_dirty = None
        git_diff_sha256 = None

    packages: dict[str, str] = {}
    for name in sorted(set(package_names)):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = "unavailable"
    return {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "git_available": git_available,
        "git_commit": git_commit,
        "git_dirty": git_dirty,
        "git_diff_sha256": git_diff_sha256,
        "packages": packages,
    }


def publish_compatibility_symlink(link_path: str | Path, target_path: str | Path) -> Path:
    """原子创建相对 compatibility symlink，既有目标只允许幂等复用。"""

    link = Path(os.path.abspath(os.fspath(link_path)))
    target = Path(target_path).resolve()
    if not target.exists():
        raise FileNotFoundError(target)
    if not link.parent.is_dir():
        raise FileNotFoundError(link.parent)
    if link.is_symlink():
        if link.resolve() == target:
            return link
        raise FileExistsError(link)
    if link.exists():
        raise FileExistsError(link)
    link.symlink_to(os.path.relpath(target, start=link.parent))
    return link


class MethodStatus(str, Enum):
    """单个方法的四种执行状态。"""

    SUCCESS = "success"
    SKIPPED_BY_USER = "skipped_by_user"
    UNAVAILABLE = "unavailable"
    FAILED = "failed"


class StageStatus(str, Enum):
    """阶段计算和人工决策的四种状态。"""

    FAILED = "FAILED"
    NEEDS_REVIEW = "NEEDS_REVIEW"
    SUCCESS_WITH_WARNINGS = "SUCCESS_WITH_WARNINGS"
    SUCCESS = "SUCCESS"


def determine_stage_status(
    required_methods: Mapping[str, MethodStatus | str],
    hard_postconditions: Mapping[str, bool],
    *,
    needs_review: bool = False,
    warnings: Sequence[str] = (),
    allow_no_required_methods: bool = False,
) -> StageStatus:
    """根据已发生的执行事实确定 stage 状态，不代替科研选择。"""

    if not hard_postconditions:
        raise ValueError("hard_postconditions must not be empty")
    if not required_methods and not allow_no_required_methods:
        raise ValueError("required_methods must not be empty unless explicitly allowed")
    method_statuses = [MethodStatus(value) for value in required_methods.values()]
    if any(status is not MethodStatus.SUCCESS for status in method_statuses):
        return StageStatus.FAILED
    if any(not passed for passed in hard_postconditions.values()):
        return StageStatus.FAILED
    if needs_review:
        return StageStatus.NEEDS_REVIEW
    if warnings:
        return StageStatus.SUCCESS_WITH_WARNINGS
    return StageStatus.SUCCESS


@dataclass(frozen=True)
class RunPaths:
    """一个不可覆盖 run 的 draft/promoted 路径。"""

    run_id: str
    run_dir: Path

    @property
    def draft_dir(self) -> Path:
        return self.run_dir / "draft"

    @property
    def promoted_dir(self) -> Path:
        return self.run_dir / "promoted"

    @property
    def manifest_path(self) -> Path:
        return self.draft_dir / "manifest.json"


_RUN_ID_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*\Z")


def prepare_run(root: str | Path, run_id: str) -> RunPaths:
    """为新 RUN_ID 建立 draft；已存在的 RUN_ID 一律拒绝覆盖。"""

    if not _RUN_ID_RE.fullmatch(run_id):
        raise ValueError("run_id may contain only letters, digits, '.', '_' and '-'")
    root_path = Path(root).resolve()
    root_path.mkdir(parents=True, exist_ok=True)
    paths = RunPaths(run_id=run_id, run_dir=root_path / run_id)
    paths.run_dir.mkdir()
    paths.draft_dir.mkdir()
    return paths


def sha256_file(path: str | Path, *, chunk_size: int = 1024 * 1024) -> str:
    """分块计算输入或 checkpoint 的 SHA256，避免将大文件载入内存。"""

    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        while chunk := stream.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_write_json(path: str | Path, payload: Mapping[str, Any], *, overwrite: bool = False) -> None:
    """在目标目录原子写 JSON；默认不覆盖任何既有证据。"""

    destination = Path(path)
    if not destination.parent.is_dir():
        raise FileNotFoundError(destination.parent)
    fd, temporary_name = tempfile.mkstemp(
        dir=destination.parent, prefix=f".{destination.name}.", suffix=".tmp"
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as stream:
            json.dump(payload, stream, ensure_ascii=False, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        if overwrite:
            os.replace(temporary, destination)
        else:
            os.link(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def _load_manifest(manifest_path: Path) -> dict[str, Any]:
    with manifest_path.open(encoding="utf-8") as stream:
        manifest = json.load(stream)
    if not isinstance(manifest, dict):
        raise ValueError("manifest must be a JSON object")
    return manifest


def validate_checkpoint(manifest_path: str | Path, checkpoint_path: str | Path | None = None) -> Path:
    """校验 manifest 声明的 checkpoint 路径和内容 hash，可用于安全恢复。"""

    manifest_file = Path(manifest_path).resolve()
    manifest = _load_manifest(manifest_file)
    checkpoint = manifest.get("checkpoint")
    if not isinstance(checkpoint, dict):
        raise ValueError("manifest checkpoint must be an object")
    relative = Path(checkpoint.get("path", ""))
    if not relative.name or relative.is_absolute():
        raise ValueError("manifest checkpoint path must be relative")
    expected = (manifest_file.parent / relative).resolve()
    if not expected.is_relative_to(manifest_file.parent):
        raise ValueError("manifest checkpoint path escapes the run directory")
    if checkpoint_path is not None and Path(checkpoint_path).resolve() != expected:
        raise ValueError("checkpoint path does not match manifest")
    if not expected.is_file():
        raise FileNotFoundError(expected)
    if sha256_file(expected) != checkpoint.get("sha256"):
        raise ValueError("checkpoint hash does not match manifest")
    return expected


def resume_run(root: str | Path, run_id: str, *, promoted: bool = False) -> RunPaths:
    """仅在 manifest 与 checkpoint 完整匹配时恢复既有 run 上下文。"""

    if not _RUN_ID_RE.fullmatch(run_id):
        raise ValueError("run_id may contain only letters, digits, '.', '_' and '-'")
    paths = RunPaths(run_id=run_id, run_dir=Path(root).resolve() / run_id)
    run_location = paths.promoted_dir if promoted else paths.draft_dir
    manifest_path = run_location / "manifest.json"
    manifest = _load_manifest(manifest_path)
    if manifest.get("run_id") != run_id:
        raise ValueError("manifest run_id does not match run paths")
    validate_checkpoint(manifest_path)
    return paths


def _valid_warning_acceptance(value: Any) -> bool:
    return isinstance(value, Mapping) and all(
        isinstance(value.get(field), str) and value[field].strip()
        for field in ("accepted_by", "accepted_at")
    )


def promote_run(paths: RunPaths, *, warning_acceptance: Mapping[str, str] | None = None) -> Path:
    """验证状态和 checkpoint 后，将整个 draft 原子提升为 promoted。"""

    draft_dir = paths.draft_dir.resolve()
    promoted_dir = paths.promoted_dir.resolve()
    manifest_path = draft_dir / "manifest.json"
    if promoted_dir.exists():
        raise FileExistsError(promoted_dir)
    manifest = _load_manifest(manifest_path)
    if manifest.get("run_id") != paths.run_id:
        raise ValueError("manifest run_id does not match run paths")
    status = StageStatus(manifest.get("stage_status"))
    if status in {StageStatus.FAILED, StageStatus.NEEDS_REVIEW}:
        raise ValueError(f"stage status {status.value} cannot be promoted")
    checkpoint = validate_checkpoint(manifest_path)
    if status is StageStatus.SUCCESS_WITH_WARNINGS:
        acceptance = (
            manifest.get("warning_acceptance") if warning_acceptance is None else warning_acceptance
        )
        if not _valid_warning_acceptance(acceptance):
            raise ValueError("warning acceptance requires accepted_by and accepted_at")
        if manifest.get("warning_acceptance") != acceptance:
            manifest["warning_acceptance"] = dict(acceptance)
            atomic_write_json(manifest_path, manifest, overwrite=True)
    relative_checkpoint = checkpoint.relative_to(draft_dir)
    os.replace(draft_dir, promoted_dir)
    promoted_checkpoint = validate_checkpoint(promoted_dir / "manifest.json")
    if promoted_checkpoint != promoted_dir / relative_checkpoint:
        raise ValueError("promoted checkpoint does not match draft checkpoint")
    return promoted_checkpoint
