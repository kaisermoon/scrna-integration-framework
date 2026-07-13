"""各阶段共用的运行状态、manifest 与 checkpoint 技术管道。"""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import math
import os
import platform
import re
import stat
import subprocess
import sys
import tempfile
import types
import warnings
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Any

_REDACTED = "<redacted>"
_SENSITIVE_PARAMETER_TOKENS = {"TOKEN", "SECRET", "PASSWORD", "CREDENTIAL", "CREDENTIALS"}
_RESEARCH_KEY_ALLOWLIST = frozenset(
    "BATCH_KEY CELL_TYPE_KEY CLUSTER_KEY DISEASE_KEY DONOR_KEY GROUP_KEY HVG_KEY LABEL_KEY SAMPLE_KEY".split()
)
_R_PACKAGE_RE = re.compile(r"[A-Za-z](?:[A-Za-z0-9.]*[A-Za-z0-9])?\Z")
_VERSION_RE = re.compile(r"[0-9]+(?:[.-][0-9]+)*\Z")


def utc_now_rfc3339() -> str:
    """生成无时区歧义的 UTC 时间，供 manifest 记录事件顺序。"""
    return datetime.now(timezone.utc).isoformat(timespec="microseconds").replace("+00:00", "Z")


def collect_r_environment(
    package_names: Sequence[str] = (), *, rscript: str | Path = "Rscript", timeout: float = 10
) -> dict[str, Any]:
    """显式采集 R 版本；不由 Python runtime provenance 自动触发。"""

    packages = sorted(set(package_names))
    if any(not isinstance(name, str) or _R_PACKAGE_RE.fullmatch(name) is None for name in packages):
        raise ValueError("R package names must be safe identifiers")
    script = (
        "cat('R\\t', as.character(getRversion()), '\\n', sep=''); "
        "for (p in commandArgs(trailingOnly=TRUE)) { "
        "v <- tryCatch(as.character(packageVersion(p)), error=function(e) 'unavailable'); "
        "cat('P\\t', p, '\\t', v, '\\n', sep='') }"
    )
    unavailable = {"available": False, "version": "unavailable", "packages": {name: "unavailable" for name in packages}}
    try:
        result = subprocess.run(
            [os.fspath(rscript), "--vanilla", "-e", script, *packages],
            check=True, capture_output=True, text=True, timeout=timeout,
        )
    except OSError:
        return unavailable | {"error": "Rscript unavailable"}
    except subprocess.TimeoutExpired:
        return unavailable | {"error": "R provenance timed out"}
    except subprocess.SubprocessError:
        return unavailable | {"error": "R provenance command failed"}
    try:
        lines = [line.split("\t") for line in result.stdout.splitlines() if line]
        if not lines or len(lines[0]) != 2 or lines[0][0] != "R":
            raise ValueError
        version = lines[0][1]
        package_lines = lines[1:]
        if (len(package_lines) != len(packages) or
                any(len(fields) != 3 or fields[0] != "P" for fields in package_lines)):
            raise ValueError
        found = {fields[1]: fields[2] for fields in package_lines}
        if (_VERSION_RE.fullmatch(version) is None or len(version) > 128 or
                set(found) != set(packages) or len(found) != len(package_lines) or
                any(value != "unavailable" and _VERSION_RE.fullmatch(value) is None or len(value) > 128
                    for value in found.values())):
            raise ValueError
    except (IndexError, ValueError):
        return unavailable | {"error": "malformed R provenance output"}
    return {"available": True, "version": version, "packages": found}


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
        git_status = subprocess.run(
            ["git", "-C", str(root), "status", "--porcelain=v1", "--untracked-files=all", "-z"],
            check=True,
            capture_output=True,
            timeout=10,
        ).stdout
        git_available = True
        git_dirty: bool | None = bool(git_status)
        git_status_sha256 = hashlib.sha256(git_status).hexdigest() if git_status else None
        git_untracked_count = sum(record.startswith(b"?? ") for record in git_status.split(b"\0"))
        git_tracked_dirty: bool | None = bool(git_diff)
        git_tracked_diff_sha256 = hashlib.sha256(git_diff).hexdigest() if git_diff else None
    except (FileNotFoundError, subprocess.SubprocessError):
        git_available = False
        git_commit = "unavailable"
        git_dirty = None
        git_status_sha256 = None
        git_untracked_count = None
        git_tracked_dirty = None
        git_tracked_diff_sha256 = None

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
        "git_status_sha256": git_status_sha256,
        "git_untracked_count": git_untracked_count,
        "git_tracked_dirty": git_tracked_dirty,
        "git_tracked_diff_sha256": git_tracked_diff_sha256,
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
_SHA256_RE = re.compile(r"[0-9a-fA-F]{64}\Z")
_GIT_COMMIT_RE = re.compile(r"(?:[0-9a-fA-F]{40}|[0-9a-fA-F]{64})\Z")
_RFC3339_UTC_RE = re.compile(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(?:\.\d+)?Z\Z")
MANIFEST_SCHEMA_VERSION = "1"
MAX_MANIFEST_FAILURE_BYTES = 16_384
_RUNTIME_PROVENANCE_KEYS = frozenset("python platform git_available git_commit git_dirty git_status_sha256 git_untracked_count git_tracked_dirty git_tracked_diff_sha256 packages".split())
_MANIFEST_COMMON_KEYS = frozenset("schema_version run_id stage stage_status started_at completed_at inputs effective_parameters runtime_provenance method_status hard_postconditions warnings artifacts".split())


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


def _stat_identity(value: os.stat_result) -> tuple[int, int, int, int, int]:
    return (value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_mode)


def _confined_path(path: str | Path, root: str | Path) -> tuple[Path, str]:
    base = Path(os.path.abspath(os.fspath(root)))
    candidate = Path(path)
    candidate = Path(os.path.abspath(os.fspath(candidate if candidate.is_absolute() else base / candidate)))
    if candidate == base or not candidate.is_relative_to(base):
        raise ValueError("path must identify an entry inside root")
    relative = candidate.relative_to(base).as_posix()
    if "\\" in relative or any(ord(char) < 32 or ord(char) == 127 for char in relative):
        raise ValueError("path must produce a safe relative manifest path")
    cursor = base
    for part in (Path(relative).parts[:-1]):
        info = os.lstat(cursor)
        if stat.S_ISLNK(info.st_mode) or not stat.S_ISDIR(info.st_mode):
            raise ValueError("path contains a symlink or non-directory component")
        cursor /= part
    info = os.lstat(cursor)
    if stat.S_ISLNK(info.st_mode) or not stat.S_ISDIR(info.st_mode):
        raise ValueError("path root contains a symlink or is not a directory")
    return candidate, relative


def _open_parent(root: str | Path, relative: str) -> tuple[int, str]:
    flags = os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(root, flags)
    try:
        for part in Path(relative).parts[:-1]:
            child = os.open(part, flags, dir_fd=descriptor)
            os.close(descriptor)
            descriptor = child
    except Exception:
        os.close(descriptor)
        raise
    return descriptor, Path(relative).parts[-1]


def _bound_file_hash_at(parent: int, name: str) -> tuple[tuple[int, ...], int, str]:
    before = os.stat(name, dir_fd=parent, follow_symlinks=False)
    if stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode):
        raise ValueError("fingerprinted path must be a regular non-symlink file")
    flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(name, flags, dir_fd=parent)
    try:
        opened = os.fstat(descriptor)
        digest = hashlib.sha256()
        with os.fdopen(descriptor, "rb", closefd=False) as stream:
            while chunk := stream.read(1024 * 1024):
                digest.update(chunk)
        after_open = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    after = os.stat(name, dir_fd=parent, follow_symlinks=False)
    identities = {_stat_identity(item) for item in (before, opened, after_open, after)}
    if len(identities) != 1:
        raise RuntimeError("file changed while fingerprinting")
    return _stat_identity(before), before.st_size, digest.hexdigest()


def _scan_tree(root: int) -> list[tuple[str, str, int, str]]:
    root_info = os.fstat(root)
    if stat.S_ISLNK(root_info.st_mode) or not stat.S_ISDIR(root_info.st_mode):
        raise ValueError("tree root must be a non-symlink directory")
    records: list[tuple[str, str, int, str]] = []
    seen: set[str] = set()

    def visit(directory: int, prefix: str = "") -> None:
        with os.scandir(directory) as entries:
            ordered = sorted(entries, key=lambda item: item.name)
        for entry in ordered:
            info = entry.stat(follow_symlinks=False)
            relative = f"{prefix}/{entry.name}" if prefix else entry.name
            folded = relative.casefold()
            if folded in seen:
                raise ValueError("tree contains duplicate or case-colliding paths")
            seen.add(folded)
            if stat.S_ISLNK(info.st_mode):
                raise ValueError("tree must not contain symlinks")
            if stat.S_ISDIR(info.st_mode):
                records.append(("directory", relative, 0, ""))
                child = os.open(
                    entry.name,
                    os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0),
                    dir_fd=directory,
                )
                try:
                    if _stat_identity(os.fstat(child)) != _stat_identity(info):
                        raise RuntimeError("tree changed while fingerprinting")
                    visit(child, relative)
                finally:
                    os.close(child)
            elif stat.S_ISREG(info.st_mode):
                _, size, checksum = _bound_file_hash_at(directory, entry.name)
                records.append(("file", relative, size, checksum))
            else:
                raise ValueError("tree must contain only directories and regular files")

    visit(root)
    return sorted(records, key=lambda record: record[1])


def fingerprint_input(
    role: str, path: str | Path, root: str | Path, *, kind: str = "file"
) -> dict[str, str]:
    """为输入生成稳定记录，并在计算期间绑定路径与内容。"""

    if not isinstance(role, str) or _RUN_ID_RE.fullmatch(role) is None:
        raise ValueError("role must be a safe identifier")
    _, relative = _confined_path(path, root)
    if kind == "file":
        parent, name = _open_parent(root, relative)
        try:
            identity, _, checksum = _bound_file_hash_at(parent, name)
        finally:
            os.close(parent)
    elif kind == "tree":
        parent, name = _open_parent(root, relative)
        try:
            tree = os.open(name, os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0),
                           dir_fd=parent)
        finally:
            os.close(parent)
        try:
            identity = _stat_identity(os.fstat(tree))
            first, second = _scan_tree(tree), _scan_tree(tree)
            if first != second:
                raise RuntimeError("tree changed while fingerprinting")
            canonical = json.dumps(first, ensure_ascii=False, separators=(",", ":")).encode()
            checksum = hashlib.sha256(canonical).hexdigest()
        finally:
            os.close(tree)
    else:
        raise ValueError("kind must be file or tree")
    verify_parent, verify_name = _open_parent(root, relative)
    try:
        verify = os.stat(verify_name, dir_fd=verify_parent, follow_symlinks=False)
    finally:
        os.close(verify_parent)
    if _stat_identity(verify) != identity:
        raise RuntimeError("path changed while fingerprinting")
    return {"role": role, "path": relative, "kind": kind, "sha256": checksum}


def artifact_record(role: str, path: str | Path, state_dir: str | Path) -> dict[str, str | int]:
    """为 run state 内的常规文件生成可校验的 artifact 记录。"""

    if not isinstance(role, str) or _RUN_ID_RE.fullmatch(role) is None:
        raise ValueError("role must be a safe identifier")
    _, relative = _confined_path(path, state_dir)
    parent, name = _open_parent(state_dir, relative)
    try:
        identity, size, checksum = _bound_file_hash_at(parent, name)
    finally:
        os.close(parent)
    verify_parent, verify_name = _open_parent(state_dir, relative)
    try:
        verify = os.stat(verify_name, dir_fd=verify_parent, follow_symlinks=False)
    finally:
        os.close(verify_parent)
    if _stat_identity(verify) != identity:
        raise RuntimeError("artifact path changed while fingerprinting")
    return {"role": role, "path": relative, "size": size, "sha256": checksum}


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


def _reject_run_state_symlinks(
    run_dir: Path, state_dir: Path, manifest_file: Path | None = None
) -> None:
    checks = (("run root", run_dir.parent), ("run directory", run_dir), ("run state", state_dir))
    if manifest_file is not None:
        checks += (("manifest", manifest_file),)
    for label, path in checks:
        if path.is_symlink():
            raise ValueError(f"{label} path must not be a symlink: {path}")


def _resolve_manifest_file(manifest_path: str | Path) -> Path:
    candidate = Path(manifest_path)
    state_dir = candidate.parent
    _reject_run_state_symlinks(state_dir.parent, state_dir, candidate)
    return candidate.resolve()


def _valid_warning_acceptance(value: Any) -> bool:
    return isinstance(value, Mapping) and value.keys() == {"accepted_by", "accepted_at"} and all(isinstance(value.get(field), str) and value[field].strip() for field in ("accepted_by", "accepted_at"))


def _required_container(manifest: Mapping[str, Any], key: str, expected: type) -> Any:
    value = manifest.get(key)
    if not isinstance(value, expected):
        raise ValueError(f"manifest {key} must be a {expected.__name__}")
    return value


def _validate_json_value(value: Any, label: str) -> None:
    if value is None or isinstance(value, (str, int, bool)):
        return
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError(f"{label} must not contain NaN or infinity")
    elif isinstance(value, list):
        for index, item in enumerate(value):
            _validate_json_value(item, f"{label}[{index}]")
    elif isinstance(value, Mapping):
        if any(not isinstance(key, str) for key in value):
            raise ValueError(f"{label} object keys must be strings")
        for key, item in value.items():
            _validate_json_value(item, f"{label}.{key}")
    else:
        raise ValueError(f"{label} contains non-JSON value {type(value).__name__}")


def _safe_relative_path(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value or "\\" in value:
        raise ValueError(f"{label} must be a safe relative path")
    path = Path(value)
    if (not path.name or path.is_absolute() or ".." in path.parts or re.match(r"[A-Za-z]:", value) or
            any(ord(char) < 32 or ord(char) == 127 for char in value)):
        raise ValueError(f"{label} must be a safe relative path")
    return path.as_posix()


def _validate_records(records: list[Any], label: str, *, require_kind: bool) -> set[str]:
    roles, paths = set(), set()
    for index, record in enumerate(records):
        if not isinstance(record, Mapping):
            raise ValueError(f"manifest {label}[{index}] must be an object")
        if record.keys() != ({"role", "path", "sha256", "kind"} if require_kind else {"role", "path", "size", "sha256"}):
            raise ValueError(f"manifest {label}[{index}] has missing or unknown keys")
        role, checksum = record.get("role"), record.get("sha256")
        if not isinstance(role, str) or _RUN_ID_RE.fullmatch(role) is None or role in roles:
            raise ValueError(f"manifest {label} roles must be unique safe identifiers")
        path = _safe_relative_path(record.get("path"), f"manifest {label}[{index}].path")
        if path in paths:
            raise ValueError(f"manifest {label} paths must be unique")
        if require_kind and record.get("kind") not in {"file", "tree"}:
            raise ValueError(f"manifest {label}[{index}].kind must be file or tree")
        if not require_kind and (type(record.get("size")) is not int or record["size"] < 0):
            raise ValueError(f"manifest {label}[{index}].size must be a non-negative integer")
        if not isinstance(checksum, str) or _SHA256_RE.fullmatch(checksum) is None:
            raise ValueError(f"manifest {label}[{index}].sha256 must be 64 hexadecimal characters")
        roles.add(role)
        paths.add(path)
    return paths


def _validate_runtime_provenance(runtime: dict[str, Any]) -> None:
    if runtime.keys() not in (_RUNTIME_PROVENANCE_KEYS, _RUNTIME_PROVENANCE_KEYS | {"r_environment"}):
        raise ValueError("manifest runtime_provenance must use the exact collector keys")
    if (not all(isinstance(runtime[key], str) and runtime[key] for key in ("python", "platform", "git_commit")) or
            type(runtime["git_available"]) is not bool or not isinstance(runtime["packages"], dict)):
        raise ValueError("manifest runtime_provenance has invalid minimum field types")
    packages = runtime["packages"]
    if any(not isinstance(key, str) or not key or not isinstance(value, str) or not value for key, value in packages.items()):
        raise ValueError("manifest runtime_provenance packages must map non-empty strings to strings")
    if "r_environment" in runtime:
        environment = runtime["r_environment"]
        available = environment.get("available") if isinstance(environment, dict) else None
        packages = environment.get("packages") if isinstance(environment, dict) else None
        common = (type(available) is bool and isinstance(packages, dict) and
                  all(isinstance(key, str) and _R_PACKAGE_RE.fullmatch(key) and isinstance(value, str)
                      for key, value in packages.items()))
        if available is True:
            valid = (environment.keys() == {"available", "version", "packages"} and common and
                     isinstance(environment.get("version"), str) and
                     _VERSION_RE.fullmatch(environment["version"]) is not None and
                     all(value == "unavailable" or _VERSION_RE.fullmatch(value) for value in packages.values()))
        else:
            error = environment.get("error") if isinstance(environment, dict) else None
            valid = (isinstance(environment, dict) and environment.keys() == {"available", "version", "packages", "error"} and common and
                     available is False and environment.get("version") == "unavailable" and
                     all(value == "unavailable" for value in packages.values()) and
                     isinstance(error, str) and error.strip() and len(error) <= 256 and
                     not any(ord(char) < 32 for char in error))
        if not valid:
            raise ValueError("manifest r_environment has invalid availability fields")
    optional = ("git_dirty", "git_tracked_dirty")
    if any(value is not None and type(value) is not bool for value in (runtime[key] for key in optional)):
        raise ValueError("manifest runtime_provenance dirty fields must be booleans or null")
    count = runtime["git_untracked_count"]
    if count is not None and (type(count) is not int or count < 0):
        raise ValueError("manifest git_untracked_count must be a non-negative integer or null")
    if runtime["git_available"]:
        if (_GIT_COMMIT_RE.fullmatch(runtime["git_commit"]) is None or count is None or
                any(runtime[key] is None for key in optional)):
            raise ValueError("available Git provenance requires commit, booleans, and count")
        for dirty_key, hash_key in (("git_dirty", "git_status_sha256"),
                                    ("git_tracked_dirty", "git_tracked_diff_sha256")):
            digest = runtime[hash_key]
            if (runtime[dirty_key] and (not isinstance(digest, str) or _SHA256_RE.fullmatch(digest) is None) or
                    not runtime[dirty_key] and digest is not None):
                raise ValueError(f"manifest {hash_key} does not match {dirty_key}")
    elif (runtime["git_commit"] != "unavailable" or count is not None or
          any(runtime[key] is not None for key in (*optional, "git_status_sha256", "git_tracked_diff_sha256"))):
        raise ValueError("unavailable Git provenance must use unavailable/null fields")


def _utc_timestamp(manifest: Mapping[str, Any], key: str) -> datetime:
    value = manifest.get(key)
    if not isinstance(value, str) or _RFC3339_UTC_RE.fullmatch(value) is None:
        raise ValueError(f"manifest {key} must be an RFC3339 UTC-Z timestamp")
    return datetime.fromisoformat(value[:-1] + "+00:00")


def validate_manifest(manifest_path_or_mapping: str | Path | Mapping[str, Any], *, allow_legacy: bool = False,
                      state: str | None = None, expected_run_id: str | None = None) -> dict[str, Any]:
    """严格校验 v1 manifest；legacy 仅能经显式临时兼容开关进入。"""

    file_context = not isinstance(manifest_path_or_mapping, Mapping)
    if not file_context:
        manifest = dict(manifest_path_or_mapping)
        if state is None:
            raise ValueError("state is required when validating a manifest mapping")
    else:
        manifest_file = _resolve_manifest_file(manifest_path_or_mapping)
        if manifest_file.name != "manifest.json" or manifest_file.parent.name not in {"draft", "promoted"}:
            raise ValueError("manifest path must be <run_id>/<draft|promoted>/manifest.json")
        inferred_state, inferred_run_id = manifest_file.parent.name, manifest_file.parent.parent.name
        if state is not None and state != inferred_state:
            raise ValueError("manifest state does not match path context")
        if expected_run_id is not None and expected_run_id != inferred_run_id:
            raise ValueError("expected run_id does not match path context")
        state, expected_run_id = inferred_state, inferred_run_id
        manifest = _load_manifest(manifest_file)
    if state not in {"draft", "promoted"}:
        raise ValueError("manifest state must be draft or promoted")
    if "schema_version" not in manifest:
        if allow_legacy:
            warnings.warn("legacy manifest compatibility is temporary; migrate to schema v1",
                          DeprecationWarning, stacklevel=2)
            return manifest
        raise ValueError("manifest schema_version is required")
    if manifest.get("schema_version") != MANIFEST_SCHEMA_VERSION:
        raise ValueError(f"manifest schema_version must be {MANIFEST_SCHEMA_VERSION!r}")
    _validate_json_value(manifest, "manifest")
    run_id, stage = manifest.get("run_id"), manifest.get("stage")
    if not isinstance(run_id, str) or _RUN_ID_RE.fullmatch(run_id) is None:
        raise ValueError("manifest run_id must be a safe ASCII identifier")
    if expected_run_id is not None and run_id != expected_run_id:
        raise ValueError("manifest run_id does not match path context")
    if not isinstance(stage, str) or _RUN_ID_RE.fullmatch(stage) is None:
        raise ValueError("manifest stage must be a safe ASCII identifier")
    status_value = manifest.get("stage_status")
    try:
        status = StageStatus(status_value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"invalid manifest stage_status: {status_value!r}") from error
    required_keys = _MANIFEST_COMMON_KEYS | ({"failure"} if status is StageStatus.FAILED else {"checkpoint"})
    allowed_keys = required_keys | ({"warning_acceptance"} if status is StageStatus.SUCCESS_WITH_WARNINGS else set())
    if not required_keys <= manifest.keys() or manifest.keys() - allowed_keys:
        raise ValueError("manifest has missing, unknown, or status-forbidden top-level keys")
    started, completed = _utc_timestamp(manifest, "started_at"), _utc_timestamp(manifest, "completed_at")
    if completed < started:
        raise ValueError("manifest completed_at must not precede started_at")
    inputs = _required_container(manifest, "inputs", list)
    if not inputs:
        raise ValueError("manifest inputs must be non-empty")
    _validate_records(inputs, "inputs", require_kind=True)
    _required_container(manifest, "effective_parameters", dict)
    runtime = _required_container(manifest, "runtime_provenance", dict)
    _validate_runtime_provenance(runtime)
    methods = _required_container(manifest, "method_status", dict)
    if not methods or any(not isinstance(key, str) or not key or not isinstance(value, str) or not value
                          for key, value in methods.items()):
        raise ValueError("manifest method_status must be a non-empty string mapping")
    try:
        [MethodStatus(value) for value in methods.values()]
    except ValueError as error:
        raise ValueError("manifest method_status values must use MethodStatus") from error
    postconditions = _required_container(manifest, "hard_postconditions", dict)
    if not postconditions or any(not isinstance(key, str) or not key or not isinstance(value, bool)
                                 for key, value in postconditions.items()):
        raise ValueError("manifest hard_postconditions must be a non-empty boolean mapping")
    warnings_value = _required_container(manifest, "warnings", list)
    if (any(not isinstance(item, str) or not item.strip() for item in warnings_value) or
            len(warnings_value) != len(set(warnings_value))):
        raise ValueError("manifest warnings must contain unique non-empty strings")
    artifacts = _required_container(manifest, "artifacts", list)
    artifact_paths = _validate_records(artifacts, "artifacts", require_kind=False)
    if status is StageStatus.FAILED:
        if state != "draft":
            raise ValueError("FAILED manifest must remain draft")
        if "checkpoint" in manifest:
            raise ValueError("FAILED manifest must not declare a checkpoint")
        failure = manifest.get("failure")
        if not isinstance(failure, Mapping) or not failure:
            raise ValueError("FAILED manifest requires a non-empty failure object")
        try:
            failure_size = len(json.dumps(failure, ensure_ascii=False, allow_nan=False).encode())
        except (TypeError, ValueError) as error:
            raise ValueError("manifest failure must be finite JSON data") from error
        if failure_size > MAX_MANIFEST_FAILURE_BYTES:
            raise ValueError("manifest failure exceeds the size limit")
    else:
        checkpoint = manifest.get("checkpoint")
        if not isinstance(checkpoint, Mapping) or checkpoint.keys() != {"path", "sha256"}:
            raise ValueError("non-FAILED manifest requires a non-empty checkpoint object")
        checkpoint_path = _safe_relative_path(checkpoint.get("path"), "manifest checkpoint.path")
        checksum = checkpoint.get("sha256")
        if not isinstance(checksum, str) or _SHA256_RE.fullmatch(checksum) is None:
            raise ValueError("manifest checkpoint.sha256 must be 64 hexadecimal characters")
        if checkpoint_path in artifact_paths:
            raise ValueError("manifest checkpoint and artifact paths must be disjoint")
    if status is StageStatus.NEEDS_REVIEW and state != "draft":
        raise ValueError("NEEDS_REVIEW manifest must remain draft")
    acceptance = manifest.get("warning_acceptance")
    if acceptance is not None and not _valid_warning_acceptance(acceptance):
        raise ValueError("warning_acceptance requires accepted_by and accepted_at")
    if status is StageStatus.SUCCESS:
        if warnings_value or acceptance is not None:
            raise ValueError("SUCCESS manifest must have no warnings or warning_acceptance")
    elif status is StageStatus.SUCCESS_WITH_WARNINGS:
        if not warnings_value:
            raise ValueError("SUCCESS_WITH_WARNINGS requires non-empty warnings")
        if state == "promoted" and not _valid_warning_acceptance(acceptance):
            raise ValueError("promoted warnings require warning_acceptance")
    elif acceptance is not None:
        raise ValueError(f"{status.value} manifest must not declare warning_acceptance")
    if acceptance is not None:
        accepted_at = _utc_timestamp(acceptance, "accepted_at")
        if accepted_at < completed:
            raise ValueError("warning accepted_at must not precede completed_at")
    return manifest


def _validate_declared_artifacts(
    manifest_file: Path, manifest: Mapping[str, Any]
) -> dict[str, Path]:
    if "artifacts" not in manifest:
        return {}
    artifacts = manifest["artifacts"]
    if not isinstance(artifacts, list) or not artifacts:
        raise ValueError("manifest artifacts must be a non-empty list when present")
    state_dir = manifest_file.parent.resolve()
    validated: dict[str, Path] = {}
    seen_paths: set[Path] = set()
    for index, artifact in enumerate(artifacts):
        if not isinstance(artifact, Mapping):
            raise ValueError(f"artifact {index} must be an object")
        role = artifact.get("role")
        if not isinstance(role, str) or _RUN_ID_RE.fullmatch(role) is None:
            raise ValueError(f"artifact {index} role must be a non-empty safe name")
        if role in validated:
            raise ValueError(f"duplicate artifact role: {role}")
        raw_path = artifact.get("path")
        if not isinstance(raw_path, str) or not raw_path.strip():
            raise ValueError(f"artifact {role} path must be a non-empty relative path")
        relative = Path(raw_path)
        if not relative.name or any(ord(char) < 32 or ord(char) == 127 for char in raw_path):
            raise ValueError(f"artifact {role} path must be a safe relative file path")
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"artifact {role} path must stay inside the run state directory")
        lexical = manifest_file.parent / relative
        resolved = lexical.resolve()
        if not resolved.is_relative_to(state_dir):
            raise ValueError(f"artifact {role} path escapes the run state directory")
        cursor = manifest_file.parent
        for part in relative.parts:
            cursor /= part
            if cursor.is_symlink():
                raise ValueError(f"artifact {role} path contains a symlink")
        if resolved in seen_paths:
            raise ValueError(f"duplicate artifact path: {relative.as_posix()}")
        if not resolved.exists():
            raise FileNotFoundError(resolved)
        if not resolved.is_file():
            raise ValueError(f"artifact {role} must be a regular file")
        expected_hash = artifact.get("sha256")
        if not isinstance(expected_hash, str) or _SHA256_RE.fullmatch(expected_hash) is None:
            raise ValueError(f"artifact {role} sha256 must be 64 hexadecimal characters")
        if sha256_file(resolved) != expected_hash.lower():
            raise ValueError(f"artifact {role} hash does not match manifest")
        validated[role] = resolved
        seen_paths.add(resolved)
    return validated


def validate_artifacts(manifest_path: str | Path) -> dict[str, Path]:
    """校验 manifest 声明的全部附属产物，返回规范化 role 到路径映射。"""

    manifest_file = _resolve_manifest_file(manifest_path)
    return _validate_declared_artifacts(manifest_file, _load_manifest(manifest_file))


def validate_checkpoint(manifest_path: str | Path, checkpoint_path: str | Path | None = None) -> Path:
    """校验 manifest 声明的 checkpoint 路径和内容 hash，可用于安全恢复。"""

    manifest_file = _resolve_manifest_file(manifest_path)
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
    _validate_declared_artifacts(manifest_file, manifest)
    return expected


def resume_run(root: str | Path, run_id: str, *, promoted: bool = False) -> RunPaths:
    """仅在 manifest 与 checkpoint 完整匹配时恢复既有 run 上下文。"""

    if not _RUN_ID_RE.fullmatch(run_id):
        raise ValueError("run_id may contain only letters, digits, '.', '_' and '-'")
    root_path = Path(root)
    if root_path.is_symlink():
        raise ValueError(f"run root path must not be a symlink: {root_path}")
    paths = RunPaths(run_id=run_id, run_dir=root_path.resolve() / run_id)
    run_location = paths.promoted_dir if promoted else paths.draft_dir
    manifest_path = _resolve_manifest_file(run_location / "manifest.json")
    manifest = _load_manifest(manifest_path)
    if manifest.get("run_id") != run_id:
        raise ValueError("manifest run_id does not match run paths")
    validate_checkpoint(manifest_path)
    return paths


def promote_run(paths: RunPaths, *, warning_acceptance: Mapping[str, str] | None = None) -> Path:
    """验证状态和 checkpoint 后，将整个 draft 原子提升为 promoted。"""

    draft_dir = paths.draft_dir.resolve()
    promoted_dir = paths.promoted_dir.resolve()
    _reject_run_state_symlinks(paths.run_dir, paths.draft_dir, paths.draft_dir / "manifest.json")
    if paths.promoted_dir.is_symlink():
        raise ValueError(f"run state path must not be a symlink: {paths.promoted_dir}")
    manifest_path = _resolve_manifest_file(draft_dir / "manifest.json")
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
    try:
        promoted_checkpoint = validate_checkpoint(promoted_dir / "manifest.json")
        if promoted_checkpoint != promoted_dir / relative_checkpoint:
            raise ValueError("promoted checkpoint does not match draft checkpoint")
    except Exception as original_error:
        try:
            if draft_dir.exists() or draft_dir.is_symlink():
                raise FileExistsError(f"cannot roll back over existing draft: {draft_dir}")
            if promoted_dir.is_symlink() or not promoted_dir.is_dir():
                raise ValueError(f"cannot roll back unsafe promoted state: {promoted_dir}")
            os.replace(promoted_dir, draft_dir)
        except Exception as rollback_error:
            raise RuntimeError(
                f"post-promotion validation failed ({original_error}); "
                f"rollback failed ({rollback_error})"
            ) from original_error
        raise
    return promoted_checkpoint
