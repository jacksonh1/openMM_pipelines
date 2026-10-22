"""RunLayout — where a run's files live when a separate storage_directory is used.

Model (see `PrepConfig.storage_directory`): `storage_directory` holds the COMPLETE record of a
run; `outdir` is a lightweight working view. Huge raw simulation files are written directly onto
storage and symlinked into `outdir` (`offload`); small files + analysis outputs are written real
in `outdir` and copied into storage (`record` + `mirror`). With `storage_directory=None`
everything is real in `outdir` (today's behaviour).

Run identity is the outdir (path-based): one outdir <-> one storage folder. The folder name is
**deterministic** (so analysis can re-find it and reruns reuse it) and **collision-proof** (a
hash of the absolute outdir path, so two runs sharing a basename never clash). An `O_EXCL`
ownership marker records which outdir owns the folder — enabling orphan detection and a safe,
targeted `clean()`.

Safety-first lifecycle: `offload` refuses to overwrite an existing output (never silently delete
an expensive trajectory); `clean()` wipes a folder only after the marker confirms it belongs to
this outdir. Writers write directly onto storage (never *through* a symlink), so there is no
cross-filesystem `rename()` and no dependence on how a reporter tests file existence.

Pure Python: no OpenMM / mdtraj import, so it is unit-testable anywhere.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import os
import shutil
from dataclasses import dataclass
from pathlib import Path

_MARKER = ".source_outdir"
_HASH_LEN = 8


def _abs(p) -> Path:
    return Path(p).expanduser().resolve()


def read_storage_directory(outdir) -> Path | None:
    """Recover the `storage_directory` a run was produced with, from `outdir/config.json`.

    Absent config.json or a null field => None (a legacy / no-storage run). Analysis uses this to
    rebuild the same `RunLayout` production used. (A missing key means an older run with no
    storage — a semantically meaningful default, not a swallowed error.)
    """
    cfg = Path(outdir) / "config.json"
    if not cfg.exists():
        return None
    sd = json.loads(cfg.read_text()).get("storage_directory")
    return Path(sd) if sd else None


@dataclass(frozen=True)
class RunLayout:
    """Decides symlink-vs-copy for a run's files. Build once per stage from (outdir, storage)."""

    outdir: Path
    storage_directory: Path | None = None

    @property
    def storage_home(self) -> Path | None:
        """The run's folder inside storage: deterministic + collision-proof, or None."""
        if self.storage_directory is None:
            return None
        abs_out = _abs(self.outdir)
        short = hashlib.sha1(str(abs_out).encode()).hexdigest()[:_HASH_LEN]
        return _abs(self.storage_directory) / f"{abs_out.name}-{short}"

    def _ensure_home(self) -> Path:
        """Create storage_home + its O_EXCL ownership marker, or verify an existing marker
        belongs to this outdir (fail loud on a foreign owner — the near-impossible hash clash)."""
        home = self.storage_home
        assert home is not None, "storage_directory is None"
        home.mkdir(parents=True, exist_ok=True)
        marker = home / _MARKER
        abs_out = str(_abs(self.outdir))
        try:
            fd = os.open(str(marker), os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
        except FileExistsError:
            recorded = json.loads(marker.read_text())["abs_outdir"]
            if recorded != abs_out:
                raise RuntimeError(
                    f"storage folder {home} is owned by a different run ({recorded}), not "
                    f"{abs_out}; pass a distinct storage_directory or clean the folder."
                )
            return home
        with os.fdopen(fd, "w") as fh:
            json.dump({"abs_outdir": abs_out,
                       "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
                       "created_utc": _dt.datetime.now(_dt.timezone.utc).isoformat()}, fh)
        return home

    def offload(self, relpath: str) -> Path:
        """Return the path to WRITE a huge raw file to. With storage: the real file lives under
        storage_home and `outdir/relpath` is a symlink to it (the writer writes directly onto
        storage). Without storage: `outdir/relpath`. Fails loud if an output already occupies
        either slot — never silently overwrite an expensive trajectory (use `clean()` to rerun).
        """
        out_path = self.outdir / relpath
        if self.storage_directory is None:
            out_path.parent.mkdir(parents=True, exist_ok=True)
            return out_path

        home = self._ensure_home()
        target = home / relpath
        if target.exists():
            raise RuntimeError(
                f"{target} already exists in storage; refusing to overwrite an existing run "
                f"output. Clean the run (RunLayout.clean) before rerunning."
            )
        target.parent.mkdir(parents=True, exist_ok=True)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        if out_path.is_symlink():
            if os.readlink(out_path) != str(target):
                raise RuntimeError(
                    f"{out_path} already symlinks elsewhere ({os.readlink(out_path)}); "
                    f"expected {target}."
                )
        elif out_path.exists():
            raise RuntimeError(
                f"{out_path} already exists as a real file; refusing to overwrite. "
                f"Clean the run (RunLayout.clean) before rerunning."
            )
        else:
            out_path.symlink_to(target)
        return target

    def record(self, relpath: str) -> Path:
        """Return the real `outdir/relpath` for a small / keep-in-outdir file."""
        out_path = self.outdir / relpath
        out_path.parent.mkdir(parents=True, exist_ok=True)
        return out_path

    def mirror(self, subdir: str) -> None:
        """Copy every non-symlink file under `outdir/subdir` into `storage_home/subdir` so
        storage is the complete record. Skips symlinks (offloaded files already live in storage).
        Idempotent; no-op without storage."""
        if self.storage_directory is None:
            return
        src_root = self.outdir / subdir
        if not src_root.exists():
            return
        home = self._ensure_home()
        for src in src_root.rglob("*"):
            if src.is_symlink() or not src.is_file():
                continue
            dest = home / src.relative_to(self.outdir)
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dest)

    def clean(self) -> None:
        """Delete this run's storage folder (and the `outdir` symlinks into it), but ONLY after
        the ownership marker confirms it belongs to this outdir — the interlock for a safe rerun
        reset. No-op without storage or if the folder is absent."""
        if self.storage_directory is None:
            return
        home = self.storage_home
        if home is None or not home.exists():
            return
        marker = home / _MARKER
        if marker.exists():
            recorded = json.loads(marker.read_text())["abs_outdir"]
            if recorded != str(_abs(self.outdir)):
                raise RuntimeError(
                    f"refusing to clean {home}: owned by {recorded}, not {_abs(self.outdir)}."
                )
        if self.outdir.exists():  # remove the now-dangling outdir symlinks into this folder
            for p in self.outdir.rglob("*"):
                if p.is_symlink() and os.readlink(p).startswith(str(home)):
                    p.unlink()
        shutil.rmtree(home)
