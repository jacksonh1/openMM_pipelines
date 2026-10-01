# openmm_pipelines

Pythonic OpenMM engine for explicit-solvent MD of **designed protein structures**.
Plain production MD today; T-REMD / REST2 planned behind a shared prepare/produce seam.

## Install

```bash
conda env create -f environment.yml   # creates openmm_env
conda activate openmm_env
pip install -e .
```

## Use

```python
from openmm_pipelines import MDConfig, prepare, md, analysis

cfg = MDConfig(pdb_in="design.pdb", outdir="run1", total_ns=100)
md.run(prepare(cfg), cfg)      # prepare = build + equilibrate; run = production
analysis.analyze(cfg.outdir)   # RMSD/Rg/RMSF/DSSP vs the input pose
```

Outputs land in `outdir/` (`prepared/`, `production/`, `analysis/` + config/run_info).

## Docs

- `SKELETON_v2.md` — design + architecture overview.
- `knowledgebase/` — `CODEBASE_SHAPE.md` (map), `DECISIONS.md` (why), `GOTCHAS.md`, `TODO.md`.
