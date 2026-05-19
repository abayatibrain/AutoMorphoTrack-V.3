# AutoMorphoTrack — Jupyter Notebooks

This repository hosts the **interactive Jupyter-notebook companion** to the
[AutoMorphoTrack](https://github.com/abayatibrain/AutoMorphoTrack) Python
package. The notebooks walk through the same pipeline (detection,
morphology, shape profiling, tracking, motility, colocalization, and the
integrated correlation summary) in a step-by-step, didactic format.

If you want the **installable package** with CLI and the MCP connector for
Claude Code, install from PyPI:

```bash
pip install automorphotrack
```

and see the [package repo](https://github.com/abayatibrain/AutoMorphoTrack)
for full documentation.

## Repository layout

```
notebooks/
  AutoMorphoTrack_Example.ipynb       Minimal end-to-end example
  AutoMorphoTrack_Full_Package.ipynb  Full pipeline run, all stages
  AutoMorphoTrack_Source_Code.ipynb   Cell-by-cell deep dive into each module
docs/
  usage_guide_v3.md                   Notebook-level usage guide
  legacy_README_v3.md                 Historical README for V3.1
```

## Quick start

```bash
# Clone the notebooks
git clone https://github.com/abayatibrain/AutoMorphoTrack-V.3.git
cd AutoMorphoTrack-V.3

# Install the package the notebooks call into
pip install automorphotrack

# Launch Jupyter
jupyter lab notebooks/
```

Open `AutoMorphoTrack_Example.ipynb` and replace the sample TIF path with
your own multichannel stack.

## Relationship to the eLife-reviewed paper

The notebooks here correspond to the "Jupyter notebook" execution mode
described in [Bayati et al., *eLife* 2026](https://elifesciences.org/reviewed-preprints/109936),
DOI [10.7554/eLife.109936.1](https://doi.org/10.7554/eLife.109936.1). The
modular package version and the AMTComparison companion script are tracked
in separate repos:

- **Package** — https://github.com/abayatibrain/AutoMorphoTrack
- **Comparison framework** — https://github.com/abayatibrain/AMTcomparison

## What changed in this cleanup

This repo previously contained a mix of monolithic `.py` modules, `.bak`
backups, `__pycache__`, sample-output PNG/CSVs, and a zip archive of the
whole tree. All of that has been removed; the canonical, maintained code
lives in the **package** repo. Only the notebooks (which are *not*
duplicated in the package repo) remain here so they can be cited and
re-executed as supplementary material.

## License

MIT (see `LICENSE` or the package repo for the full text).
