# Quarto paper setup

This directory mirrors the Quarto paper layout used in `eco-ben/ADRIA-reef-indicators`.

## Included

- `paper.qmd`: manuscript source with Elsevier-style Quarto front matter
- `references.bib`: bibliography database
- `global-change-biology.csl`: CSL file used by the reference setup
- `_extensions/`: vendored `authors-block` and Elsevier journal extensions

## Rendering

Install Quarto first, then render from this directory:

```powershell
quarto render paper.qmd
```

Preview while writing:

```powershell
quarto preview paper.qmd
```

## Julia environment

The manuscript is configured with `exeflags: ["--project=.."]` so Julia code chunks use the
study environment from the parent directory, where `Project.toml` and `Manifest.toml` live.

## Notes

- The current manuscript file is a scaffold, not copied prose from the reference paper.
- The output format follows the reference repo's `elsevier-docx` setup.
- If you later want PDF output as well, the vendored Elsevier assets are already present.