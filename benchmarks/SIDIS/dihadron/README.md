# Dihadron correlation benchmark

Multi-trigger dihadron C(Δφ) per x_Bj decade (`analysis/ePIC_multikt2.cxx`), run on
campaign RECO files. There is no simulation or reconstruction step: the input is a list
of RECO file URLs. Moved from https://github.com/rahmans1/DihadronAnalysis (c36b8c6);
the macro is unchanged apart from one comment.

| Stage | Target | Does |
|---|---|---|
| map | `results/dihadron/<sample>/<mode>/job_<job>.root` | Runs the macro on `results/dihadron/<sample>/lists/list_<job>.txt` |
| merge | `results/dihadron/<sample>/<mode>/merged.root` | hadd of every listed job |
| finalize | `results/dihadron/<sample>/<mode>/plots` | One PNG per x decade with triggers, `hists.json` (bin contents of every histogram), tests |

`<sample>` starts with the collision system (`ep`, `eAu` or `eAg`), which selects the
theory file for the saturation weight. `<mode>` is `truth` or `reco`.

## Run locally

Inside eic-shell, on the nine-file list in `lists/` (26.07.3 ep 9×130, x 10⁻²–10⁻¹):

```
snakemake --cores 1 dihadron_run_locally
```

The theory files `Cth_ep.root`, `Cth_eAu.root`, `Cth_eAg.root` (about 55 MB each) are
kept with Git LFS in [eic/epic-data](https://github.com/eic/epic-data)
(`physics_benchmarks/SIDIS/dihadron`). The benchmark downloads the ones a sample needs
at a fixed commit and checks each sha256. To use a local copy instead, add
`--config DIHADRON_CTH_DIR=/path/to/TH`.

## Production

Running over a whole campaign sample is production-scale work, so it is not done in CI.
CI compiles the macro and runs the nine-file smoke test above. epicprod runs whole
samples from the card.

`card.yml` describes the three stages for a production system: it writes each job's
file list, runs the map target per job, then merge and finalize once.
