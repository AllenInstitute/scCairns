# Adapting to your data (non-mouse / non-SNS)

The pipeline engine is **dataset-agnostic**. What's mouse- and sympathetic-nervous-
system-specific is a handful of *default values* — and every one is overridable from
your config, with no code changes. This page is the single reference for those
override points, plus a recommendation on whether to fork.

## Every domain-specific default and its override

| What's SNS/mouse-specific | Where the default lives | Override with |
|---|---|---|
| **QC gene patterns** (`mt-`, `rps`/`rpl`, `hb`) — mouse lowercase regex | `qc.*` in the config | `qc.mt_gene_patterns`, `qc.ribo_gene_patterns`, `qc.hb_gene_pattern` — set to your organism's symbol case (e.g. `["^MT-"]` for human) |
| **Marker sets** for the dotplot (`SNS_MARKER_SETS`) | `code/pipeline_config.py`, surfaced as `inspection.marker_sets` in the resolved config | `inspection.markers_json` (a JSON file `{"Group": ["GENE", ...]}`) or `--markers` |
| **Contamination panels** (endothelial/hepatic mouse genes) | `code/flag_contamination.py` `DEFAULT_CONTAM_PANELS` | `inspection.contamination.panels: {label: [gene, ...]}` |
| **Neuronal-fraction heuristic** (`neuronal_markers`, `neuronal_cutoff`) | `inspection.*` in the config | `inspection.neuronal_markers`, `inspection.neuronal_cutoff` (or ignore — it only annotates a "likely neuronal" hint) |
| **`pca_color_gene`** (`Snap25`) | `inspection.pca_color_gene` | any gene in your object |
| **Script name** `integrate_scvi.py` | filename only | cosmetic — the code is generic |

Nothing above changes the model; they change which genes QC, markers, and
contamination scoring look at.

## A worked non-mouse config

Start from **[`examples/pipeline_generic.yml`](../examples/pipeline_generic.yml)**,
which has every knob above marked `# <-- EDIT`. The essentials for, say, a human
dataset:

```yaml
data:
  species: human
  batch_key: sample
  gene_symbol_case: preserve        # or 'upper' to normalize symbols

qc:
  mt_gene_patterns: ["^MT-"]        # human mito genes are MT-*
  ribo_gene_patterns: ["^RPS", "^RPL"]
  hb_gene_pattern: "^HB[^(P)]"

inspection:
  markers_json: my_human_markers.json     # your cell-type panels
  contamination:
    panels:                                # your off-target lineages, or null
      endothelial: [FLT1, KDR, PECAM1, CDH5]
      erythroid:   [HBB, HBA1, HBA2, ALAS2]
```

A markers JSON is just:

```json
{
  "T cells":  ["CD3D", "CD3E", "TRAC"],
  "B cells":  ["MS4A1", "CD79A"],
  "Myeloid":  ["LYZ", "CD14", "FCGR3A"]
}
```

If contamination scoring isn't relevant to your data, you don't have to disable
anything — just don't reference `contamination_flagged_cells.csv` in your
`decisions.yaml`, and ignore the contamination plots.
