# CAD Clustering using Latent Class Analysis

This repository provides all derivation, functions, and validation steps for assigning cardiovascular patient clusters using a Latent Class Analysis (LCA) model derived from the SWEDEHEART registry.

## Structure

- `derivation/`: Cluster derivation in the SWEDEHEART dataset, including model selection and interpretation.
- `functions/`: Fully hardcoded and documented cluster assignment function using pre-specified LCA probabilities.
- `validation/`: Evaluation of function performance on external data (SMART) compared with `poLCA.posterior`.
- `figures/`: Visuals used across HTML outputs (e.g., elbow plot, correlation heatmaps).
- `docs/`: HTML-rendered reports for web access via GitHub Pages.

## GitHub Pages

Rendered HTML outputs are accessible at:  
## 📄 Supplementary Materials (Interactive HTML)

- 🔹 [Derivation of Clusters (SWEDEHEART)](https://jorisholt.github.io/cad-clustering-lca/derivation_clusters_SWEDEHEART.html)  
  Full annotated code for variable preprocessing, model selection, and latent class analysis (LCA) derivation.

- 🔹 [Validation of Custom Function](https://jorisholt.github.io/cad-clustering-lca/Function-clustermembership.html)  
  Demonstrates that a hardcoded implementation reproduces cluster assignments and posterior probabilities from `poLCA.posterior()

> No individual-level data are shared. All results and figures are fully reproducible from the included scripts and functions.
