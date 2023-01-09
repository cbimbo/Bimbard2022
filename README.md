# Code accompanying Bimbard et al., Nat. Neuroscience, 2023
Code to reproduce the figures from [Bimbard et al., 2023, "Behavioral origin of sound-evoked activity in mouse visual cortex"](https://www.nature.com/articles/s41593-022-01227-x).

Companion data is available on [Figshare](https://doi.org/10.6084/m9.figshare.21371247.v2), and follows the [ONE filename convention](https://int-brain-lab.github.io/iblenv/one_docs/one_reference.html#).

## How to use

The main analyses can be performed using function `audioVisProcessing`. Intermediary files are generated for each analysis.
Figures can be generated using function `audioVisFigures`.

## Dependencies
- rastermap: https://github.com/MouseLand/rastermap
- glmnet-matlab: https://github.com/junyangq/glmnet-matlab
- allenCCF: https://github.com/cortex-lab/allenCCF
