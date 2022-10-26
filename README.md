# Code accompanying Bimbard et al., Nat. Neuroscience, 2022
Code to reproduce the figures from Bimbard et al., 2022, "Behavioral origin of sound-evoked activity in mouse visual cortex".

Companion data is available on [Figshare](https://figshare.com/articles/dataset/Dataset_for_Bimbard_et_al_2022/21371247), and follows the [ONE filename convention](https://int-brain-lab.github.io/iblenv/one_docs/one_reference.html#).

## How to use

The main analyses can be performed using function `audioVisProcessing`. Intermediary files are generated for each analysis.
Figures can be generated using function `audioVisFigures`.

## Dependencies
- rastermap: https://github.com/MouseLand/rastermap
- glmnet-matlab: https://github.com/junyangq/glmnet-matlab
- allenCCF: https://github.com/cortex-lab/allenCCF
