## Plot Tree from Tetrad

These scripts allow you to plot the tree created with [tetrad](https://github.com/eaton-lab/tetrad), using the multispecies coalescent.

**Note:** as of 18 July 2022, tetrad requires a diff to install properly. The analyses in our paper can be recreated by:

1. Installing [tetrad](https://github.com/eaton-lab/tetrad)
2. Fixing the [diff](https://github.com/eaton-lab/tetrad/issues/5#issuecomment-873135542) for the working version as of 18 July 2022.
3. [Converting our raw vcf to HDF5](https://ipyrad.readthedocs.io/en/master/API-analysis/cookbook-vcf2hdf5.html)
4. Running tetrad on the HDF5 with 1,000 bootstraps and 3 million quartets to ensure all quartets are sampled.

That should generate a consensus tree similar to the one included here.
