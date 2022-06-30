#!/usr/bin/env python

# Loop through a bunch of subsets of one dataset, run Structure on all of them for a specified number of Ks

# Before running ths script with `python run-all-structure-plots-again.py` change line 12, 16, 19, and 23.
# The script will automatically run 20 runs of Structure for the Ks you provide with 100k burn in and replications, save outputs for all Ks, and save the Evanno table.

import ipyrad.analysis as ipa
import pandas as pd

# Specify the dataset
data = "/media/tangled/raw_data_backup/20210600_MATE-combined-3rad/88clust-111inds_outfiles/88clust-111inds.snps.hdf5"

# Specify to the metadata subsets that you want to analyze
dfs = ["/media/tangled/raw_data_backup/20210600_MATE-combined-3rad/rad-seq-data-analysis/structure/all-86-inds/metadata-86-inds.csv",  "/media/tangled/raw_data_backup/20210600_MATE-combined-3rad/rad-seq-data-analysis/structure/no-app-suwannee/metadata-no-suwannee-or-appalachicolae.csv", "/media/tangled/raw_data_backup/20210600_MATE-combined-3rad/rad-seq-data-analysis/structure/no-suwannee"]

# Specify a list of names in the same order as dfs
names = ["all-83-inds", "no-suwannee", "no-app-suwannee"]

# Specify a range of Ks following the order as above. Add one on the end because Python is 0 based
k_ranges = [range(1,9), range(1,8), range(1,8)]

# Loop through the structure process
for i in range(len(dfs)):
	df = pd.read_csv(dfs[i], sep = ",")
	print("Now working on", names[i])

	# Create the imap dictionary of populations and sample names
	imap = df.groupby('dummypop')['ind'].apply(list).to_dict()

	# Set parameters
	minmap = {i: 0.01 for i in imap}

	# Set up structure run
	struct = ipa.structure(name = names[i], data = data, imap = imap, minmap = minmap, mincov = 0.8)

	# Set structure parametrs
	struct.mainparams.burnin = 100000
	struct.mainparams.numreps = 100000

	# Set k range
	k = k_ranges[i]

	# Run!
	struct.run(nreps = 20, kpop = k, auto = True)

	print("Done running Structure! Saving files now...")

	# Save Evanno Table
	struct.get_evanno_table(k).to_csv("%s_evanno.csv" % names[i])

	print("Evanno table saved!")

	# Save each Structure output, starting from k = 2, going to max k
	for z in range(2, max(k)+1):
		struct.get_clumpp_table(int(z)).to_csv("%s_structure_resultsk%s.csv" % (names[i], z))
