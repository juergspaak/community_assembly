Computer code from
"Continuous assembly required: perpetual species turnover in two trophic level ecosystems"
J. W. Spaak, P. B. Adler and S. P. Ellner

Functions
assembly_functions.py
	Contains the functions to create the species traits, LV community parameters and run the community assembly

convert_biotime.py
	Loads biotime dataset and converts and compute the presence absence matrices for each dataset
	This script loads two datasets which are not part of this repository, rather these can be downloaded
	directly from biotime (the files are too large for a git repository)
	Generates the file "biotime_converted.npz" with the presence absence matrices

functions_for_plotting.py
	Contains functions to plot the various figures

various_competition_kernels.py
	Contains the functions for competition with non-gaussian competition kernels

Figure plotting
plot_*.py
	Creates the figure Figure_*.pdf