## Quick guide

To generate manuscript, knit the .Rmd file. This will produce the pdf document and latex files. <p>
	
**Requires** that the full repo is downloaded locally and that the data folder (see main ReadMe) has been added to the top level of the repo folder. <p>

## To obtain data for one task

Run relevant .R file for the given task (e.g. AB.R -> attentional blink)
This file will call efilids_functions.R which contains all the main functions for the bootstrapping analyses (see methods section of manuscript). The outputs of this set of operation (RData files) then need to be moved to a folder called 'IMM[task-name]' which should then be compressed as a zip file. <p>

To obtain the information called from the document .Rmd file (see doc folder), get_statistics.R should be run with the relevant settings adjusted. <p>

## Useful things to know

The bootsrapping code carries two ways of sampling, the method reported in the manuscript and a further method where we sampled a parent N (e.g. 13) and bootstapped 1000 surrogate experiments from that N, and repeated that x 1000 (resulting in 1000^2 samples). This was part of an exploration of the best way to control for reducing homogeneity confounds as N approached N max. However, with exploration and simulations, we found this latter method is innappriopriate - hence most code relating to this method is commented out in the code. <p>

## To generate plots

For each task, plot settings are created using 'create_plot_settings.R', the resulting saved settings are called by the code 'main_figure_per_plot.R' which calls the functions in 'plotting.R'. 'main_figure_per_plot.R' requires some changed settings for each task, and running of select code blocks for each task. Should largely be discernable from the notes/comments. <p>
	
## Other

The remaining, unmentioned files contain functions for performing inferential stats on the full data set, as well as inferential stats on the output of get_statistics.R (the latter of which was ultimately not used). <p>


