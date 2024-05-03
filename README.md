Once all the scripts and landmarks are uploaded, everything should work quite easily!

all_analyses.R will run all the analyses conducted in the study. For the script to work you need:
1) All the landmark files (the .fcsv files inside the all_landmarks folder)
2) A table with that contains each landmark file, the species the file belongs to and it's diet (files_species_diet.csv)
3) A phylogeny of all the species used (phylogeny.tre)

You may need to install packages in R, but these are all listed on the top of the all_analyses.R script

All graphs will be created as .svg file (these are vector files so better scalability) in a folder ./plots

Becuase ggplot does not *yet* work with dispRity, you will have to create a separate.csv file that contains the disparity data (disparity_data.csv) in order to create the graph in lines 696 - 716
