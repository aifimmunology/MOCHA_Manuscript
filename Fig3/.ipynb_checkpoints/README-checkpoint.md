# Figure 3 Readme 

The code and datasets to generate the results for figure are included
and organized in this directory as follows. 

## generate Data
The generate_correlationsAndDifferentials_CD16.R script
generates the tile matrix and differential objects using
for all downstream analyses in Figure3. Given file size limitations, 
the tileMatrix is saved as a tarfile. 

## Panel A: Venn Diagram 

Panel A includes the R script and output graph that shows the Venn Diagram for the R. 

## Panel B: Stability of Analysis 

The Fig3 panelB.R script generates the plots for the N-1 perturbation analyses, 
while the .csv files are the outputs of the N-1 experiments for ArchR, MOCHA, and Signac. 

The remaining R scripts are used to generate the csv outputs for each method. 

## Panel C: RunTime Comparison

The fig3_panelC.R generates the line plots comparing the runtimes for each method. The csv files 
are the runtimes for each method, and the remaining R scripts were used to generate the 
runtimes fpr each method. 

## Panel D, E: Unique Regions

The panelD-D_comparison.R script generates the graphs for panels D and E. 
Panel D displays the median intensity difference against the difference of 0s 
for all unique hits for each method. The panelE_XXX.png show the top 5 hits for each
method, displaying the histograms of intensities across conditions. 

## Panel F-G: Zero Inflated Correlation  

The panelF,G_correlations.R script generates the scatter plots across both panels showing the 
overall agreement between the unadjusted spearman correlation (S) and the Zero-inflated Spearman (ZI-S) across all peaks (panel F) and across individual pairs of peaks (panel G). The ZI_vs_spearman.png is panel (f), and a the remaining .png were used to illustrate examples of disagreements for panel (g). 
