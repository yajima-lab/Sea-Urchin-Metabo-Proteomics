# MetabolismManuscript
R code for Metabolism Manuscript.

Most important R scripts:
1.) global_metab_maps.R
2.) micnm_all.R
3.) timecourse_all.R
4.) timecourse_clustering_metabolism.R

Script #1 transforms all the metabolomics data (mic/nm, mic/nm inhibitors, metabolomic timecourse) and generates KEGG global maps used in the main and supplemental figures of the paper

Script #2 used the mic/nm proteomics dataset to generate a PCA plot, Heatmap of Significant Proteins, and pathway analysis plots (from Over-representation analysis and GSEA)

Script #3 uses the timecourse proteomics dataset to generate a PCA plot, Heatmap of Signficant Proteins, Fuzzy-c clustering analysis of significant proteins (enzymes only and all signfiicant proteins), and Pathway Analysis of clusters presented in fuzzy-c clustering analysis.

Script #4 uses the timecourse metabolomics dataset to generate a fuzzy-c clustering analysis of the top 50 metabolites int he dataset. Pathway Analysis involved using clusters presented in fuzzy-c clustering analysis.

Notes for further use:
Add explanation of certain parts of Mfuzz code.
Add the location of output to the corresponding supplementary figure.
edit readme for runtimes.


