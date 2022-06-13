# MetabolismManuscript
R code for Metabolism Manuscript.

## Scripts

### Most important R scripts:
    1.) global_metab_maps.R
    2.) micnm_all.R
    3.) timecourse_all.R
    4.) timecourse_clustering_metabolism.R
    5.) timecourse_eggscaling.R

### General Description of Scripts:
Script #1: 
<br />
transforms all the metabolomics data (mic/nm, mic/nm inhibitors, metabolomic timecourse) and generates KEGG global maps used in the main and supplemental figures of the paper
<br />
Script #2: 
<br />
used the mic/nm proteomics dataset to generate a PCA plot, Heatmap of Significant Proteins, and pathway analysis plots (from Over-representation analysis and GSEA)
<br />
Script #3: 
<br />
uses the timecourse proteomics dataset to generate a PCA plot, Heatmap of Signficant Proteins, Fuzzy-c clustering analysis of significant proteins (enzymes only and all signfiicant proteins), and Pathway Analysis of clusters presented in fuzzy-c clustering analysis.
<br />
Script #4: 
<br />
uses the timecourse metabolomics dataset to generate a fuzzy-c clustering analysis of the top 50 metabolites int he dataset. Pathway Analysis involved using clusters presented in fuzzy-c clustering analysis.
<br />
Script #5: 
<br />
uses the timecourse proteomics and metabolomics data and performs scaling to egg expression levels to output files to visualize on the timecourse app.

#### One Important Remark about Script #3 and Script #4:
<br />
There is a possibiliy of reaching a neighboring local minimum instead of a global minimum when doing fuzzy c-means clustering. If this occurs the R script will output a self-imposed error of "You detected a local minimum, retry the calculation". Retry the calculation until this error doesn't appear. Additionally a histogram with a plot of the frequency obtaining the local vs approx. global minimum will appear.

#### Important Note about Script #1:
<br />
Running the pathview (global metabolism) plots to obtain outputs in the "/Output/metab_global/Figures" folder will take approximately 2-3 hours per figure.
    
## Outputs

### data_files_for_app
    Contains csv files that the mic/nm app and timecourse app requires to run.

### metab_global
    shows 2 folders:
        1.) Figures folder contains all the outputs of the global metabolism maps using the pathview package in script #1.
        2.) Files folder contains ID conversion tables of various datasets and includes a name_map_all_datasets_manu.csv to show every metabolite conversion that was made using the Metaboanalyst online server + manual annotation.

### metab_timecourse_figures
    shows plots of: 
        1.) Fuzzy c-means clustering of timecourse metabolomics data.
        2.) Overrepresentation analysis of 5 clusters found using KEGG annotations

### micnm_supplementary_figures
    shows plots of:
        1.) PCA of the mic/nm data
        2.) Heatmap of the mic/nm data
        3.) GSEA of mic/nm data using FUNC, FUNC_ver2, GO, and KEGG annotations
        4.) Overrepresentation Analysis of mic/nm data using FUNC, FUNC_ver2, GO, and KEGG annotations

### timecourse_supplementary_figures
    shows plots of:
        1.) PCA of timecourse data
        2.) Heatmap of timecourse data
        3.) Fuzzy c-means clustering of all timecourse data
        4.) Fuzzy c-means clustering of only annotated enzymes in the timecourse data
        5.) Overrepresentation analysis of the fuzzy c-means clustering done   in (3) or (4) using FUNC, GO, and KEGG annotations.

More specific details related to each plot are located within the methods section of the manuscript.
        







