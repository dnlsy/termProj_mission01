# termProj_mission01
bioinformatics_practice1 (2021 Spring, SNU) TermProject01  

Exercise making **Figure5B** in **Cho, Jun, et al. "LIN28A is a suppressor of ER-associated translation in embryonic stem cells." *Cell* 151.4 (2012): 765-777**.



1. *CoLab_termProj_provided.ipynb*

   Data generation to investigate the correlation between CLIP enrichment and RPF density change in LIN28A-knockdown cells.

   Acquired the transcript count table from from provided data.  

    

2. *mission01.Rmd* and *mission01.html*

   Imported data: the transcript count table and protein localization information for genes

   Removed extremes and adjusted the cut-off for a better visualization

   Merged subcellular location information to the transcript count table    

   Optimized the graph