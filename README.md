This project features the comparison of H. sapiens (human) and S. purpuratus (sea urchin) scRNA-seq datasets using the SAMap algorithm. In our analysis, we used the publicly available scRNA-seq datasets for the adult S. purpuratus (GSE246430), fetal H. sapiens ovaries (GSE118127), (GSE143380), and adult H. sapiens ovaries. Our analysis indicated strong elements of conservation between muscle and immune cell types and ortholog gene pairs with overlapping expression.

We used the S. purpuratus V5 transcriptome (from https://www.echinobase.org/echinobase/displayJBrowse.do?data=data/sp5_0) and the human v47 transcriptome (from https://www.gencodegenes.org/human/) for the reciprocal blast. Details for the reciprocal blast can be found at https://github.com/atarashansky/SAMap. Our code for changing the gene names of the Sp and Hs datasets to be compatible with the scRNA-seq datasets is in hs_sp_name_converison.R. 

The code for the individual processing for the Hs and Sp samples are located in hs_scrnaseq.R and sp_scrnaseq.R. The saved h5ad files were used for subsequent SAMap analysis.

The SAMap analysis.py file contains the code for the cross-species SAMap analysis comparing the individual scRNA-seq datasets. This code generates overlapping UMAP plots, mapping tables with SAMap scores, and top gene pairs in each cluster mapping. Additional information on how SAMap objects can be constructed and processed is available in the SAMap vignette: (https://github.com/atarashansky/SAMap/blob/main/SAMap_vignette.ipynb) 

The SAMap_vis.R generates sankey plots from an input folder containing mapping tables. The sp_integration contains for our comparison of the Sp adult ovary with Sp embryonic cells. 






