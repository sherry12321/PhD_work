fGSEA <- function(input_ranking, geneset, title, outbase){
    library(fgsea)
    library(msigdbr)
    library(data.table)
    # set up the rankings in an appropriate format
    df_gsea_sorted0 = read.csv(input_ranking)
    #print(head(df_gsea_sorted0))
    df_gsea_sorted = df_gsea_sorted0$score
    names(df_gsea_sorted) = df_gsea_sorted0$primerid
            
    # read the gene set into an apprpriate format
    gset = gage::readList(geneset)
    #print(gset)
        
    fgseaRes <- fgsea(pathways = gset, stats = df_gsea_sorted, minSize  = 10, maxSize  = 500, eps = 0.0)

    # Get positively enriched pathways
    topPathwaysUp <- fgseaRes[ES > 0][order(padj), ]

    # Get negatively enriched pathways
    topPathwaysDown <- fgseaRes[ES < 0][order(padj), ]

    # to save
    file_name_up = paste0(outbase, title, '_up_regulated.csv')
    file_name_down = paste0(outbase, title, '_down_regulated.csv')
        
    fwrite(topPathwaysUp, file=file_name_up, sep=",", sep2=c("", " ", ""))
    fwrite(topPathwaysDown, file=file_name_down, sep=",", sep2=c("", " ", ""))
}
