#library(pheatmap)
#library(RColorBrewer)
#library(ggplot2)
#library(ggrepel)
#library(gridExtra)
#library(grid)
#library(lemon)
print('we are in')


if(FALSE)
{
	#BiocManager::install('DropletUtils', lib = '/tmp2/b04401068/R/R/4.0/')
	install.packages('Matrix', lib = '/tmp2/b04401068/lib/R/')
}

#Fig. 1 all breast cancer
if(FALSE)
{
    library(ggplot2)
    library(ggrepel)
    library(cowplot)
    pdf("figure/Fig_1_lower.pdf", onefile = TRUE, width = 20, height = 10)
    #all breast gene priority demonstration
    mat <- read.csv( paste('R_data/all_gene_prob_thres_drug_mean.csv'),sep=',',row.names = 1);
    mat$type <- factor(mat$assignment)
    graph <- ggplot(mat, aes(y = median_prob, x = perc)) +
        #geom_smooth(method='lm')+
        #geom_point(alpha = 0.6, size = 4) +
        geom_point( shape = 1, aes(colour = type), alpha = 0.6) +
        scale_colour_manual( values = c('CommonEssCore' = 'red', 'lower 10%' = 'blue', 'priority'='black'), name = '', breaks=c('CommonEssCore','lower 10%', 'priority'), labels = c('Common Essentials', 'Essential<10%', 'priority'))+
        #scale_color_manual( values = c('other' = 'azure2', 'Low' = 'cornflowerblue', 'Medium' = 'coral', 'High' = 'red') )+
        #scale_fill_brewer(palette = 'Dark2')+
        #scale_x_continuous(limits = c(4.5,7))+
        #scale_y_continuous(limits = c(0,0.05))+
        #guides( color=FALSE)+
        #geom_text_repel( aes(label = redox), max.overlaps = 50)+
        labs(x = 'fraction of cell line', y = 'median probability dependency score')
    p1 <- graph + theme_classic(base_size = 18)+theme(legend.position = c(0.2,0.9))
    #p1 <- graph + theme_classic(base_size = 18)

    #label further prioritization
    mat <- read.csv( paste('R_data/all_gene_prob_thres_drug_mean.csv'),sep=',',row.names = 1);
    #mat <- read.csv( paste('R_data/all_gene_PIK3CA_2.csv'),sep=',',row.names = 1);
    print(head(mat))
    mat$sel <- factor(mat$sel)
    graph <- ggplot(mat, aes(y = median_prob, x = perc)) +
        #geom_smooth(method='lm')+
        #geom_point(alpha = 0.6, size = 4) +
        #geom_point( shape = 1, aes(colour = sel), alpha = 0.6) +
        geom_point(  aes(colour = sel), alpha = 0.6) +
        scale_colour_manual( name = '', values = c('below' = 'azure2', 'above' = 'blue', 'exist_drug' = 'red'), breaks = c('below','above','exist_drug'), labels = c('lower dep','higher dep','existing drug') )+
        #scale_color_manual( values = c('other' = 'azure2', 'Low' = 'cornflowerblue', 'Medium' = 'coral', 'High' = 'red') )+
        #scale_fill_brewer(palette = 'Dark2')+
        #scale_x_continuous(limits = c(4.5,7))+
        #scale_y_continuous(limits = c(0,0.05))+
        #guides( color=FALSE)+
        geom_text_repel( aes(label = anno), max.overlaps = 50)+
        labs(x = 'fraction of cell line', y = 'median probability dependency score')
    p2 <- graph + theme_classic(base_size = 18)+theme(legend.position = c(0.2,0.9))
    #p2 <- graph + theme_classic(base_size = 18)
    graph <- plot_grid(p1,p2, labels = c('A','B'))
    print(graph)
    dev.off()
}
#combine enrichment
if(FALSE)
{
	library(clusterProfiler)
	library(org.Hs.eg.db)
	library(enrichplot)
	library(ggplot2)
	library(ggrepel)
	mat <- read.csv( paste('R_data/all_breast_gene_pri.csv'),sep=',',row.names = 1);
    	#pdf("figure/breast_cancer_GOenrich_all.pdf", onefile = TRUE, width = 12, height = 9)
    	pdf("figure/Fig_1_with_enrich.pdf", onefile = TRUE, width = 16, height = 16)
	print('Fig_1')
	diff <- mat[ mat[,'assignment'] == 'priority',]
	eg <- bitr( diff$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')
	eg <- eg[!duplicated(eg$SYMBOL),]
	rownames(diff) <- diff$gene
	rownames(eg) <- eg$SYMBOL
	diff <- diff[eg$SYMBOL,]
	diff <- cbind(diff,eg)
	gene <- diff[,'ENTREZID']
	bp <- enrichGO(gene          = gene,
			OrgDb         = org.Hs.eg.db,
			ont           = "MF",
			pAdjustMethod = "BH",
			pvalueCutoff  = 0.01,
			qvalueCutoff  = 0.05,
		readable      = TRUE)
	#head(ego)
	bp <- pairwise_termsim(bp)
	bp2 <- simplify(bp, cutoff=0.7, by="p.adjust", select_fun=min)
	#p1 <- emapplot(bp)
	p3 <- emapplot(bp2)
    mat <- read.csv( paste('R_data/all_breast_gene_pri.csv'),sep=',',row.names = 1);
    mat$type <- factor(mat$assignment)
    graph <- ggplot(mat, aes(y = -median_eff, x = perc)) +
        #geom_smooth(method='lm')+
        #geom_point(alpha = 0.6, size = 4) +
        geom_point( shape = 1, aes(colour = type), alpha = 0.6) +
        scale_colour_manual( values = c('CommonEssCore' = 'red', 'lower 10%' = 'blue', 'priority'='black'), name = '', breaks=c('CommonEssCore','lower 10%', 'priority'), labels = c('Common Essentials', 'Essential<10%', 'priority'))+
        #scale_color_manual( values = c('other' = 'azure2', 'Low' = 'cornflowerblue', 'Medium' = 'coral', 'High' = 'red') )+
        #scale_fill_brewer(palette = 'Dark2')+
        #scale_x_continuous(limits = c(4.5,7))+
        #scale_y_continuous(limits = c(0,0.05))+
        #guides( color=FALSE)+
        #geom_text_repel( aes(label = redox), max.overlaps = 50)+
        labs(x = 'fraction of cell line', y = 'median genetic effect score')
    p1 <- graph + theme_classic(base_size = 18)+theme(legend.position = c(0.2,0.9))
    #p1 <- graph + theme_classic(base_size = 18)

    #label further prioritization
    mat <- read.csv( paste('R_data/all_breast_gene_PIK3CA.csv'),sep=',',row.names = 1);
    #mat <- read.csv( paste('R_data/all_gene_PIK3CA_2.csv'),sep=',',row.names = 1);
    print(head(mat))
    mat$sel <- factor(mat$sel)
    graph <- ggplot(mat, aes(y = -median_eff, x = perc)) +
        #geom_smooth(method='lm')+
        #geom_point(alpha = 0.6, size = 4) +
        #geom_point( shape = 1, aes(colour = sel), alpha = 0.6) +
        geom_point(  aes(colour = sel), alpha = 0.6) +
        scale_colour_manual( name = '', values = c('below' = 'azure2', 'above' = 'blue', 'exist_drug' = 'red'), breaks = c('below','above','exist_drug'), labels = c('lower dep','higher dep','existing drug') )+
        #scale_color_manual( values = c('other' = 'azure2', 'Low' = 'cornflowerblue', 'Medium' = 'coral', 'High' = 'red') )+
        #scale_fill_brewer(palette = 'Dark2')+
        #scale_x_continuous(limits = c(4.5,7))+
        #scale_y_continuous(limits = c(0,0.05))+
        #guides( color=FALSE)+
        geom_text_repel( aes(label = anno), max.overlaps = 50)+
        labs(x = 'fraction of cell line', y = 'median genetic effect score' )
    p2 <- graph + theme_classic(base_size = 18)+theme(legend.position = c(0.2,0.9))
    #p2 <- graph + theme_classic(base_size = 18)
    upper <- cowplot::plot_grid(p1,p2, labels = c('A','B'), label_size = 25)
    #graph <- cowplot::plot_grid(upper, p3, labels = c('','D'), rel_heights = c(1,2), ncol = 1)
    graph <- cowplot::plot_grid(upper, p3, labels = c('','C'), label_size = 25, ncol = 1)

    print(graph)
	#graph <- cowplot::plot_grid(p1, p2, ncol=2, labels = LETTERS[1:2])
	#p2 <- p2 + labs(title = 'breast')
	#print(p2)
	dev.off()
	#saveRDS(kk,file = 'working_data/BT474_6hr_oraKEGG_down_more.rds') 
}

#Fig. 2 subtype
#Breast cancer subtype percent-median gene effect plot
if(FALSE)
{
    library(ggplot2)
    library(ggrepel)
    library(cowplot)
    print('Fig_2')
    pdf("figure/Fig_2.pdf", onefile = TRUE, width = 17, height = 17)
    mat <- read.csv( paste('R_data/all_gene_prob_lowthres.csv'),sep=',',row.names = 1);
    #mat <- read.csv( paste('R_data/all_gene.csv'),sep=',',row.names = 1);
    print(head(mat))
    mat$type <- factor(mat$assignment)
    mat$molecular <- factor(mat$molecular)
    graph <- ggplot(mat, aes(y = median_prob, x = perc)) +
	facet_wrap(~molecular, nrow = 1)+
        #geom_smooth(method='lm')+
        #geom_point(alpha = 0.6, size = 4) +
        geom_point( shape = 1, aes(colour = type), alpha = 0.6) +
        #scale_colour_manual( values = c('CommonEssCore' = 'red', 'lower 10%' = 'blue', 'priority'='black'), name = '', breaks=c('CommonEssCore','lower 10%', 'priority'), labels = c('Common Essentials', 'Essential<10%', 'priority'))+
        scale_colour_manual( values = c('CommonEssCore' = 'red', 'not priority' = 'blue', 'priority'='black'), name = '', breaks=c('CommonEssCore','not priority', 'priority'), labels = c('Common Essentials', 'not priority', 'priority'))+
        #scale_color_manual( values = c('other' = 'azure2', 'Low' = 'cornflowerblue', 'Medium' = 'coral', 'High' = 'red') )+
        #scale_fill_brewer(palette = 'Dark2')+
        #scale_x_continuous(limits = c(4.5,7))+
        #scale_y_continuous(limits = c(0,0.05))+
        #guides( color=FALSE)+
        #geom_text_repel( aes(label = redox), max.overlaps = 50)+
        labs(x = 'fraction of cell line', y = 'median dependency score')
    p1 <- graph + theme_classic(base_size = 18)+theme(legend.position = c(0.1,0.9))
    #graph <- graph + theme_classic(base_size = 18)+theme(legend.position = c(0.1,0.9))

    #mat <- read.csv( paste('R_data/breast_priority_open_target.csv'),sep=',',row.names = 1);
    mat <- read.csv( paste('R_data/all_gene_prob_thres_drug_lowthres.csv'),sep=',',row.names = 1);
    #mat <- read.csv( paste('R_data/all_gene_PIK3CA_2.csv'),sep=',',row.names = 1);
    print(head(mat))
    mat$molecular <- factor(mat$molecular)
    mat <- mat[mat$sel != 'below',]
    mat$sel <- factor(mat$sel)
    graph <- ggplot(mat, aes(y = median_prob, x = perc)) +
	facet_wrap(~molecular, nrow = 1)+
        #geom_smooth(method='lm')+
        #geom_point(alpha = 0.6, size = 4) +
        #geom_point( shape = 1, aes(colour = sel), alpha = 0.6) +
        geom_point(  aes(colour = sel), alpha = 0.6) +
        #scale_colour_manual( name = '', values = c('below' = 'azure2', 'above' = 'blue', 'exist_drug' = 'red'), breaks = c('below','above','exist_drug'), labels = c('lower priority','top priority','existing drug') )+
        scale_colour_manual( name = '', values = c('above' = 'blue', 'exist_drug' = 'red'), breaks = c('above','exist_drug'), labels = c('top priority','existing drug') )+
        #scale_color_manual( values = c('other' = 'azure2', 'Low' = 'cornflowerblue', 'Medium' = 'coral', 'High' = 'red') )+
        #scale_fill_brewer(palette = 'Dark2')+
        scale_x_continuous(limits = c(0,1))+
        scale_y_continuous(limits = c(0,1))+
        #guides( color=FALSE)+
        geom_text_repel( aes(label = anno), max.overlaps = 50)+
        labs(x = 'fraction of cell line', y = 'median dependency score')
    p2 <- graph + theme_classic(base_size = 18)+theme(legend.position = c(0.1,0.9))
    graph <- plot_grid(p1,p2, labels = c('B','C'), label_size = 30,ncol = 1)
    print(graph)
    dev.off()
}

#Fig. 3
if(TRUE)
{
	library(clusterProfiler)
	library(org.Hs.eg.db)
	library(enrichplot)
	library(ggplot2)
	#mat <- read.csv( paste('R_data/all_gene_prob.csv'),sep=',',row.names = 1);
        mat <- read.csv( paste('R_data/all_gene_prob_thres_drug_lowthres.csv'),sep=',',row.names = 1);
	#mat <- read.csv( paste('R_data/all_gene_PIK3CA.csv'),sep=',',row.names = 1);
    	#pdf("figure/Fig_3.pdf", onefile = TRUE, width = 25, height = 20)
    	pdf("figure/Fig_3_lowthres_tmp.pdf", onefile = TRUE, width = 17, height = 25)
	print('Fig 3')
	tmp <- list()
	j <- 1
	for( i in c('ER+','HER2+','TNBC'))
	{
		diff <- mat[ mat[,'molecular'] == i,]
		eg <- bitr( diff$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')
		eg <- eg[!duplicated(eg$SYMBOL),]
		rownames(diff) <- diff$gene
		rownames(eg) <- eg$SYMBOL
		diff <- diff[eg$SYMBOL,]
		diff <- cbind(diff,eg)
		gene <- diff[,'ENTREZID']
		bp <- enrichGO(gene          = gene,
				OrgDb         = org.Hs.eg.db,
				ont           = "MF",
				pAdjustMethod = "BH",
				pvalueCutoff  = 0.01,
				qvalueCutoff  = 0.05,
			readable      = TRUE)
		#head(ego)
		bp <- pairwise_termsim(bp)
		bp2 <- simplify(bp, cutoff=0.7, by="p.adjust", select_fun=min)
		#p1 <- emapplot(bp)
		p2 <- emapplot(bp2)
		#graph <- cowplot::plot_grid(p1, p2, ncol=2, labels = LETTERS[1:2])
		write.table( bp2 , file = paste0('R_data/GO_result_',i,'_new_lowthres.tsv'),quote = F,row.names = T,col.names = T,sep = '\t')
		tmp[[j]] <- p2
		j <- j + 1
		#print(p2)
	}
    	graph <- cowplot::plot_grid(tmp[[1]],tmp[[2]], labels = c('A','B'), label_size = 30, ncol = 2)
    	graph_2 <- cowplot::plot_grid(tmp[[3]],NULL, labels = c('C',''), label_size = 30, ncol = 2)
    	#graph <- cowplot::plot_grid(tmp[[1]],tmp[[2]],tmp[[3]], labels = c('A','B','C'), label_size = 30, ncol = 3)
    library(ComplexHeatmap)
    library(circlize)
    cor <- read.csv( paste0('R_data/priority_heatmap_data.csv'), header = TRUE, row.names = 1, sep = ',')
    anno_row <- read.csv( paste0('R_data/priority_heatmap_anno_row.csv'), header = TRUE, sep = ',')
    anno_col <- read.csv( paste0('R_data/priority_heatmap_anno_col.csv'), header = TRUE,row.names = 1, sep = ',')
    rownames(anno_row) <- anno_row[,'name']
    cor <- as.matrix(cor)
    anno_row[,'molecular_type'] = as.factor(anno_row[,'molecular_type'])
    #column_ha = HeatmapAnnotation( ER = anno_col[,'ER'], HER2 = anno_col[,'HER2'], TNBC = anno_col[,'TNBC'], heatmap_legend_param = list( show_legend = FALSE, col = list( ER  = c("True" = "red","False" = "white"), HER2 = c("True" = "green","False" = "white"), TNBC = c("True" = "black", "False" = "white") ) ))
    column_ha = HeatmapAnnotation( ER = anno_col[,'ER'], HER2 = anno_col[,'HER2'], TNBC = anno_col[,'TNBC'], show_legend = FALSE,  col = list( ER  = c("True" = "red","False" = "white"), HER2 = c("True" = "green","False" = "white"), TNBC = c("True" = "black", "False" = "white") ) )
    print(head(anno_row))
    row_ha = rowAnnotation( molecular = anno_row[,'molecular_type'], col = list( molecular = c("ER+" = "red","HER2+" = "green", "TNBC" = "black")))
    #row_ha = rowAnnotation( molecular = anno_row[,'molecular_type'])
    #col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
    col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
    #graph_h <- Heatmap(cor, show_row_names = TRUE, show_column_names = TRUE, col =col_fun,heatmap_legend_param = list(title = 'dep', at = c(-2,0,2), labels = c(2,0,-2)), name = 'priority_heatmap', column_title = 'HIG_subtype',top_annotation = column_ha, left_annotation = row_ha)
    graph_h <- Heatmap(cor, show_row_names = TRUE, show_column_names = TRUE, col =col_fun,heatmap_legend_param = list(title = 'dep', at = c(0,0.5,1), labels = c(0,0.5,1)), name = 'priority_heatmap', column_title = 'top priority subtype',top_annotation = column_ha, left_annotation = row_ha)
    gb <- grid.grabExpr(draw(graph_h))
    #gb <- cowplot::plot_grid(gb,NULL, labels = c('D',''), label_size = 30, ncol = 2, rel_widths = c(2,1))
    #gb <- cowplot::plot_grid(gb,NULL, labels = c('D',''), label_size = 30, ncol = 2,)
   graph <- cowplot::plot_grid(graph,graph_2,gb, labels = c('','','D'), label_size = 30, ncol = 1, rel_heights = c(0.8,0.8,1))
	print(graph)
	dev.off()
	#saveRDS(kk,file = 'working_data/BT474_6hr_oraKEGG_down_more.rds') 
}

#Fig. 4
if( FALSE )
{
	print('Fig 4')
    library(ggplot2)
    library(ggrepel)
    library(cowplot)
    pdf("figure/Fig_4.pdf", onefile = TRUE, width = 17, height = 25)
   
    mat <- read.csv( paste('R_data/mutation_association_lowthres.txt'),sep='\t');
    print(head(mat))
    mat$sel <- factor(mat$assign)
    graph <- ggplot(mat, aes(y = fdr, x = cohen_d)) +
	facet_wrap(~molecular, nrow = 1)+
        #geom_smooth(method='lm')+
        #geom_point(alpha = 0.6, size = 4) +
        #geom_point( shape = 1, aes(colour = sel), alpha = 0.6) +
        geom_point(  aes(colour = sel), alpha = 0.6) +
        scale_colour_manual( guide = FALSE, values = c('not_sig' = 'azure2', 'sig' = 'red') )+
        #scale_color_manual( values = c('other' = 'azure2', 'Low' = 'cornflowerblue', 'Medium' = 'coral', 'High' = 'red') )+
        #scale_fill_brewer(palette = 'Dark2')+
        scale_x_continuous(limits = c(-5,5))+
        #scale_y_continuous(limits = c(0,0.05))+
        guides( colour= "none")+
        geom_text_repel( aes(label = anno), max.overlaps = 50)+
        labs(x = 'Cohen_d', y = '-log_10(FDR)', title = 'MUT' )
    #graph <- graph + theme_classic(base_size = 18)+theme(legend.position = c(0.1,0.9))
    p2 <- graph + theme_classic(base_size = 18)
    
    mat <- read.csv( paste('R_data/cnv_association_derive_lowthres.txt'),sep='\t');
    print(head(mat))
    mat$sel <- factor(mat$assign)
    graph <- ggplot(mat, aes(y = fdr, x = cohen_d)) +
	facet_wrap(~molecular, nrow = 1)+
        #geom_smooth(method='lm')+
        #geom_point(alpha = 0.6, size = 4) +
        #geom_point( shape = 1, aes(colour = sel), alpha = 0.6) +
        geom_point(  aes(colour = sel), alpha = 0.6) +
        scale_colour_manual( guide = FALSE, values = c('not_sig' = 'azure2', 'sig' = 'red') )+
        #scale_color_manual( values = c('other' = 'azure2', 'Low' = 'cornflowerblue', 'Medium' = 'coral', 'High' = 'red') )+
        #scale_fill_brewer(palette = 'Dark2')+
        #scale_x_continuous(limits = c(4.5,7))+
        #scale_y_continuous(limits = c(0,0.05))+
        scale_x_continuous(limits = c(-5,5))+
        guides( colour= "none")+
        #geom_text_repel( aes(label = name),box.padding = 0.5, max.overlaps = 50)+
        geom_text_repel( aes(label = anno), max.overlaps = Inf, force = 10)+
        labs(x = 'Cohen_d', y = '-log_10(FDR)', title = 'CNV' )
    #graph <- graph + theme_classic(base_size = 18)+theme(legend.position = c(0.1,0.9))
    p3 <- graph + theme_classic(base_size = 18)
    #p1 <- plot_grid(p1,NULL, labels = c('A',''), label_size = 30,ncol = 2, rel_widths = c(2,1))
    graph <- plot_grid(p2,p3, labels = c('A','B'), label_size = 30,ncol = 1)
    print(graph)
    dev.off()
}



#Fig. 4
if( FALSE )
{
	print('Fig 4')
    library(ggplot2)
    library(ggrepel)
    library(cowplot)
    pdf("figure/Fig_4.pdf", onefile = TRUE, width = 17, height = 25)
    mat <- read.csv( paste('R_data/mut_association_breast.csv'),sep=',');
    print(head(mat))
    mat[ mat[,'g2_status'] == 'cnv', 'g2_status'] <- 'CNV'
    mat[ mat[,'g2_status'] == 'mut', 'g2_status'] <- 'MUT'
    mat$sel <- factor(mat$sel)
    graph <- ggplot(mat, aes(y = fdr, x = -cohen_d)) +
	facet_wrap(~g2_status, nrow = 1)+
        #geom_smooth(method='lm')+
        #geom_point(alpha = 0.6, size = 4) +
        #geom_point( shape = 1, aes(colour = sel), alpha = 0.6) +
        geom_point(  aes(colour = sel), alpha = 0.6) +
        scale_colour_manual( guide = FALSE, values = c('not_sig' = 'azure2', 'sig' = 'red') )+
        #scale_color_manual( values = c('other' = 'azure2', 'Low' = 'cornflowerblue', 'Medium' = 'coral', 'High' = 'red') )+
        #scale_fill_brewer(palette = 'Dark2')+
        #scale_x_continuous(limits = c(4.5,7))+
        #scale_y_continuous(limits = c(0,0.05))+
        guides( colour= "none")+
        geom_text_repel( aes(label = name))+
        labs(x = 'cohen_d', y = '-log_10(FDR)' )
    #p1 <- graph + theme_classic(base_size = 18)+theme(legend.position = c(0.1,0.9))
    p1 <- graph + theme_classic(base_size = 18)
    
    mat <- read.csv( paste('R_data/mut_association_breast_type_specific.csv'),sep=',');
    print(head(mat))
    mat <- mat[ mat[,'g2_status'] == 'mut',]
    mat$sel <- factor(mat$sel)
    graph <- ggplot(mat, aes(y = fdr, x = -cohen_d)) +
	facet_wrap(~molecular, nrow = 1)+
        #geom_smooth(method='lm')+
        #geom_point(alpha = 0.6, size = 4) +
        #geom_point( shape = 1, aes(colour = sel), alpha = 0.6) +
        geom_point(  aes(colour = sel), alpha = 0.6) +
        scale_colour_manual( guide = FALSE, values = c('not_sig' = 'azure2', 'sig' = 'red') )+
        #scale_color_manual( values = c('other' = 'azure2', 'Low' = 'cornflowerblue', 'Medium' = 'coral', 'High' = 'red') )+
        #scale_fill_brewer(palette = 'Dark2')+
        #scale_x_continuous(limits = c(4.5,7))+
        #scale_y_continuous(limits = c(0,0.05))+
        guides( colour= "none")+
        geom_text_repel( aes(label = name), max.overlaps = 50)+
        labs(x = 'cohen_d', y = '-log_10(FDR)', title = 'MUT' )
    #graph <- graph + theme_classic(base_size = 18)+theme(legend.position = c(0.1,0.9))
    p2 <- graph + theme_classic(base_size = 18)
    
    mat <- read.csv( paste('R_data/mut_association_breast_type_specific.csv'),sep=',');
    print(head(mat))
    mat <- mat[ mat[,'g2_status'] == 'cnv',]
    mat$sel <- factor(mat$sel)
    graph <- ggplot(mat, aes(y = fdr, x = -cohen_d)) +
	facet_wrap(~molecular, nrow = 1)+
        #geom_smooth(method='lm')+
        #geom_point(alpha = 0.6, size = 4) +
        #geom_point( shape = 1, aes(colour = sel), alpha = 0.6) +
        geom_point(  aes(colour = sel), alpha = 0.6) +
        scale_colour_manual( guide = FALSE, values = c('not_sig' = 'azure2', 'sig' = 'red') )+
        #scale_color_manual( values = c('other' = 'azure2', 'Low' = 'cornflowerblue', 'Medium' = 'coral', 'High' = 'red') )+
        #scale_fill_brewer(palette = 'Dark2')+
        #scale_x_continuous(limits = c(4.5,7))+
        #scale_y_continuous(limits = c(0,0.05))+
        guides( colour= "none")+
        #geom_text_repel( aes(label = name),box.padding = 0.5, max.overlaps = 50)+
        geom_text_repel( aes(label = name), max.overlaps = 50)+
        labs(x = 'cohen_d', y = '-log_10(FDR)', title = 'CNV' )
    #graph <- graph + theme_classic(base_size = 18)+theme(legend.position = c(0.1,0.9))
    p3 <- graph + theme_classic(base_size = 18)
    p1 <- plot_grid(p1,NULL, labels = c('A',''), label_size = 30,ncol = 2, rel_widths = c(2,1))
    graph <- plot_grid(p1,p2,p3, labels = c('','B','C'), label_size = 30,ncol = 1)
    print(graph)
    dev.off()
}

#Breast cancer subtype percent-median gene effect plot
if(FALSE)
{
    library(ggplot2)
    library(ggrepel)
    pdf("figure/breast_cancer_gene_Fig_1.pdf", onefile = TRUE, width = 9, height = 9)
    mat <- read.csv( paste('R_data/er_gene.csv'),sep=',',row.names = 1);
    print(head(mat))
    mat$type <- factor(mat$assignment)
    graph <- ggplot(mat, aes(y = -median_eff, x = perc)) +
        #geom_smooth(method='lm')+
        #geom_point(alpha = 0.6, size = 4) +
        geom_point( shape = 1, aes(colour = type), alpha = 0.6) +
        scale_colour_manual( values = c('CommonEssCore' = 'red', 'lower 10%' = 'blue', 'priority'='black') )+
        #scale_color_manual( values = c('other' = 'azure2', 'Low' = 'cornflowerblue', 'Medium' = 'coral', 'High' = 'red') )+
        #scale_fill_brewer(palette = 'Dark2')+
        #scale_x_continuous(limits = c(4.5,7))+
        #scale_y_continuous(limits = c(0,0.05))+
        #guides( color=FALSE)+
        #geom_text_repel( aes(label = redox), max.overlaps = 50)+
        labs(title = 'ER positive', x = 'percentage of cell line dep > 0.5', y = 'median gene effect score' )
    graph <- graph + theme_classic(base_size = 18)
            theme(legend.position = c(0.9,0.9))
    print(graph)
    mat <- read.csv( paste('R_data/her_gene.csv'),sep=',',row.names = 1);
    print(head(mat))
    mat$type <- factor(mat$assignment)
    graph <- ggplot(mat, aes(y = -median_eff, x = perc)) +
        #geom_smooth(method='lm')+
        #geom_point(alpha = 0.6, size = 4) +
        geom_point( shape = 1, aes(colour = type), alpha = 0.6) +
        scale_colour_manual( values = c('CommonEssCore' = 'red', 'lower 10%' = 'blue', 'priority'='black') )+
        #scale_color_manual( values = c('other' = 'azure2', 'Low' = 'cornflowerblue', 'Medium' = 'coral', 'High' = 'red') )+
        #scale_fill_brewer(palette = 'Dark2')+
        #scale_x_continuous(limits = c(4.5,7))+
        #scale_y_continuous(limits = c(0,0.05))+
        #guides( color=FALSE)+
        #geom_text_repel( aes(label = redox), max.overlaps = 50)+
        labs(title = 'HER2 positive', x = 'percentage of cell line dep > 0.5', y = 'median gene effect score' )
    graph <- graph + theme_classic(base_size = 18)
            theme(legend.position = c(0.9,0.9))
    print(graph)
    mat <- read.csv( paste('R_data/tnbc_gene.csv'),sep=',',row.names = 1);
    print(head(mat))
    mat$type <- factor(mat$assignment)
    graph <- ggplot(mat, aes(y = -median_eff, x = perc)) +
        #geom_smooth(method='lm')+
        #geom_point(alpha = 0.6, size = 4) +
        geom_point( shape = 1, aes(colour = type), alpha = 0.6) +
        scale_colour_manual( values = c('CommonEssCore' = 'red', 'lower 10%' = 'blue', 'priority'='black') )+
        #scale_color_manual( values = c('other' = 'azure2', 'Low' = 'cornflowerblue', 'Medium' = 'coral', 'High' = 'red') )+
        #scale_fill_brewer(palette = 'Dark2')+
        #scale_x_continuous(limits = c(4.5,7))+
        #scale_y_continuous(limits = c(0,0.05))+
        #guides( color=FALSE)+
        #geom_text_repel( aes(label = redox), max.overlaps = 50)+
        labs(title = 'TNBC', x = 'percentage of cell line dep > 0.5', y = 'median gene effect score' )
    graph <- graph + theme_classic(base_size = 18)
            theme(legend.position = c(0.9,0.9))
    print(graph)
    dev.off()
}



#Breast cancer subtype percent-median gene effect plot
if(FALSE)
{
    library(ggplot2)
    library(ggrepel)
    pdf("figure/breast_cancer_gene_Fig_1_facet.pdf", onefile = TRUE, width = 17, height = 9)
    mat <- read.csv( paste('R_data/all_gene.csv'),sep=',',row.names = 1);
    print(head(mat))
    mat$type <- factor(mat$assignment)
    mat$molecular <- factor(mat$molecular)
    graph <- ggplot(mat, aes(y = -median_eff, x = perc)) +
	facet_wrap(~molecular, nrow = 1)+
        #geom_smooth(method='lm')+
        #geom_point(alpha = 0.6, size = 4) +
        geom_point( shape = 1, aes(colour = type), alpha = 0.6) +
        scale_colour_manual( values = c('CommonEssCore' = 'red', 'lower 10%' = 'blue', 'priority'='black') )+
        #scale_color_manual( values = c('other' = 'azure2', 'Low' = 'cornflowerblue', 'Medium' = 'coral', 'High' = 'red') )+
        #scale_fill_brewer(palette = 'Dark2')+
        #scale_x_continuous(limits = c(4.5,7))+
        #scale_y_continuous(limits = c(0,0.05))+
        #guides( color=FALSE)+
        #geom_text_repel( aes(label = redox), max.overlaps = 50)+
        labs(x = 'fraction of cell line dep > 0.5', y = 'median gene effect score' )
    graph <- graph + theme_classic(base_size = 18)+theme(legend.position = c(0.1,0.9))
    print(graph)
    dev.off()
}


#Breast cancer subtype percent-median gene effect plot_ median cut-off
if(FALSE)
{
    library(ggplot2)
    library(ggrepel)
    pdf("figure/breast_cancer_gene_median_cut_Fig_2_facet.pdf", onefile = TRUE, width = 17, height = 9)
    mat <- read.csv( paste('R_data/all_gene_PIK3CA.csv'),sep=',',row.names = 1);
    print(head(mat))
    mat$molecular <- factor(mat$molecular)
    mat$sel <- factor(mat$sel)
    graph <- ggplot(mat, aes(y = -median_eff, x = perc)) +
	facet_wrap(~molecular, nrow = 1)+
        #geom_smooth(method='lm')+
        #geom_point(alpha = 0.6, size = 4) +
        geom_point( shape = 1, aes(colour = sel), alpha = 0.6) +
        scale_colour_manual( values = c('below' = 'red', 'above' = 'blue') )+
        #scale_color_manual( values = c('other' = 'azure2', 'Low' = 'cornflowerblue', 'Medium' = 'coral', 'High' = 'red') )+
        #scale_fill_brewer(palette = 'Dark2')+
        #scale_x_continuous(limits = c(4.5,7))+
        #scale_y_continuous(limits = c(0,0.05))+
        #guides( color=FALSE)+
        geom_text_repel( aes(label = anno), max.overlaps = 50)+
        labs(x = 'fraction of cell line dep > 0.5', y = 'median gene effect score' )
    graph <- graph + theme_classic(base_size = 18)+theme(legend.position = c(0.1,0.9))
    print(graph)
    dev.off()
}

#Breast cancer subtype percent-median gene effect plot
if(FALSE)
{
    library(ggplot2)
    library(ggrepel)
    pdf("figure/breast_cancer_gene_Fig_1_all.pdf", onefile = TRUE, width = 10, height = 9)
    mat <- read.csv( paste('R_data/all_breast_gene_pri.csv'),sep=',',row.names = 1);
    print(head(mat))
    mat$type <- factor(mat$assignment)
    graph <- ggplot(mat, aes(y = -median_eff, x = perc)) +
        #geom_smooth(method='lm')+
        #geom_point(alpha = 0.6, size = 4) +
        geom_point( shape = 1, aes(colour = type), alpha = 0.6) +
        scale_colour_manual( values = c('CommonEssCore' = 'red', 'lower 10%' = 'blue', 'priority'='black') )+
        #scale_color_manual( values = c('other' = 'azure2', 'Low' = 'cornflowerblue', 'Medium' = 'coral', 'High' = 'red') )+
        #scale_fill_brewer(palette = 'Dark2')+
        #scale_x_continuous(limits = c(4.5,7))+
        #scale_y_continuous(limits = c(0,0.05))+
        #guides( color=FALSE)+
        #geom_text_repel( aes(label = redox), max.overlaps = 50)+
        labs(x = 'fraction of cell line dep > 0.5', y = 'median gene effect score', title = 'breast' )
    graph <- graph + theme_classic(base_size = 18)+theme(legend.position = c(0.2,0.9))
    print(graph)
    dev.off()
}

#Breast cancer subtype percent-median gene effect plot_ median cut-off
if(FALSE)
{
    library(ggplot2)
    library(ggrepel)
    pdf("figure/breast_cancer_gene_median_cut_Fig_2_facet_2.pdf", onefile = TRUE, width = 17, height = 9)
    #mat <- read.csv( paste('R_data/breast_priority_open_target.csv'),sep=',',row.names = 1);
    mat <- read.csv( paste('R_data/all_gene_PIK3CA_2.csv'),sep=',',row.names = 1);
    print(head(mat))
    mat$molecular <- factor(mat$molecular)
    mat$sel <- factor(mat$sel)
    graph <- ggplot(mat, aes(y = -median_eff, x = perc)) +
	facet_wrap(~molecular, nrow = 1)+
        #geom_smooth(method='lm')+
        #geom_point(alpha = 0.6, size = 4) +
        #geom_point( shape = 1, aes(colour = sel), alpha = 0.6) +
        geom_point(  aes(colour = sel), alpha = 0.6) +
        scale_colour_manual( values = c('below' = 'azure2', 'above' = 'blue', 'exist_drug' = 'red') )+
        #scale_color_manual( values = c('other' = 'azure2', 'Low' = 'cornflowerblue', 'Medium' = 'coral', 'High' = 'red') )+
        #scale_fill_brewer(palette = 'Dark2')+
        #scale_x_continuous(limits = c(4.5,7))+
        #scale_y_continuous(limits = c(0,0.05))+
        #guides( color=FALSE)+
        geom_text_repel( aes(label = anno), max.overlaps = 50)+
        labs(x = 'fraction of cell line dep > 0.5', y = 'median gene effect score' )
    graph <- graph + theme_classic(base_size = 18)+theme(legend.position = c(0.1,0.9))
    print(graph)
    dev.off()
}

#Breast cancer percent-median gene effect plot_ median cut-off
if(FALSE)
{
    library(ggplot2)
    library(ggrepel)
    pdf("figure/breast_cancer_gene_median_cut_all_Fig_2_facet_2.pdf", onefile = TRUE, width = 10, height = 9)
    #mat <- read.csv( paste('R_data/all_breast_gene_PIK3CA_2.csv'),sep=',',row.names = 1);
    mat <- read.csv( paste('R_data/all_gene_PIK3CA_2.csv'),sep=',',row.names = 1);
    print(head(mat))
    mat$sel <- factor(mat$sel)
    graph <- ggplot(mat, aes(y = -median_eff, x = perc)) +
        #geom_smooth(method='lm')+
        #geom_point(alpha = 0.6, size = 4) +
        #geom_point( shape = 1, aes(colour = sel), alpha = 0.6) +
        geom_point(  aes(colour = sel), alpha = 0.6) +
        scale_colour_manual( values = c('below' = 'azure2', 'above' = 'blue', 'exist_drug' = 'red') )+
        #scale_color_manual( values = c('other' = 'azure2', 'Low' = 'cornflowerblue', 'Medium' = 'coral', 'High' = 'red') )+
        #scale_fill_brewer(palette = 'Dark2')+
        #scale_x_continuous(limits = c(4.5,7))+
        #scale_y_continuous(limits = c(0,0.05))+
        #guides( color=FALSE)+
        geom_text_repel( aes(label = anno), max.overlaps = 50)+
        labs(x = 'fraction of cell line dep > 0.5', y = 'median gene effect score', title = 'breast' )
    graph <- graph + theme_classic(base_size = 18)+theme(legend.position = c(0.2,0.9))
    print(graph)
    dev.off()
}



#Breast cancer mutation association
if(FALSE)
{
    library(ggplot2)
    library(ggrepel)
    pdf("figure/breast_cancer_gene_mutation_all_cell_Fig_3_facet.pdf", onefile = TRUE, width = 17, height = 9)
    mat <- read.csv( paste('R_data/mut_association_breast.csv'),sep=',');
    print(head(mat))
    mat$sel <- factor(mat$sel)
    graph <- ggplot(mat, aes(y = fdr, x = -cohen_d)) +
	facet_wrap(~g2_status, nrow = 1)+
        #geom_smooth(method='lm')+
        #geom_point(alpha = 0.6, size = 4) +
        #geom_point( shape = 1, aes(colour = sel), alpha = 0.6) +
        geom_point(  aes(colour = sel), alpha = 0.6) +
        scale_colour_manual( values = c('not_sig' = 'azure2', 'sig' = 'red') )+
        #scale_color_manual( values = c('other' = 'azure2', 'Low' = 'cornflowerblue', 'Medium' = 'coral', 'High' = 'red') )+
        #scale_fill_brewer(palette = 'Dark2')+
        #scale_x_continuous(limits = c(4.5,7))+
        #scale_y_continuous(limits = c(0,0.05))+
        guides( colour= "none")+
        geom_text_repel( aes(label = name), max.overlaps = 50)+
        labs(x = 'cohen_d', y = '-log_10(FDR)' )
    graph <- graph + theme_classic(base_size = 18)+theme(legend.position = c(0.1,0.9))
    print(graph)
    dev.off()
}

#Breast cancer mutation association
if(FALSE)
{
    library(ggplot2)
    library(ggrepel)
    pdf("figure/breast_cancer_gene_mutation_type_speicifc_Fig_3_facet.pdf", onefile = TRUE, width = 17, height = 9)
    mat <- read.csv( paste('R_data/mut_association_breast_type_specific.csv'),sep=',');
    print(head(mat))
    mat <- mat[ mat[,'g2_status'] == 'mut',]
    mat$sel <- factor(mat$sel)
    graph <- ggplot(mat, aes(y = fdr, x = -cohen_d)) +
	facet_wrap(~molecular, nrow = 1)+
        #geom_smooth(method='lm')+
        #geom_point(alpha = 0.6, size = 4) +
        #geom_point( shape = 1, aes(colour = sel), alpha = 0.6) +
        geom_point(  aes(colour = sel), alpha = 0.6) +
        scale_colour_manual( values = c('not_sig' = 'azure2', 'sig' = 'red') )+
        #scale_color_manual( values = c('other' = 'azure2', 'Low' = 'cornflowerblue', 'Medium' = 'coral', 'High' = 'red') )+
        #scale_fill_brewer(palette = 'Dark2')+
        #scale_x_continuous(limits = c(4.5,7))+
        #scale_y_continuous(limits = c(0,0.05))+
        guides( colour= "none")+
        geom_text_repel( aes(label = name), max.overlaps = 50)+
        labs(x = 'cohen_d', y = '-log_10(FDR)', title = 'mut' )
    graph <- graph + theme_classic(base_size = 18)+theme(legend.position = c(0.1,0.9))
    print(graph)
    mat <- read.csv( paste('R_data/mut_association_breast_type_specific.csv'),sep=',');
    print(head(mat))
    mat <- mat[ mat[,'g2_status'] == 'cnv',]
    mat$sel <- factor(mat$sel)
    graph <- ggplot(mat, aes(y = fdr, x = -cohen_d)) +
	facet_wrap(~molecular, nrow = 1)+
        #geom_smooth(method='lm')+
        #geom_point(alpha = 0.6, size = 4) +
        #geom_point( shape = 1, aes(colour = sel), alpha = 0.6) +
        geom_point(  aes(colour = sel), alpha = 0.6) +
        scale_colour_manual( values = c('not_sig' = 'azure2', 'sig' = 'red') )+
        #scale_color_manual( values = c('other' = 'azure2', 'Low' = 'cornflowerblue', 'Medium' = 'coral', 'High' = 'red') )+
        #scale_fill_brewer(palette = 'Dark2')+
        #scale_x_continuous(limits = c(4.5,7))+
        #scale_y_continuous(limits = c(0,0.05))+
        guides( colour= "none")+
        #geom_text_repel( aes(label = name),box.padding = 0.5, max.overlaps = 50)+
        geom_text_repel( aes(label = name), max.overlaps = 50)+
        labs(x = 'cohen_d', y = '-log_10(FDR)', title = 'cnv' )
    graph <- graph + theme_classic(base_size = 18)+theme(legend.position = c(0.1,0.9))
    print(graph)
    dev.off()
}




if(FALSE)
{
	library(clusterProfiler)
	library(org.Hs.eg.db)
	library(enrichplot)
	mat <- read.csv( paste('R_data/all_gene_PIK3CA.csv'),sep=',',row.names = 1);
	for( i in c('ER+','HER2+','TNBC'))
	{
		diff <- mat[ mat[,'molecular'] == i,]
		eg <- bitr( diff$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')
		eg <- eg[!duplicated(eg$SYMBOL),]
		rownames(diff) <- diff$gene
		rownames(eg) <- eg$SYMBOL
		diff <- diff[eg$SYMBOL,]
		diff <- cbind(diff,eg)
		gene <- diff[,'ENTREZID']
		bp <- enrichGO(gene          = gene,
				OrgDb         = org.Hs.eg.db,
				ont           = "MF",
				pAdjustMethod = "BH",
				pvalueCutoff  = 0.01,
				qvalueCutoff  = 0.05,
			readable      = TRUE)
		bp <- pairwise_termsim(bp)
		bp2 <- simplify(bp, cutoff=0.7, by="p.adjust", select_fun=min)
		bpread <- setReadable(bp2, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
		write.table( bpread , file = paste0('R_data/GO_result_',i,'.tsv'),quote = F,row.names = T,col.names = T,sep = '\t')
	}
}


if(FALSE)
{
	library(clusterProfiler)
	library(org.Hs.eg.db)
	library(enrichplot)
	library(ggplot2)
	mat <- read.csv( paste('R_data/all_breast_gene_pri.csv'),sep=',',row.names = 1);
	diff <- mat[ mat[,'assignment'] == 'priority',]
	eg <- bitr( diff$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')
	eg <- eg[!duplicated(eg$SYMBOL),]
	rownames(diff) <- diff$gene
	rownames(eg) <- eg$SYMBOL
	diff <- diff[eg$SYMBOL,]
	diff <- cbind(diff,eg)
	gene <- diff[,'ENTREZID']
	bp <- enrichGO(gene          = gene,
			OrgDb         = org.Hs.eg.db,
			ont           = "MF",
			pAdjustMethod = "BH",
			pvalueCutoff  = 0.01,
			qvalueCutoff  = 0.05,
		readable      = TRUE)
	#head(ego)
	bp <- pairwise_termsim(bp)
	bp2 <- simplify(bp, cutoff=0.7, by="p.adjust", select_fun=min)
	bpread <- setReadable(bp2, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
	write.table( bpread , file = paste0('R_data/GO_result_breast.tsv'),quote = F,row.names = T,col.names = T,sep = '\t')
	#p1 <- emapplot(bp)
	#saveRDS(kk,file = 'working_data/BT474_6hr_oraKEGG_down_more.rds') 
}

if(FALSE)
{
	library(clusterProfiler)
	library(org.Hs.eg.db)
	library(enrichplot)
	library(ggplot2)
	mat <- read.csv( paste('R_data/all_breast_gene_pri.csv'),sep=',',row.names = 1);
    	pdf("figure/breast_cancer_GOenrich_all.pdf", onefile = TRUE, width = 12, height = 9)
	diff <- mat[ mat[,'assignment'] == 'priority',]
	eg <- bitr( diff$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')
	eg <- eg[!duplicated(eg$SYMBOL),]
	rownames(diff) <- diff$gene
	rownames(eg) <- eg$SYMBOL
	diff <- diff[eg$SYMBOL,]
	diff <- cbind(diff,eg)
	gene <- diff[,'ENTREZID']
	bp <- enrichGO(gene          = gene,
			OrgDb         = org.Hs.eg.db,
			ont           = "MF",
			pAdjustMethod = "BH",
			pvalueCutoff  = 0.01,
			qvalueCutoff  = 0.05,
		readable      = TRUE)
	#head(ego)
	bp <- pairwise_termsim(bp)
	bp2 <- simplify(bp, cutoff=0.7, by="p.adjust", select_fun=min)
	#p1 <- emapplot(bp)
	p2 <- emapplot(bp2)
	#graph <- cowplot::plot_grid(p1, p2, ncol=2, labels = LETTERS[1:2])
	p2 <- p2 + labs(title = 'breast')
	print(p2)
	dev.off()
	#saveRDS(kk,file = 'working_data/BT474_6hr_oraKEGG_down_more.rds') 
}


if(FALSE)
{
	library(clusterProfiler)
	library(org.Hs.eg.db)
	library(enrichplot)
	library(ggplot2)
	mat <- read.csv( paste('R_data/all_gene_PIK3CA.csv'),sep=',',row.names = 1);
    	pdf("figure/breast_cancer_GOenrich.pdf", onefile = TRUE, width = 12, height = 9)
	for( i in c('ER+','HER2+','TNBC'))
	{
		diff <- mat[ mat[,'molecular'] == i,]
		eg <- bitr( diff$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')
		eg <- eg[!duplicated(eg$SYMBOL),]
		rownames(diff) <- diff$gene
		rownames(eg) <- eg$SYMBOL
		diff <- diff[eg$SYMBOL,]
		diff <- cbind(diff,eg)
		gene <- diff[,'ENTREZID']
		bp <- enrichGO(gene          = gene,
				OrgDb         = org.Hs.eg.db,
				ont           = "MF",
				pAdjustMethod = "BH",
				pvalueCutoff  = 0.01,
				qvalueCutoff  = 0.05,
			readable      = TRUE)
		#head(ego)
		bp <- pairwise_termsim(bp)
		bp2 <- simplify(bp, cutoff=0.7, by="p.adjust", select_fun=min)
		#p1 <- emapplot(bp)
		p2 <- emapplot(bp2)
		#graph <- cowplot::plot_grid(p1, p2, ncol=2, labels = LETTERS[1:2])
		p2 <- p2 + labs(title = i)
		print(p2)
	}
	dev.off()
	#saveRDS(kk,file = 'working_data/BT474_6hr_oraKEGG_down_more.rds') 
}


if(FALSE)
{
    library(ComplexHeatmap)
    library(circlize)
    cor <- read.csv( paste0('R_data/priority_heatmap_data.csv'), header = TRUE, row.names = 1, sep = ',')
    anno_row <- read.csv( paste0('R_data/priority_heatmap_anno_row.csv'), header = TRUE, sep = ',')
    anno_col <- read.csv( paste0('R_data/priority_heatmap_anno_col.csv'), header = TRUE,row.names = 1, sep = ',')
    rownames(anno_row) <- anno_row[,'name']
    cor <- as.matrix(cor)
    anno_row[,'molecular_type'] = as.factor(anno_row[,'molecular_type'])
    #column_ha = HeatmapAnnotation( ER = anno_col[,'ER'], HER2 = anno_col[,'HER2'], TNBC = anno_col[,'TNBC'], heatmap_legend_param = list( show_legend = FALSE, col = list( ER  = c("True" = "red","False" = "white"), HER2 = c("True" = "green","False" = "white"), TNBC = c("True" = "black", "False" = "white") ) ))
    column_ha = HeatmapAnnotation( ER = anno_col[,'ER'], HER2 = anno_col[,'HER2'], TNBC = anno_col[,'TNBC'], show_legend = FALSE,  col = list( ER  = c("True" = "red","False" = "white"), HER2 = c("True" = "green","False" = "white"), TNBC = c("True" = "black", "False" = "white") ) )
    print(head(anno_row))
    row_ha = rowAnnotation( molecular = anno_row[,'molecular_type'], col = list( molecular = c("ER+" = "red","HER2+" = "green", "TNBC" = "black")))
    #row_ha = rowAnnotation( molecular = anno_row[,'molecular_type'])
    col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
    pdf("figure/priority_heatmap.pdf", onefile = TRUE, width = 13, height = 10)
    graph <- Heatmap(cor, show_row_names = TRUE, show_column_names = TRUE, col =col_fun,heatmap_legend_param = list(title = 'dep', at = c(-2,0,2), labels = c(2,0,-2)), name = 'priority_heatmap', column_title = 'priority_subtype',top_annotation = column_ha, left_annotation = row_ha)
    print(graph)
    dev.off()
}

if(FALSE)
{
	#library(clusterProfiler)
	library(pathview)
	gene <- readRDS(file = 'working_data/glist_MDAMB468_24hr.rds')
	#print(kk2$ID)
	for( i in c('hsa04710','hsa03030','hsa04110','hsa03040','hsa03050','hsa03410'))
	{
		print(i)
		pathview(gene.data  = gene,
                     pathway.id = i,
                     species    = "hsa",
                     limit      = list(gene=max(abs(gene)), cpd=1),
		     #kegg.dir = './enrich_figure/BT474_6hr_ora/')
		     kegg.dir = './enrich_figure/MDAMB468_24hr_gsea/')
	}
	#pdf( paste0("figure/pathview_test.pdf"), onefile = TRUE, width = 10, height = 7)
	#pdf( paste0("figure/pathview_test.pdf"), onefile = TRUE)
	#print(hsa03010)
	#dev.off()
}


if(FALSE)
{
	library(enrichplot)
	library(clusterProfiler)
	kk2 <-readRDS(file = 'working_data/BT474_24hr_gseaKEGG_all_stat.rds') 
	geneList <- readRDS(file = 'working_data/glist_BT474_24hr.rds')
	set.seed(123)

	pdf( paste0("enrich_figure/enrichplot_BT474_24hr_gsea_cnetplot.pdf"), onefile = TRUE, width = 10, height = 10)
	edo <- pairwise_termsim(kk2)
	#print(edo)
	p1 <- emapplot(edo, showCategory = 41, cex_label_category = 0.8)
	print(p1)
	print('success')
	dev.off()
}

if(FALSE)
{
	library(Seurat)
	sc.data <- Read10X(data.dir = './working_data/aggr/')
	sc <- CreateSeuratObject( counts = sc.data, project = 'sc')
	print('Finish creating')
	saveRDS(sc, file = './working_data/sc_create.rds')
	sc[['RNA']]@counts
	sc[['percent.mt']] <- PercentageFeatureSet(sc, pattern = '^MT-')
	#cell_num <- table(strsplit2(rownames(sc@meta.data),split = '-')[,2])
	#attributes(cell_num)$dimnames[[1]] <- c('M4340','M4341','M4342','M4343')
	#sample_name <- ulist(apply(as.data.frame(cell_num), 1, function(x){rep(x[1], x[2])}))
	#sc <- AddMetaData( object = sc, metadata = sample_name, col.name = 'sample')
	sc <- subset(sc, subset = percent.mt < 10 & nCount_RNA > 5000)
	sc <- NormalizeData(sc, normalization.method = 'LogNormalize', scale.factor = 10000)
	sc[['RNA']]@data
	saveRDS(sc, file = './working_data/sc_normal.rds')
	print('Finish normalization')
	#all.genes <- rownames(sc)
	#sc <- ScaleData(sc, features = all.genes)
	#saveRDS(sc, file = './working_data/sc_scale.rds')
	#print('Finish scaling')
}

if(FALSE)
{
	#library(RcisTarget)
	#data(motifAnnotations_hgnc)
	#print('load annotations')
	#tmp <- motifAnnotations_hgnc[(motif == 'cisbp__M2125'),]
	#print(tmp)
	#tmp <- motifAnnotations_hgnc[(motif == 'transfac_pro__M07935'),]
	#print(tmp)
	library(Seurat)
	sc <- readRDS(file = 'working_data/sc_create.rds')
	sc[['RNA']]@counts
	sc[['percent.mt']] <- PercentageFeatureSet(sc, pattern = '^MT-')
	sc <- subset(sc, subset = percent.mt < 10 & nCount_RNA > 5000)

	name = colnames(x = sc)
	name = grep("-2|-3|-4",x = name)
	print(head(name))
	sc <- subset(sc, cells = name)

	#sc <- NormalizeData(sc, normalization.method = 'RC', scale.factor = 1000000)
	#motif <- GetAssayData(object = sc, assay = 'RNA', slot = 'data')
	motif <- as.data.frame(sc@assays$RNA@counts)
	rownames(motif) <- rownames(x = sc)

	#library(feather)
	print(dim(motif))
	print(head(motif))
	#print(motif$geneSet)
	saveRDS(motif, file = 'working_data/sc_count_1-3d.rds')
	#write_feather(motif, 'working_data/sc_count.feather')
	#print(head(motif[(geneSet == '100'),]))

}

if(FALSE)
{
	library(DESeq2)
	#library(feather)
	#count <- read_feather('working_data/sc_count.feather')
	count <- readRDS(file = 'working_data/sc_count.rds')
	print('reading')
	sf <- estimateSizeFactorsForMatrix(count)
	print('size factor')
	ncount <- t(t(count)/sf)
	print('count')

	saveRDS(ncount, file = 'working_data/sc_count_dseq.rds')
	print('write')
	print(head(rownames(ncount)))
}


if(FALSE)
{
	#library(DESeq2)
	#library(feather)
	library(matrixStats)
	count <- readRDS(file = 'working_data/sc_count.rds')
	#count <- read_feather('working_data/sc_countdseq.feather')
	print(head(rownames(count)))
	print(head(colnames(count)))
	means <- rowMeans(count)
	count <- as.matrix(count)
	vars <- rowVars(count)
	cv <- vars / means / means

	pdf("deseq.pdf", onefile = TRUE)
	plot( NULL, xaxt = 'n', yaxt = 'n', log = 'xy',
	     xlim = c( 1e-1, 3e5), ylim = c( 0.005, 8),
	     xlab = 'average normalized read count',
	     ylab = 'cv2')
        axis( 1, 10^(-1:5), c('0.1', '1', '10', '100', '1000', 
			   expression(10^4), expression(10^5)))
        axis( 2, 10^(-2:1), c('0.01', '0.1', '1', '10'), las = 2)
        abline( h=10^(-2:1), v = 10^(-1:5), col = '#D0D0D0', lwd=2)
        points( means, cv, pch = 20, cex = 0.2) 
	dev.off()
}


if(FALSE)
{
	library(Seurat)
	library(DropletUtils)
	sc <- readRDS(file = './working_data/sc_normal.rds')
	#print(sc)
	#print(head(rownames(sc@meta.data)))
	#print(rownames(sc@meta.data)[9000])
	#print( head(scrna))
	scrna <- GetAssayData(object = sc, assay = 'RNA', slot = 'counts')
	write10xCounts(x = scrna, path = './working_data/all_counts/')
	print('counts finish')
	scrna <- GetAssayData(object = sc, assay = 'RNA', slot = 'data')
	write10xCounts(x = scrna, path = './working_data/all_normalize/')
	print('normalize finish')
	#scrna <- GetAssayData(object = sc, assay = 'RNA', slot = 'scale.data')
	#saveRDS(scrna, file = './working_data/sc_scale_data.rds')
	#write10xCounts(x = scrna, path = './working_data/all_scale/')
	#print('scale finish')
	#print( head(scrna))
}

if(FALSE)
{
	#library(DropletUtils)
	#library(Matrix)
	sc <- readRDS(file = './working_data/sc_scale_data.rds')
	#write10xCounts(x = sc, path = './working_data/all_scale/')
	print(typeof(sc))
	print(dim(sc))
	write.table( sc,file = 'sc_scale_data.csv',quote = F,row.names = T,col.names = T,sep = '\t')
	#print(head(colnames(sc)))
	#print(head(rownames(sc)))
}

if(FALSE)
{
	library(Seurat)
	sc <- readRDS(file = './working_data/sc_normal.rds')
	#print(colnames(x=sc))
	#sc <- subset(sc, cells = colnames(x = sc)[8550:dim(sc)[2]])
	#name = colnames(x=sc)
	#name = grep("-2|-3|-4",x = name)
	#print(head(name))
	#sc <- subset(sc, cells = name)
	#print(typeof(sc))
	cell <- read.csv( paste0('working_data/cell.csv'), header = TRUE, row.names = 1, sep = ',')
	cell <- rownames(cell)
	sc <- subset(sc, cells = cell)
	print(dim(sc))
	#motif <- as.data.frame(sc@assays$RNA@data)
	#print(head(rownames(motif)))
	#rownames(motif) <- rownames(x = sc)
	#print(head(rownames(motif)))
	#saveRDS(motif, file = './working_data/sc_normal_matrix.rds')
	sc <- FindVariableFeatures(sc, selection.method = 'vst', nfeatures = 500)
	print('find variable')
	print(sc@assays$var.genes)
	write.table( sc@assays$RNA@var.features,file = 'working_data/sc_vargene_500.csv',quote = F,row.names = T,col.names = T,sep = '\t')
	sc <- FindVariableFeatures(sc, selection.method = 'vst', nfeatures = 1000)
	print('find variable')
	print(sc@assays$var.genes)
	write.table( sc@assays$RNA@var.features,file = 'working_data/sc_vargene_1000.csv',quote = F,row.names = T,col.names = T,sep = '\t')
	sc <- FindVariableFeatures(sc, selection.method = 'vst', nfeatures = 1500)
	print('find variable')
	print(sc@assays$var.genes)
	write.table( sc@assays$RNA@var.features,file = 'working_data/sc_vargene_1500.csv',quote = F,row.names = T,col.names = T,sep = '\t')
	sc <- FindVariableFeatures(sc, selection.method = 'vst', nfeatures = 2000)
	print('find variable')
	print(sc@assays$var.genes)
	write.table( sc@assays$RNA@var.features,file = 'working_data/sc_vargene_2000.csv',quote = F,row.names = T,col.names = T,sep = '\t')
	#sc <- ScaleData(sc, features = rownames(sc))
	#print('scaling')
	#s.genes <- cc.genes$s.genes
	#g2m.genes <- cc.genes$g2m.genes
	#sc <- CellCycleScoring(sc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
	#print('cell cycle')
	#head(sc[[]])
	#saveRDS(sc, file = './working_data/sc_scale_2.rds')
	#write.table( sc[[]],file = 'sc_scale_data_3.csv',quote = F,row.names = T,col.names = T,sep = '\t')
}

if(FALSE)
{
	library(Seurat)
	library(DropletUtils)
	sc <- readRDS(file = './working_data/sc_scale_2.rds')
	print(typeof(sc))
	print(dim(sc))
	sc <- ScaleData(sc, vars.to.regress = c('S.Score','G2M.Score'), features = rownames(sc))
	print('scale data')
	saveRDS(sc, file = './working_data/sc_scale_regress.rds')
	print('save regress')
	scrna <- GetAssayData(object = sc, assay = 'RNA', slot = 'scale.data')
	write.table( scrna ,file = 'working_data/sc_scale_data_regress.csv',quote = F,row.names = T,col.names = T,sep = '\t')
	#scrna <- GetAssayData(object = sc, assay = 'RNA', slot = 'scale.data')
	#write10xCounts(x = scrna, path = './working_data/all_scale/')
}

if(FALSE)
{
	library(Seurat)
	sc <- readRDS(file = './working_data/sc_scale_regress.rds')
	#write10xCounts(x = scrna, path = './working_data/all_scale/')
}
#gsh metabolism
if( FALSE )
{
	library(corrplot)
	pdf("1-1_ferro_correlation.pdf", onefile = TRUE)
	dep <- read.csv( 'working_data/ferro_gene_cor_1.csv', row.names = 1);
	dep <- as.matrix(dep)
	corrplot(dep, type = 'lower', diag = FALSE, order = 'hclust', addrect = 3,  tl.cex = 0.5, method = 'circle', tl.col = 'black', col = colorRampPalette(rev(brewer.pal(n=10, name = 'RdBu')))(200)) 
	dep <- read.csv( 'working_data/ferro_gene_cor_2.csv', row.names = 1);
	dep <- as.matrix(dep)
	corrplot(dep, type = 'lower', diag = FALSE, order = 'hclust', addrect = 3,  tl.cex = 0.5, method = 'circle', tl.col = 'black', col = colorRampPalette(rev(brewer.pal(n=10, name = 'RdBu')))(200)) 
	dep <- read.csv( 'working_data/ferro_gene_cor_3.csv', row.names = 1);
	dep <- as.matrix(dep)
	corrplot(dep, type = 'lower', diag = FALSE, order = 'hclust', addrect = 3,  tl.cex = 0.5, method = 'circle', tl.col = 'black', col = colorRampPalette(rev(brewer.pal(n=10, name = 'RdBu')))(200)) 
	dep <- read.csv( 'working_data/ferro_gene_cor_4.csv', row.names = 1);
	dep <- as.matrix(dep)
	corrplot(dep, type = 'lower', diag = FALSE, order = 'hclust', addrect = 3,  tl.cex = 0.5, method = 'circle', tl.col = 'black', col = colorRampPalette(rev(brewer.pal(n=10, name = 'RdBu')))(200)) 
	#corrplot.mixed(dep, order = 'AOE') 
	#corrplot(dep, order = 'hclust', addrect = 3, method = 'circle', tl.col = 'black',  tl.cex = 0.5, col = colorRampPalette(rev(brewer.pal(n=10, name = 'RdBu')))(200)) 
	#corrplot(dep, is.corr = FALSE, method = 'circle', na.label = 'NA', cl.pos = 'n', tl.cex = 3, tl.col = 'black', col = colorRampPalette(rev(brewer.pal(n=10, name = 'RdBu')))(200) ) 
	dev.off()
}

#geini3
if(FALSE)
{
	library(GENIE3)
	tf <- read.csv( '3grn_6/tf_3.csv', header = FALSE)
	target <- read.csv( '3grn_6/target_3.csv', header = FALSE)
	tf <- as.vector(tf[,1])
	target <- as.vector(target[,1])
	print(head(target))
	set.seed(123)
	exp <- read.csv( '3grn_6/exp_3.csv', row.names = 1)
	exp <- as.matrix(exp)
	#print(head(exp))
	print('finish reading')
	weightMat <- GENIE3(exp, regulators = tf, targets = target, verbose = TRUE, nCores = 20, nTrees = 50)
	write.table( weightMat,file = 'genie3_3_6_weight.csv',quote = F,row.names = T,col.names = T,sep = ',')
}

#Rcistarget
if(FALSE)
{
	genelist <- read.table('data/TF_target_YAP.csv',sep = '\t',header = TRUE,stringsAsFactors = FALSE)[,1]
	geneLists <- list(yap = genelist)
	files <- list.files(path = 'working_data/1-3d_std', pattern = '*', full.names = TRUE, recursive = FALSE)
	print(head(files))
	#print(substring(head(files),26))
	print(substring(head(files),23))
	if(FALSE)
	{
	lapply(files, function(x){
		       geneLists[[substring(x,18)]] <- read.table(x, sep = '\t', header = TRUE, stringsAsFactors = FALSE)[,1]
})
	}
	for (x in files)
	{
	       geneLists[[substring(x,23)]] <- read.table(x, sep = '\t', header = TRUE, stringsAsFactors = FALSE)[,1]
	}
	print(head(geneLists))
	library(RcisTarget)
	data(motifAnnotations_hgnc)
	print('load annotations')
	motifRankings <- importRankings('working_data/hg38_500_100.feather')
	print('load rankings')
	motifEnrichmentTable_wGenes <- cisTarget( geneLists, motifRankings,
						 motifAnnot = motifAnnotations_hgnc)
	print('motif enrich')
	saveRDS(motifEnrichmentTable_wGenes, file = 'working_data/1-3d_std_motif.rds')
}
if(FALSE)
{
	#library(RcisTarget)
	#data(motifAnnotations_hgnc)
	#print('load annotations')
	#tmp <- motifAnnotations_hgnc[(motif == 'cisbp__M2125'),]
	#print(tmp)
	#tmp <- motifAnnotations_hgnc[(motif == 'transfac_pro__M07935'),]
	#print(tmp)
	motif <- readRDS(file = 'working_data/1-3d_std_motif.rds')
	library(feather)
	print(dim(motif))
	#print(motif$geneSet)
	write_feather(motif, 'working_data/1-3d_std_motif.feather')
	#print(head(motif[(geneSet == '100'),]))

}
if(FALSE)
{
	genelist <- read.table('working_data/ferro_gene_add.txt',sep = '\t',header = TRUE,stringsAsFactors = FALSE)[,1]
	print(head(genelist))
	library(RcisTarget)
	geneLists <- list(geneListName = genelist)
	geneLists[['another']] <- c('SOX4','APP','TARDBP','ARNTL','NRF2L2','HIF1A','PTEN','YAP1','GJA1','AQP11','TEAD1','TEAD4','EPAS1','ARID5B','JUN')
	data(motifAnnotations_hgnc)
	print('load annotations')
	motifRankings <- importRankings('working_data/hg38_500_100.feather')
	print('load rankings')
	motifEnrichmentTable_wGenes <- cisTarget( geneLists, motifRankings,
						 motifAnnot = motifAnnotations_hgnc)
	print('motif enrich')
	saveRDS(motifEnrichmentTable_wGenes, file = 'working_data/test_motif.rds')
}
if(FALSE)
{
	library(RcisTarget)
	data(motifAnnotations_hgnc)
	print('load annotations')
	print(motifAnnotations_hgnc[(TF %in% c('YAP1')),])
	#tmp <- motifAnnotations_hgnc[(motif == 'cisbp__M2125'),]
	#print(tmp)
	#tmp <- motifAnnotations_hgnc[(motif == 'transfac_pro__M07935'),]
	#print(tmp)
	#motif <- readRDS(file = 'working_data/test_motif.rds')
	#print(dim(motif))
	#print(head(motif[geneSet == 'another']))

}

#removing doublets
if( FALSE )
{
	exp_tmp <- readRDS(file = './working_data/sc_count_1-3d.rds')
	exp <- as.matrix(exp_tmp)
	rm(exp_tmp)
	cell <- read.csv( paste0('working_data/cell.csv'), header = TRUE, row.names = 1, sep = ',')
	cell <- rownames(cell)
	exp <- exp[,cell]
	saveRDS(exp, file = './working_data/sc_count_1-3d.rds')
}

#scenic initiation
if( FALSE )
{
	exp <- readRDS(file = './working_data/sc_count_1-3d.rds')
	exp <- as.matrix(exp)
	print(head(colnames(exp)))
	print(head(rownames(exp)))
	print(dim(exp))
	print('finish loading')

	#dir.create('./int')
	library(SCENIC)
	org <- 'hgnc'
	dbDir <- '/tmp2/b04401068/working_data/hg19'
	scenicOptions <- initializeScenic(org = org, dbDir = dbDir, nCores = 20)
	saveRDS(scenicOptions, file = 'int/scenicOptions.Rds')

	genesKept <- geneFiltering( exp, scenicOptions = scenicOptions)
	interestingGenes <- c('YAP1','NFE2L2','ACSL4','TFRC')
	print(interestingGenes[which(!interestingGenes %in% genesKept)])
	exp_filter <- exp[genesKept,]
	print(dim(exp_filter))
	##exp_filter for corrrelation
	exp_filter <- log2(exp_filter+1)
	rm(exp)
}

#exporting geneskept
if(FALSE)
{
	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.rds')
	#tmp <- getIntName(scenicOptions, 'genie3wm')
	#print(tmp)
	genesKept <- loadInt( scenicOptions, 'genesKept')
	print(head(genesKept))
	write.table( genesKept,file = 'working_data/genesKept.csv',quote = F,row.names = T,col.names = T,sep = ',')
}

#scenic runGenie3
if( FALSE )
{
	exp_tmp <- readRDS(file = './working_data/sc_count_1-3d.rds')
	exp <- as.matrix(exp_tmp)
	rm(exp_tmp)
	cell <- read.csv( paste0('working_data/cell.csv'), header = TRUE, row.names = 1, sep = ',')
	cell <- rownames(cell)
	exp <- exp[,cell]
	print(head(colnames(exp)))
	print(head(rownames(exp)))
	print(dim(exp))
	print('finish loading')

	#dir.create('./int')
	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	#tmp <- getIntName(scenicOptions, 'genie3wm')
	#print(tmp)
	genesKept <- loadInt( scenicOptions, 'genesKept')
	exp_filter_tmp <- exp[genesKept,]
	exp_filter <- log2(exp_filter_tmp+1)
	rm(exp)
	rm(exp_filter_tmp)
	runGenie3(exp_filter, scenicOptions, nParts = 10)
}

#runGenie3 check parrellel on 10
if( FALSE )
{
	dir.create('./working_data/para/')
	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	regulator <- getDbTfs(scenicOptions)
	genekept <- readRDS(file = './int/1.1_genesKept.Rds')
	exp <- readRDS(file = './working_data/sc_count_1-3d.rds')
	exp <- exp[genekept,]
	exp <- log2(exp+1)
	regulator <- regulator[regulator %in% rownames(exp)]
	print(length(regulator))
	exp <- t(exp)
	exp <- as.data.frame(exp)
	print(dim(exp))
	write.table( regulator, file = paste0('working_data/para/regulator.csv'), quote = F,row.names = T,col.names = T,sep = '\t')
	library(feather)
	increment <- 950
	for( i in c(1:10) )
	{
		start <- (i-1)*increment + 1
		end <- start + increment - 1
		if( end > length(genekept) )
			end <- length(genekept)
		target <- genekept[c(start:end)]
		print(i)
		print(length(target))
		output <- union( regulator, target)
		write.table( target, file = paste0('working_data/para/target_', toString(i),'.csv'), quote = F,row.names = T,col.names = T,sep = '\t')
		write_feather( exp[,output], paste0('working_data/para/exp_', toString(i),'.feather'))
	}
}

#runGenie3 retrieve result (rerun option, failed)
if( FALSE )
{
	exp_tmp <- readRDS(file = './working_data/sc_count_1-3d.rds')
	exp <- as.matrix(exp_tmp)
	rm(exp_tmp)
	print(dim(exp))
	print('finish loading')

	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	genesKept <- loadInt( scenicOptions, 'genesKept')
	exp_filter_tmp <- exp[genesKept,]
	print(dim(exp_filter_tmp))
	exp_filter <- log2(exp_filter_tmp+1)
	rm(exp)
	rm(exp_filter_tmp)
	runGenie3(exp_filter, scenicOptions, nParts = 10, resumePreviousRun = TRUE)
}

#copying code for runGenie3
if(FALSE)
{
	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	print(getSettings(scenicOptions, 'modules/weightThreshold'))
	weight <- readRDS(file = 'int/1.3_GENIE3_weightMatrix_part_1.Rds')
	print(head(colnames(weight)))
	print(head(rownames(weight)))
	print(dim(weight))
	weight <- as.matrix(weight)
	linklist <- GENIE3::getLinkList(weight, threshold=getSettings(scenicOptions, 'modules/weightThreshold'))
	print('finish link list')
	rm(weight)
	colnames(linklist) <- c('TF','Target','weight')
	linklist <- linklist[order(linklist[,'weight'],decreasing=TRUE),]
	linklist <- unique(linklist)
	rownames(linklist) <- NULL
	saveRDS(linklist, file=getIntName(scenicOptions,'genie3ll'))
}

#run coexnetwork
if( FALSE )
{
	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	scenicOptions@settings$verbose <- TRUE
	scenicOptions@settings$seed <- 123
	print('starting')
	start_time <- Sys.time()
	scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
	end_time <- Sys.time()
	print('finishing')
	print(start_time - end_time)
	saveRDS(scenicOptions, file = 'int/scenicOptions.Rds')
}

#run regulon
if( FALSE )
{
	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	scenicOptions@settings$verbose <- TRUE
	scenicOptions@settings$seed <- 123
	scenicOptions@settings$nCores <- 20
	print('starting')
	start_time <- Sys.time()
	scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
	end_time <- Sys.time()
	print('finishing')
	print(start_time - end_time)
	saveRDS(scenicOptions, file = 'int/scenicOptions.Rds')
}

#run cell
if( FALSE )
{
	exp_tmp <- readRDS(file = './working_data/sc_count_1-3d.rds')
	exp <- as.matrix(exp_tmp)
	rm(exp_tmp)
	print(dim(exp))
	print('finish loading')

	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	genesKept <- loadInt( scenicOptions, 'genesKept')
	exp_filter_tmp <- exp[genesKept,]
	exp_filter <- log2(exp_filter_tmp+1)
	rm(exp)
	rm(exp_filter_tmp)

	scenicOptions@settings$verbose <- TRUE
	scenicOptions@settings$seed <- 123
	scenicOptions@settings$nCores <- 5
	print('starting')
	start_time <- Sys.time()
	scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exp_filter)
	end_time <- Sys.time()
	print('finishing')
	print(start_time - end_time)
	saveRDS(scenicOptions, file = 'int/scenicOptions.Rds')
}

#plot heatmap
.openDev <- function(fileName, devType, ...)
{
  if(devType=="pdf")
    pdf(paste0(fileName, ".pdf"), ...)
  
  if(devType=="png")
    png(paste0(fileName, ".png", type="cairo"), ...)
  
  if(devType=="cairo_pfd") # similar to Cairo::CairoPDF?
    grDevices::cairo_pdf(paste0(fileName, ".pdf"), ...)
}

.openDevHeatmap <- function(fileName, devType)
{
  if(devType!="pdf") 
  {
    if(devType=="png") .openDev(fileName=fileName, devType=devType, width=1200,height=1200)
    if(devType!="png") .openDev(fileName=fileName, devType=devType)
    fileName <- NA
  }else{
    fileName <- paste0(fileName,".pdf")
  }
  return(fileName)
}

.closeDevHeatmap <- function(devType)
{
  if(devType!="pdf") 
  {
    dev.off()
  }
}

if(FALSE)
{
	library(NMF)
	library(SCENIC)
	library(AUCell)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	skipHeatmap <- FALSE
	regulonAUC <- readRDS(file = getIntName(scenicOptions, 'aucell_regulonAUC'))
msg <- paste0(format(Sys.time(), "%H:%M"), "\tFinished running AUCell.")
  if(getSettings(scenicOptions, "verbose")) message(msg)
  
  if(!skipHeatmap){
    msg <- paste0(format(Sys.time(), "%H:%M"), "\tPlotting heatmap...")
    if(getSettings(scenicOptions, "verbose")) message(msg)
    
    nCellsHeatmap <- min(500, ncol(regulonAUC))
    cells2plot <- sample(colnames(regulonAUC), nCellsHeatmap)
    
    cellInfo <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "cellInfo"), ifNotExists="null")   #TODO check if exists, if not... create/ignore?
    if(!is.null(cellInfo)) cellInfo <- data.frame(cellInfo)[cells2plot,,drop=F]
    colVars <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "colVars"), ifNotExists="null")
    
    fileName <- getOutName(scenicOptions, "s3_AUCheatmap")
    
    fileName <- .openDevHeatmap(fileName=fileName, devType=getSettings(scenicOptions, "devType"))
    NMF::aheatmap(getAUC(regulonAUC)[,cells2plot],
                  annCol=cellInfo,
                  annColor=colVars,
                  main="AUC",
                  sub=paste("Subset of",nCellsHeatmap," random cells"),
                  filename=fileName)
    .closeDevHeatmap(devType=getSettings(scenicOptions, "devType"))
  }
}

#plot Tsne
if( FALSE)
{
	library(Rtsne)
	library(R2HTML)
	library(SCENIC)
	library(AUCell)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	skipTsne <- FALSE
	regulonAUC <- readRDS(file = getIntName(scenicOptions, 'aucell_regulonAUC'))

	exp_tmp <- readRDS(file = './working_data/sc_count_1-3d.rds')
	exp <- as.matrix(exp_tmp)
	rm(exp_tmp)
	print(dim(exp))
	print('finish loading')

	genesKept <- loadInt( scenicOptions, 'genesKept')
	exp_filter_tmp <- exp[genesKept,]
	exprMat <- log2(exp_filter_tmp+1)
	rm(exp)
	rm(exp_filter_tmp)
if(!skipTsne){
    msg <- paste0(format(Sys.time(), "%H:%M"), "\tPlotting t-SNEs...")
    if(getSettings(scenicOptions, "verbose")) message(msg)
    
    tSNE_fileName <- tsneAUC(scenicOptions, aucType="AUC", onlyHighConf=FALSE) # default: nPcs, perpl, seed, tsne prefix
    tSNE <- readRDS(tSNE_fileName)
    
    # AUCell (activity) plots with the default tsne, as html: 
    fileName <- getOutName(scenicOptions, "s3_AUCtSNE_colAct")
    plotTsne_AUCellHtml(scenicOptions, exprMat, fileName, tSNE) #open the resulting html locally

    # Plot cell properties:
    sub <- ""; if("type" %in% names(tSNE)) sub <- paste0("t-SNE on ", tSNE$type)
    cellInfo <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "cellInfo"), ifNotExists="null") 
    colVars <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "colVars"), ifNotExists="null")
    pdf(paste0(getOutName(scenicOptions, "s3_AUCtSNE_colProps"),".pdf"))
    plotTsne_cellProps(tSNE$Y, cellInfo=cellInfo, colVars=colVars, cex=1, sub=sub)
    dev.off()
  }
}

#run binary
if( FALSE )
{
	exp_tmp <- readRDS(file = './working_data/sc_count_1-3d.rds')
	exp <- as.matrix(exp_tmp)
	rm(exp_tmp)
	print(dim(exp))
	print('finish loading')

	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	genesKept <- loadInt( scenicOptions, 'genesKept')
	exp_filter_tmp <- exp[genesKept,]
	exp_filter <- log2(exp_filter_tmp+1)
	rm(exp)
	rm(exp_filter_tmp)
	
	scenicOptions@settings$verbose <- TRUE
	scenicOptions@settings$seed <- 123
	scenicOptions@settings$nCores <- 5

	aucellApp <- plotTsne_AUCellApp(scenicOptions, exp_filter)
	savedSelections <- shiny::runApp(aucellApp, port = 4694)

	newthres <- savedSelections$thresholds
	scenicOptions@fileNames$int['aucell_thresholds',1] <- 'int/newThresholds.Rds'
	#saveRDS(newthres, file=getIntName(scenicOptions, 'aucell_thresholds'))
	#saveRDS(scenicOptions, file = 'int/scenicOptions.Rds')
}

#run binary
if( FALSE )
{
	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	thres = readRDS(file=getIntName(scenicOptions, 'aucell_thresholds'))
	print(thres)
}


#check regulon
if( FALSE )
{
	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	regulons <- loadInt(scenicOptions, 'regulons')
	print(regulons[c('ATF4')])
}

#output for python
if( FALSE )
{
	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	#regulon <- loadInt(scenicOptions, 'aucell_regulonAUC')
	regulon <- readRDS('int/4.1_binaryRegulonActivity.Rds')
	print(dim(regulon))
	print(head(rownames(regulon)))
	print(head(colnames(regulon)))
	write.table( rownames(regulon), file = paste0('working_data/regulon_auc_bin_row.csv'), quote = F,row.names = F,col.names = F,sep = '\t')
	print('write name')
	library(feather)
	#regulon <- as.data.frame(regulon@assays@data@listData$AUC)
	regulon <- as.data.frame(regulon)
	print(dim(regulon))
	write_feather( regulon, paste0('working_data/regulon_auc_bin.feather'))
}
#plot heatmap try
if( FALSE )
{
	library(SCENIC)
	set.seed(123)
	x <- rbind(rnorm(10),rnorm(10),rnorm(10))
	ridx <- c(1,1,1,1,1,10,10,10,10,10)
	print(x)
	NMF::aheatmap(x, scale = 'none', revC = TRUE, Rowv = T, Colv = ridx, color = c('white','black'), filename = 'test_order.pdf')
}

#plot heatmap
if( FALSE )
{
	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	regulon <- readRDS('int/4.1_binaryRegulonActivity.Rds')
	#rhcr <- hclust(dist(regulon))
	#chrc <- hclust(dist(t(regulon)))
	#print('finish clustering')

	anno <- read.table('working_data/rna_time.csv', sep = ',', check.names = FALSE, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
	anno <- anno[,colnames(regulon)]
	anno <- t(anno)
	anno <- as.data.frame(anno)
	#regulon <- t(t(regulon) * anno[,'day'])
	cindex <- 4-anno$day
	anno$day <- factor(anno$day)
	#anno[,'latent'] <- as.numeric(anno[,'latent'])
	print(dim(anno))
	#regulon[regulon == 0] <- -1
	#regulon <- t(t(regulon) * anno[,'latent'])
	#regulon[regulon < 0] <- -1
	#print('finish transformation')
	
	
	#ann_colors <- list(day = c('darkgreen','navyblue','oldlace'), latent = colorRampPalette(c('navy','white','firebrick3'))(50))
	ann_colors <- list(day = c('navy','white','firebrick3'), latent = colorRampPalette(c('navy','white','firebrick3'))(50))
	print('color annotation')
	
	#change threshold for some regulons
	if( TRUE )
	{
		temp <- loadInt(scenicOptions, 'aucell_regulonAUC')
		temp <- as.data.frame(temp@assays@data@listData$AUC)
		temp <- temp[,colnames(regulon)]
		regulon['TEAD1 (106g)',] <- ifelse(temp['TEAD1 (106g)',]>0.29, 1, 0)
		#regulon['TEAD1_extended (107g)',] <- ifelse(temp['TEAD1_extended (107g)',]>0.3, 1, 0)
		#regulon['ZEB1_extended (401g)',] <- ifelse(temp['ZEB1_extended (401g)',]>0.1, 1, 0)
		#regulon['HIF1A_extended (166g)',] <- ifelse(temp['HIF1A_extended (166g)',]>0.21, 1, 0)
		print('adj thres')
	}
	
	#regulon[regulon==0] <- NA
	#regulon <- regulon[rhcr$order, chrc$order]
	#anno <- anno[chrc$order,]
	print(dim(regulon))

	#pdf('cluster_output/R_regulon_bin_all_time_wide_color_0.29_adj.pdf', width = 15, height = 7, onefile = TRUE)
	#NMF::aheatmap(regulon, scale = 'none', revC = TRUE, Rowv = T, Colv = T, color = c('white','black'), filename = 'cluster_output/R_regulon_bin_all_adj.pdf')
	#NMF::aheatmap(regulon, scale = 'none', annCol = anno, revC = TRUE, Colv = T, Rowv = T, color = c('white','black'), filename = 'cluster_output/R_regulon_bin_all_time.pdf')
	NMF::aheatmap(regulon, scale = 'none', annCol = anno, annColors = ann_colors, revC = T, Colv = cindex, Rowv = T, color = c('white','black'), filename = 'R_closer.png')
	#library(NMF)
	#NMF::aheatmap(regulon, scale = 'none', annCol = anno, annColors = ann_colors, revC = TRUE, Colv = T, Rowv = T, color = colorRampPalette(c('white','navy','grey80','firebrick3'))(4))
	#NMF::aheatmap(regulon, scale = 'none', annCol = anno, annColors = ann_colors, revC = TRUE, Colv = NA, Rowv = NA, color = colorRampPalette(c('white','navy','grey80','firebrick3'))(4))
	#breaks <- seq( -1,1,length.out = 101)
	#NMF::aheatmap(regulon, scale = 'none', breaks = breaks, annCol = anno, annColors = ann_colors, revC = TRUE, Colv = T, Rowv = T, color = colorRampPalette(c('white','yellow','navy','grey80','firebrick3'))(200))
	#NMF::aheatmap(regulon, scale = 'none', breaks = breaks, annCol = anno, annColors = ann_colors, Colv = NA, Rowv = NA, color = colorRampPalette(c('white','navy','grey80','firebrick3'))(200))
	#NMF::aheatmap(regulon, scale = 'none', annCol = anno, annColors = ann_colors, Colv = NA, Rowv = NA, color = colorRampPalette(c('white','navy','grey80','firebrick3'))(200))
	#dev.off()
}


#extracting cell
if( FALSE )
{
	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	regulon <- readRDS('int/4.1_binaryRegulonActivity.Rds')
	#rhcr <- hclust(dist(regulon))
	#change threshold for some regulons
	if( TRUE )
	{
		temp <- loadInt(scenicOptions, 'aucell_regulonAUC')
		temp <- as.data.frame(temp@assays@data@listData$AUC)
		temp <- temp[,colnames(regulon)]
		regulon['TEAD1 (106g)',] <- ifelse(temp['TEAD1 (106g)',]>0.29, 1, 0)
		print('adj thres')
	}

	set.seed(123)
	chrc <- hclust(dist(t(regulon)))
	print('finish clustering')

	treecuts <- cutree(chrc, k = 3)
	saveRDS(treecuts,file = 'working_data/cellcluster.rds')
	print(treecuts)
	treecuts<-data.frame(treecuts)
	write.table( treecuts, file = paste0('working_data/cellcluster.csv'), quote = F,row.names = F,col.names = F,sep = '\t')

	if(FALSE)
	{
		anno <- read.table('working_data/rna_time.csv', sep = ',', check.names = FALSE, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
		anno <- anno[,colnames(regulon)]
		anno <- t(anno)
		anno <- as.data.frame(anno)
		cindex <- 4-anno$day
		anno$day <- factor(anno$day)
		print(dim(anno))
		#print('finish transformation')
		
		
		ann_colors <- list(day = c('navy','white','firebrick3'), latent = colorRampPalette(c('navy','white','firebrick3'))(50))
		print('color annotation')
		
			
		#regulon <- regulon[rhcr$order, chrc$order]
		#anno <- anno[chrc$order,]
		print(dim(regulon))

		#NMF::aheatmap(regulon, scale = 'none', annCol = anno, annColors = ann_colors, revC = T, Colv = cindex, Rowv = T, color = c('white','black'), filename = 'R_closer.png')
		#dev.off()
	}
}

if( FALSE )
{
	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	set.seed(123456)
	regulon <- readRDS('int/4.1_binaryRegulonActivity.Rds')

	if( TRUE )
	{
		temp <- loadInt(scenicOptions, 'aucell_regulonAUC')
		temp <- as.data.frame(temp@assays@data@listData$AUC)
		temp <- temp[,colnames(regulon)]
		regulon['TEAD1 (106g)',] <- ifelse(temp['TEAD1 (106g)',]>0.29, 1, 0)
		print('adj thres')
	}

	chrc <- hclust(dist(t(regulon)))
	print('finish clustering')

	treecuts <- cutree(chrc, k = 4)
	treecuts<-data.frame(treecuts, row.names = colnames(regulon))
	print(head(treecuts))
	write.table( treecuts, file = paste0('working_data/cellcluster.csv'), quote = F,sep = '\t')

}

#####
#output correlation matrix
if( FALSE )
{
	exp <- readRDS(file = './working_data/sc_count_1-3d.rds')
	print(dim(exp))
	print('finish loading')

	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	genesKept <- loadInt( scenicOptions, 'genesKept')
	exp <- exp[genesKept,]
	exp <- t(exp)
	exp <- as.data.frame(exp)
	print(dim(exp))
	library(feather)
	write_feather( exp, paste0('working_data/exp_1-3d_python.feather'))
}

if( FALSE )
{
	#library(arrow)
	library(SCENIC)
	cor <- read.csv( paste0('working_data/exp_spear_cor_tmp.csv'), header = TRUE, row.names = 1, sep = ',')
	#cor <- arrow::read_feather('working_data/exp_spear_cor_tmp.feather')
	#rownames(cor) <- colnames(cor)
	print(head(rownames(cor)))
	print(dim(cor))
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	saveRDS( cor, file = getIntName(scenicOptions, "corrMat"))
}

if(FALSE)
{
	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	cor <- readRDS(file = getIntName(scenicOptions, "corrMat"))
	old <- c('IQCJ.SCHIP1','NKX2.5','HLA.F','HLA.A','HLA.E','HLA.C','HLA.B','NKX3.1','BDNF.AS')
	real <- c('IQCJ-SCHIP1','NKX2-5','HLA-F','HLA-A','HLA-E','HLA-C','HLA-B','NKX3-1','BDNF-AS')
	for( i in c(1:9))
	{
		names(cor)[names(cor) == old[i]] <- real[i]
	}
	cor <- as.matrix(cor)
	print(dim(cor))
	print(head(rownames(cor)))
	saveRDS( cor, file = getIntName(scenicOptions, "corrMat"))
}

####SCENIC checking
#runGenie3 check parrellel on 10
if( FALSE )
{
	library(feather)
	exp <- readRDS(file = './int/1.3_GENIE3_weightMatrix_part_103.Rds')
	print(length(colnames(exp)))
	print(length(rownames(exp)))
	i <- 2
	regulator <- read.csv( paste0('working_data/para/regulator.csv'), header = TRUE, sep = '\t')
	target <- read.csv( paste0('working_data/para/target_', toString(i),'.csv'), header = TRUE, sep = '\t')
	exp <- read_feather( paste0('working_data/para/exp_', toString(i),'.feather'))
	exp <- t(exp)
	print(head(regulator))
	print(head(target))
	print(dim(exp))
	print(head(rownames(exp)))
	print(head(rownames(exp[target[,'x'],])))
}


#retrive parallel
if( FALSE )
{
	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	for( i in c(1,2,3,4,5,6,7,8,9,10))
	{
		mat <- read.csv( paste0('working_data/genie3_1-3d_weightmatrix_',toString(i),'.csv'), header = TRUE, row.names = 1, sep = ',')
		print(head(rownames(mat)))
		print(dim(mat))
		saveRDS( mat, file = paste0("1.3_GENIE3_weightMatrix_part_",toString(i),'.Rds'))
	}
}

#rerun 
if( FALSE )
{
	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	genesKept <- loadInt( scenicOptions, 'genesKept')
	final <- list()
	if(TRUE)
	{
		final <- readRDS( file = paste0("int/1.3_GENIE3_weightMatrix_part_",toString(1),'.Rds'))
		genesDone <- colnames(final)
		if( FALSE)
		{
		for( i in c(2,3,4,5,6,7,8,9,10))
		{
			mat <- readRDS( file = paste0("int/1.3_GENIE3_weightMatrix_part_",toString(i),'.Rds'))
			final <- cbind(final, mat)
			#print(intersect(genesDone,colnames(mat)))
			#genesDone <- union(colnames(mat), genesDone)
		}
		}
	}
	print(length(genesDone))
	print(length(genesKept))
	genes <- setdiff( genesKept, genesDone)
	#genes <- intersect( genesKept, colnames(final))
	#final <- final[,genes]
	#print(dim(final))
	#saveRDS(final, file = '1.3_GENIE3_weightMatrix_part_1.Rds')
	print(genes)
}

if( FALSE )
{
	final <- readRDS(file = '1.3_GENIE3_weightMatrix_part_1.Rds')
	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	regulator <- getDbTfs(scenicOptions)
	genesKept <- loadInt( scenicOptions, 'genesKept')
	exp <- readRDS(file = './working_data/sc_count_1-3d.rds')
	exp <- exp[genesKept,]
	exp <- log2(exp+1)
	regulator <- regulator[regulator %in% rownames(exp)]
	print(length(regulator))
	exp <- t(exp)
	exp <- as.data.frame(exp)
	print(dim(exp))
	write.table( regulator, file = paste0('regulator.csv'), quote = F,row.names = T,col.names = T,sep = '\t')
	target <- setdiff( colnames(exp), colnames(final))
	print(length(target))
	output <- union( regulator, target)
	print(length(output))
	i<-1
	write.table( target, file = paste0('target_test_', toString(i),'.csv'), quote = F,row.names = T,col.names = T,sep = '\t')
	library(feather)
	write_feather( exp[,output], paste0('exp_test_', toString(i),'.feather'))
	if( FALSE )
	{
		increment <- 20
		for( i in c(1:1) )
		{
			start <- (i-1)*increment + 1
			end <- start + increment - 1
			if( end > length(genekept) )
				end <- length(genekept)
			target <- genekept[c(start:end)]
			print(i)
			print(length(target))
			output <- union( regulator, target)
			write.table( target, file = paste0('target_test_', toString(i),'.csv'), quote = F,row.names = T,col.names = T,sep = '\t')
			write_feather( exp[,output], paste0('exp_test_', toString(i),'.feather'))
		}
	}
}

if( FALSE )
{
	library(SCENIC)
	final <- readRDS(file = '1.3_GENIE3_weightMatrix_part_1.Rds')
	scenicOptions <- readRDS(file = 'int/scenicOptions.Rds')
	genesKept <- loadInt( scenicOptions, 'genesKept')
	mat <- read.csv( paste0('genie3_1-3d_weightmatrix_test_1.csv'), header = TRUE, row.names = 1, sep = ',')
	colnames(mat) <- c('IQCJ-SCHIP1','NKX2-5','HLA-F','HLA-A','HLA-E','HLA-C','HLA-B','NKX3-1','BDNF-AS')
	print(intersect(genesKept,colnames(mat)))
	final <- cbind(final, mat)
	saveRDS( final, file = paste0("int/1.3_GENIE3_weightMatrix_part_1.Rds"))
}

#export regulon
if( FALSE )
{
	library(SCENIC)
	scenicOptions <- readRDS(file = 'int/scenicOptions.rds')
	temp <- loadInt(scenicOptions, 'aucell_regulonAUC')
	temp <- as.data.frame(temp@assays@data@listData$AUC)
	saveRDS(temp, file = 'regulon_auc.rds')
}

#calculate ferro sensitivity using AUC
if( FALSE )
{
	exp <- readRDS(file = './working_data/sc_count_1-3d.rds')
	print(dim(exp))
	print('finish loading')
	set.seed(123)
	library(SCENIC)
	library(AUCell)
	scenicOptions <- readRDS(file = 'int/scenicOptions.rds')
	genesKept <- loadInt( scenicOptions, 'genesKept')
	exp_filter_tmp <- exp[genesKept,]
	exp_filter <- log2(exp_filter_tmp+1)
	rm(exp)
	rm(exp_filter_tmp)
	cell_rankings <- AUCell_buildRankings(exp_filter, nCores = 10)
	saveRDS(cell_rankings, file = 'working_data/cell_rank')
}

if(FALSE)
{
	library(AUCell)
	anti <- read.table('data/antiferroptosis.txt',sep = '\t',header = TRUE,stringsAsFactors = FALSE)[,1]
	pro <- read.table('data/proferroptosis.txt',sep = '\t',header = TRUE,stringsAsFactors = FALSE)[,1]
	gene <- list(anti = anti, pro = pro)
	cell_rankings <- readRDS('working_data/cell_rank')
	cells_AUC <- AUCell_calcAUC(gene, cell_rankings, aucMaxRank = 5000, nCores = 10)
	saveRDS(cells_AUC, file = 'working_data/auc_ferro')
}

if(FALSE)
{
	library(AUCell)
	obj <- readRDS('working_data/auc_ferro')
    	mat <- getAUC(obj)[,]
	mat <- t(mat)
	mat <- data.frame(mat)
	mat$all <- mat$anti - mat$pro
	print(head(mat))
	saveRDS(mat, file = 'working_data/auc_ferro_mat')
}

#plot branched on NMF
if(FALSE)
{

	time <- readRDS(file = 'time')
	print(dim(time))
	print(head(rownames(time)))
	row_dist <- as.dist((1-cor(Matrix::t(time)))/2)
	time <- (time- rowMeans(time)) / apply(time,1,sd)
	print(time[,100])
	time[time>3] <- 3
	time[time<(-3)] <- (-3)
	#time[,1:100] <- (time[,1:100] - rowMeans(time[,1:100])) / apply(time[,1:100],1,sd)
	#time[,101:200] <- (time[,101:200] - rowMeans(time[,101:200])) / apply(time[,101:200],1,sd)
	#NMF::aheatmap(time, breaks = seq( -3,3,length.out = 300), 
	NMF::aheatmap(time, 
		      Rowv = row_dist, 
		      hclustfun = 'ward', Colv = NA, filename = 'test_3.pdf')
}

#LCD cell lines
if(FALSE)
{
    library(ggplot2)
    library(ggrepel)
    mat <- read.csv( paste('working_data/gr_gsh_persister.txt'),sep='\t', row.names = 1);
    print(head(mat))
    #mat[,'LCD'] <- factor(mat[,'LCD'])
    #mat <- mat_all[(mat_all$c >= i),]
    pdf("figure/gr_gsh_persister.pdf", onefile = TRUE, width = 10, height = 9)

    graph <- ggplot(mat, aes(y = proliferate, x = glutathione.reduced)) +
        geom_smooth(method='lm')+
        #geom_point( aes(fill = LCD), colour = 'black', pch = 21, alpha = 0.6, size = 4) +
        geom_point( aes(fill = persister), colour = 'black', pch = 21, alpha = 0.6, size = 4) +
        #scale_fill_brewer(palette = 'Dark2')+
        scale_x_continuous(limits = c(4.5,7))+
        scale_y_continuous(limits = c(0,0.07))+
        #guides( fill=FALSE)+
        #geom_label_repel( aes(label = gene_pair))+
        labs()
    graph <- graph + theme_classic()+
            theme(legend.position = c(0.9,0.9))
    print(graph)
    #mat <- mat[mat[,'LCD'] == 'True',]
    mat <- mat[mat[,'persister'] == 'True',]
    graph <- ggplot(mat, aes(y = proliferate, x = glutathione.reduced)) +
        #geom_point( colour = 'black', pch = 21, alpha = 0.6, size = 4) +
        geom_point( alpha = 0.6, size = 4) +
        geom_smooth(method='lm')+
        #scale_fill_brewer(palette = 'Dark2')+
        scale_x_continuous(limits = c(4.5,7))+
        scale_y_continuous(limits = c(0,0.07))+
        #guides( fill=FALSE)+
        geom_label_repel( aes(label = name))+
        labs()
    graph <- graph + theme_classic()
    print(graph)
    dev.off()
}


#media plot previous
if(FALSE)
{
    library(ggplot2)
    library(ggrepel)
    mat <- read.csv( paste('working_data/gr_gsh_medium.csv'),sep=',', row.names = 1);
    print(head(mat))
    #mat[,'LCD'] <- factor(mat[,'LCD'])
    #mat <- mat_all[(mat_all$c >= i),]
    pdf("figure/gr_gsh_LCD.pdf", onefile = TRUE, width = 10, height = 9)

    graph <- ggplot(mat, aes(y = proliferation, x = GSH)) +
        #geom_smooth(method='lm')+
        geom_point( aes(fill = LCD), colour = 'black', pch = 21, alpha = 0.6, size = 4) +
        #geom_point( aes(fill = culture), colour = 'black', pch = 21, alpha = 0.7, size = 3) +
        #guides( fill=FALSE)+
        #geom_label_repel( aes(label = gene_pair))+
        labs()
    graph <- graph + theme_classic()+
            theme(legend.position = c(0.9,0.9))
    print(graph)
    dev.off()
}

#media plot previous
if(FALSE)
{
    library(ggplot2)
    library(ggrepel)
    mat <- read.csv( paste('working_data/gr_gsh_medium.csv'),sep=',', row.names = 1);
    print(head(mat))
    #mat[,'LCD'] <- factor(mat[,'LCD'])
    #mat <- mat_all[(mat_all$c >= i),]
    pdf("figure/gr_gsh_medium.pdf", onefile = TRUE, width = 10, height = 9)

    graph <- ggplot(mat, aes(y = proliferation, x = GSH)) +
        #geom_smooth(method='lm')+
        #geom_point( aes(fill = LCD), colour = 'black', pch = 21, alpha = 0.6, size = 4) +
        geom_point( aes(fill = culture), colour = 'black', pch = 21, alpha = 0.7, size = 3) +
        scale_fill_brewer(palette = 'Dark2')+
        scale_x_continuous(limits = c(4.5,7))+
        scale_y_continuous(limits = c(0,0.07))+
        #guides( fill=FALSE)+
        #geom_label_repel( aes(label = gene_pair))+
        labs()
    graph <- graph + theme_classic()+
            theme(legend.position = c(0.9,0.9))
    print(graph)
    #mat <- mat[mat[,'LCD'] == 'True',]
    for( i in c('other','RPMI + 10% FBS','DMEM + 10% FBS','EMEM + 10% FBS','RPMI + 20% FBS','McCoy\'s 5A + 10% FBS') )
    {
        mat_2 <- mat
        mat_2[mat[,'culture'] != i,'culture'] <- ''
        graph <- ggplot(mat_2, aes(y = proliferation, x = GSH)) +
            geom_point( aes(fill = culture), colour = 'black', pch = 21, alpha = 0.6, size = 3) +
            scale_fill_manual(values = c("black","red"))+
            scale_x_continuous(limits = c(4.5,7))+
            scale_y_continuous(limits = c(0,0.07))+
            labs()
        graph <- graph + theme_classic()+
            theme(legend.position = c(0.9,0.9))
        print(graph)
    }
    dev.off()
}

#media plots
if(FALSE)
{
    library(ggplot2)
    library(ggrepel)
    mat <- read.csv( paste('working_data/meta_media.csv'),sep='\t', row.names = 1);
    print(head(mat))
    #mat[,'LCD'] <- factor(mat[,'LCD'])
    #mat <- mat_all[(mat_all$c >= i),]
    pdf("figure/gr_gsh_medium_meta.pdf", onefile = TRUE, width = 10, height = 9)

    graph <- ggplot(mat, aes(y = proliferate, x = GSH)) +
        #geom_smooth(method='lm')+
        #geom_point( aes(fill = LCD), colour = 'black', pch = 21, alpha = 0.6, size = 4) +
        geom_point( aes(fill = culture), colour = 'black', pch = 21, alpha = 0.7, size = 3) +
        scale_fill_brewer(palette = 'Dark2')+
        scale_x_continuous(limits = c(4.5,7))+
        scale_y_continuous(limits = c(0,0.05))+
        #guides( fill=FALSE)+
        #geom_label_repel( aes(label = gene_pair))+
        labs()
    graph <- graph + theme_classic()+
            theme(legend.position = c(0.9,0.9))
    print(graph)
    #mat <- mat[mat[,'LCD'] == 'True',]
    col <- unique(mat[,'culture'])
    for( i in col )
    {
        mat_2 <- mat
        mat_2[mat[,'culture'] != i,'culture'] <- ''
        graph <- ggplot(mat_2, aes(y = proliferate, x = GSH)) +
            geom_point( aes(fill = culture), colour = 'black', pch = 21, alpha = 0.6, size = 3) +
            scale_fill_manual(values = c("black","red"))+
            scale_x_continuous(limits = c(4.5,7))+
            scale_y_continuous(limits = c(0,0.05))+
            labs()
        graph <- graph + theme_classic()+
            theme(legend.position = c(0.9,0.9))
        print(graph)
    }
    dev.off()
}

#FBS20 plots
if(FALSE)
{
    library(ggplot2)
    library(ggbeeswarm)
    mat <- read.csv( paste('working_data/RPMI_20FBS.txt'),sep='\t', row.names = 1);
    print(head(mat))
    mat[,'culture'] <- factor(mat[,'culture'])
    #mat <- mat_all[(mat_all$c >= i),]
    pdf("figure/RPMI_20FBS_boxplot.pdf", onefile = TRUE, width = 10, height = 9)

    graph <- ggplot(mat, aes(y = glutathione.reduced, x = culture, color = culture)) +
        geom_beeswarm()+
        geom_boxplot(alpha = 0, width = 0.3)+
        #geom_label_repel( aes(label = gene_pair))+
        labs()
    graph <- graph + theme_classic(base_size = 20)+
            theme(legend.position = c(0.9,0.9))
    print(graph)
    graph <- ggplot(mat, aes(y = glutathione.reduced, x = culture, color = culture)) +
        geom_beeswarm()+
        geom_violin(alpha = 0.4, width = 0.5)+
        #geom_label_repel( aes(label = gene_pair))+
        labs()
    graph <- graph + theme_classic(base_size = 20)+
            theme(legend.position = c(0.9,0.9))
    print(graph)
    dev.off()
}

#LCD cell lines new GSH & media
if(FALSE)
{
    library(ggplot2)
    library(ggrepel)
    mat <- read.csv( paste('table/LCD_cell_GSH_info_3.csv'),sep=',', row.names = 1);
    print(head(mat))
    #mat[,'LCD'] <- factor(mat[,'LCD'])
    #mat <- mat_all[(mat_all$c >= i),]
    pdf("figure/gr_gsh_LCD_2.pdf", onefile = TRUE, width = 13, height = 9)

    graph <- ggplot(mat, aes(y = proliferation, x = GSH)) +
        geom_smooth(method='lm')+
        #geom_point( aes(fill = LCD), colour = 'black', pch = 21, alpha = 0.6, size = 4) +
        geom_point( aes(fill = culture), colour = 'black', pch = 21, alpha = 0.6, size = 4) +
        #scale_fill_brewer(palette = 'Dark2')+
        scale_x_continuous(limits = c(4.5,7))+
        scale_y_continuous(limits = c(0,0.05))+
        #guides( fill=FALSE)+
        geom_text_repel( aes(label = cell), max.overlaps = 50)+
        labs()
    graph <- graph + theme_classic(base_size = 18)
    #graph <- graph + theme_classic(base_size = 18)+
    #        theme(legend.position = c(0.9,0.9))
    print(graph)
    
    dev.off()
}


#LCD cell lines
if(FALSE)
{
    library(ggplot2)
    library(ggrepel)
    mat <- read.csv( paste('working_data/gr_gsh_lcd_new.txt'),sep='\t', row.names = 1);
    print(head(mat))
    #mat[,'LCD'] <- factor(mat[,'LCD'])
    #mat <- mat_all[(mat_all$c >= i),]
    pdf("figure/gr_gsh_LCD_new.pdf", onefile = TRUE, width = 12, height = 9)

    graph <- ggplot(mat, aes(y = proliferate, x = glutathione.reduced)) +
        geom_point( aes(fill = LCD), colour = 'black', pch = 21, alpha = 0.6, size = 4) +
        #scale_fill_brewer(palette = 'Dark2')+
        scale_x_continuous(limits = c(4.75,6.75))+
        scale_y_continuous(limits = c(0,0.06))+
        #guides( fill=FALSE)+
        #geom_label_repel( aes(label = gene_pair))+
        labs()
    graph <- graph + theme_classic(base_size = 18)+
            theme(legend.position = c(0.9,0.9))
    print(graph)

    graph <- ggplot(mat, aes(y = xCT, x = glutathione.reduced)) +
        geom_point( aes(fill = LCD), colour = 'black', pch = 21, alpha = 0.6, size = 4) +
        #scale_fill_brewer(palette = 'Dark2')+
        scale_x_continuous(limits = c(4.75,6.75))+
        scale_y_continuous(limits = c(0,7.8))+
        #guides( fill=FALSE)+
        #geom_label_repel( aes(label = gene_pair))+
        labs()
    graph <- graph + theme_classic(base_size = 18)+
            theme(legend.position = c(0.9,0.9))
    print(graph)

    graph <- ggplot(mat, aes(y = xCT, x = SLC1A5)) +
        geom_point( aes(fill = LCD), colour = 'black', pch = 21, alpha = 0.6, size = 4) +
        #scale_fill_brewer(palette = 'Dark2')+
        scale_x_continuous(limits = c(4,11))+
        scale_y_continuous(limits = c(0,7.8))+
        #guides( fill=FALSE)+
        #geom_label_repel( aes(label = gene_pair))+
        labs()
    graph <- graph + theme_classic(base_size = 18)+
            theme(legend.position = c(0.9,0.9))
    print(graph)

    graph <- ggplot(mat, aes(y = xCT, x = SLC7A5)) +
        geom_point( aes(fill = LCD), colour = 'black', pch = 21, alpha = 0.6, size = 4) +
        #scale_fill_brewer(palette = 'Dark2')+
        scale_x_continuous(limits = c(2,10.5))+
        scale_y_continuous(limits = c(0,7.8))+
        #guides( fill=FALSE)+
        #geom_label_repel( aes(label = gene_pair))+
        labs()
    graph <- graph + theme_classic(base_size = 18)+
            theme(legend.position = c(0.9,0.9))
    print(graph)

    graph <- ggplot(mat, aes(y = xCT, x = SLC38A2)) +
        geom_point( aes(fill = LCD), colour = 'black', pch = 21, alpha = 0.6, size = 4) +
        #scale_fill_brewer(palette = 'Dark2')+
        scale_x_continuous(limits = c(4,11))+
        scale_y_continuous(limits = c(0,7.8))+
        #guides( fill=FALSE)+
        #geom_label_repel( aes(label = gene_pair))+
        labs()
    graph <- graph + theme_classic(base_size = 18)+
            theme(legend.position = c(0.9,0.9))
    print(graph)
    mat <- mat[mat[,'LCD'] == 'True',]
    #mat <- mat[mat[,'persister'] == 'True',]
    graph <- ggplot(mat, aes(y = proliferate, x = glutathione.reduced)) +
        #geom_point( colour = 'black', pch = 21, alpha = 0.6, size = 4) +
        geom_point( alpha = 0.6, size = 4) +
        #scale_fill_brewer(palette = 'Dark2')+
        scale_x_continuous(limits = c(4.75,6.75))+
        scale_y_continuous(limits = c(0,0.06))+
        #guides( fill=FALSE)+
        geom_label_repel( aes(label = name), max.overlaps = 50)+
        labs()
    graph <- graph + theme_classic(base_size = 18)
    print(graph)
    graph <- ggplot(mat, aes(y = xCT, x = glutathione.reduced)) +
        #geom_point( colour = 'black', pch = 21, alpha = 0.6, size = 4) +
        geom_point( alpha = 0.6, size = 4) +
        #scale_fill_brewer(palette = 'Dark2')+
        scale_x_continuous(limits = c(4.75,6.75))+
        scale_y_continuous(limits = c(0,7.8))+
        #guides( fill=FALSE)+
        geom_label_repel( aes(label = name), max.overlaps = 50)+
        labs()
    graph <- graph + theme_classic(base_size = 18)
    print(graph)

    graph <- ggplot(mat, aes(y = xCT, x = SLC1A5)) +
        geom_point( aes(fill = LCD), colour = 'black', pch = 21, alpha = 0.6, size = 4) +
        #scale_fill_brewer(palette = 'Dark2')+
        scale_x_continuous(limits = c(4,11))+
        scale_y_continuous(limits = c(0,7.8))+
        geom_label_repel( aes(label = name), max.overlaps = 50)+
        #guides( fill=FALSE)+
        #geom_label_repel( aes(label = gene_pair))+
        labs()
    graph <- graph + theme_classic(base_size = 18)+
            theme(legend.position = c(0.9,0.9))
    print(graph)

    graph <- ggplot(mat, aes(y = xCT, x = SLC7A5)) +
        geom_point( aes(fill = LCD), colour = 'black', pch = 21, alpha = 0.6, size = 4) +
        #scale_fill_brewer(palette = 'Dark2')+
        scale_x_continuous(limits = c(2,10.5))+
        scale_y_continuous(limits = c(0,7.8))+
        geom_label_repel( aes(label = name), max.overlaps = 50)+
        #guides( fill=FALSE)+
        #geom_label_repel( aes(label = gene_pair))+
        labs()
    graph <- graph + theme_classic(base_size = 18)+
            theme(legend.position = c(0.9,0.9))
    print(graph)

    graph <- ggplot(mat, aes(y = xCT, x = SLC38A2)) +
        geom_point( aes(fill = LCD), colour = 'black', pch = 21, alpha = 0.6, size = 4) +
        #scale_fill_brewer(palette = 'Dark2')+
        scale_x_continuous(limits = c(4,11))+
        geom_label_repel( aes(label = name), max.overlaps = 50)+
        scale_y_continuous(limits = c(0,7.8))+
        #guides( fill=FALSE)+
        #geom_label_repel( aes(label = gene_pair))+
        labs()
    graph <- graph + theme_classic(base_size = 18)+
            theme(legend.position = c(0.9,0.9))
    print(graph)
    dev.off()
}

#LCD cell lines
if(FALSE)
{
    library(ggplot2)
    library(ggrepel)
    mat <- read.csv( paste('working_data/gr_gsh_lcd_2.txt'),sep='\t', row.names = 1);
    print(head(mat))
    #mat[,'LCD'] <- factor(mat[,'LCD'])
    #mat <- mat_all[(mat_all$c >= i),]
    pdf("figure/gr_gsh_LCD_22.pdf", onefile = TRUE, width = 10, height = 9)

    graph <- ggplot(mat, aes(y = proliferation, x = glutathione.reduced)) +
        geom_smooth(method='lm')+
        #geom_point( aes(fill = LCD), colour = 'black', pch = 21, alpha = 0.6, size = 4) +
        geom_point( aes(fill = LCD), colour = 'black', pch = 21, alpha = 0.6, size = 4) +
        #scale_fill_brewer(palette = 'Dark2')+
        scale_x_continuous(limits = c(4.5,7))+
        scale_y_continuous(limits = c(0,0.05))+
        #guides( fill=FALSE)+
        #geom_label_repel( aes(label = gene_pair))+
        labs()
    graph <- graph + theme_classic(base_size = 18)+
            theme(legend.position = c(0.9,0.9))
    print(graph)
    mat <- mat[mat[,'LCD'] == 'True',]
    #mat <- mat[mat[,'persister'] == 'True',]
    graph <- ggplot(mat, aes(y = proliferation, x = glutathione.reduced)) +
        #geom_point( colour = 'black', pch = 21, alpha = 0.6, size = 4) +
        geom_point( alpha = 0.6, size = 4) +
        geom_smooth(method='lm')+
        #scale_fill_brewer(palette = 'Dark2')+
        scale_x_continuous(limits = c(4.5,7))+
        scale_y_continuous(limits = c(0,0.05))+
        #guides( fill=FALSE)+
        geom_label_repel( aes(label = name), max.overlaps = 50)+
        labs()
    graph <- graph + theme_classic(base_size = 18)
    print(graph)
    dev.off()
}

#media plots_2
if(FALSE)
{
    library(ggplot2)
    library(ggrepel)
    mat <- read.csv( paste('table/LCD_cell_GSH_info_3.csv'),sep=',', row.names = 1);
    print(head(mat))
    #mat[,'LCD'] <- factor(mat[,'LCD'])
    #mat <- mat_all[(mat_all$c >= i),]
    pdf("figure/gr_gsh_medium_LCD_meta.pdf", onefile = TRUE, width = 10, height = 9)

    if(FALSE)
    {
    graph <- ggplot(mat, aes(y = proliferation, x = GSH)) +
        #geom_smooth(method='lm')+
        #geom_point( aes(fill = LCD), colour = 'black', pch = 21, alpha = 0.6, size = 4) +
        geom_point( aes(fill = culture), colour = 'black', pch = 21, alpha = 0.7, size = 3) +
        scale_fill_brewer(palette = 'Dark2')+
        scale_x_continuous(limits = c(4.5,7))+
        scale_y_continuous(limits = c(0,0.05))+
        #guides( fill=FALSE)+
        #geom_label_repel( aes(label = gene_pair))+
        labs()
    graph <- graph + theme_classic(base_size = 18)+
            theme(legend.position = c(0.9,0.9))
    print(graph)
    }
    #mat <- mat[mat[,'LCD'] == 'True',]
    col <- unique(mat[,'culture'])
    for( i in col )
    {
        if( i != '')
        {
        mat_2 <- mat
        mat_2[mat[,'culture'] != i,'culture'] <- ''
        mat_2[mat[,'culture'] != i,'cell'] <- ''
        graph <- ggplot(mat_2, aes(y = proliferation, x = GSH)) +
            geom_point( aes(fill = culture), colour = 'black', pch = 21, alpha = 0.6, size = 3) +
            scale_fill_manual(values = c("black","red"))+
            scale_x_continuous(limits = c(4.5,7))+
            scale_y_continuous(limits = c(0,0.05))+
            geom_text_repel( aes(label = cell))+
            labs()
        graph <- graph + theme_classic(base_size = 18)+
            theme(legend.position = c(0.2,0.9))
        print(graph)
        }
    }
    dev.off()
}

#persister
if(FALSE)
{
	library(Seurat)
	counts <- read.csv(file = 'data/persister_data.csv', header = TRUE, row.names = 1)
	print('Finish reading')
	sc <- CreateSeuratObject( counts = counts, project = 'persister')
	print('Finish creating')
	saveRDS(sc, file = './working_data/persister_create.rds')
	#sc[['RNA']]@counts
	#sc[['percent.mt']] <- PercentageFeatureSet(sc, pattern = '^MT-')
	#cell_num <- table(strsplit2(rownames(sc@meta.data),split = '-')[,2])
	#attributes(cell_num)$dimnames[[1]] <- c('M4340','M4341','M4342','M4343')
	#sample_name <- ulist(apply(as.data.frame(cell_num), 1, function(x){rep(x[1], x[2])}))
	#sc <- AddMetaData( object = sc, metadata = sample_name, col.name = 'sample')
	#sc <- subset(sc, subset = percent.mt < 10 & nCount_RNA > 5000)
	#sc <- NormalizeData(sc, normalization.method = 'LogNormalize', scale.factor = 10000)
	#sc[['RNA']]@data
	#saveRDS(sc, file = './working_data/sc_normal.rds')
	#all.genes <- rownames(sc)
	#sc <- ScaleData(sc, features = all.genes)
	#saveRDS(sc, file = './working_data/sc_scale.rds')
	#print('Finish scaling')
}

if(FALSE)
{
	library(Seurat)
	info <- read.csv(file = 'data/persister_info.csv', header = TRUE, row.names = 1)
	print(head(info))
	sc = readRDS(file = './working_data/persister_create.rds')
	name <- colnames( x = sc)
	name <- gsub('\\.','-',name)
	sc <- RenameCells( object = sc, new.names = name)
	sc <- AddMetaData( object = sc, metadata = info)
	print(head(sc[[]]))
	saveRDS(sc, file = './working_data/persister_comb.rds')
	print('Finish combining')
}


if(FALSE)
{
	library(Seurat)
	sc = readRDS(file = './working_data/persister_comb.rds')
	sc <- FindVariableFeatures( sc, selection.method = 'vst', nFeatures = 2000)
	all.genes <- rownames(sc)
	sc <- ScaleData(sc, features = all.genes)
	sc <- RunPCA( sc, features = VariableFeatures(object = sc))
	sc <- RunUMAP( sc, dims = 1:10)
	saveRDS( sc, file = './working_data/persister_umap.rds')
	DimPlot(sc, reduction = 'umap', group.by = 'cell_line')
	DimPlot(sc, reduction = 'umap', group.by = 'sample_type')
}

if(FALSE)
{
	library(Seurat)
	sc = readRDS(file = './working_data/persister_umap.rds')
	pdf("figure/YPEL3.pdf", onefile = TRUE)
	plt <- DimPlot(sc, reduction = 'umap', group.by = 'cell_line')
	print(plt)
	plt <- DimPlot(sc, reduction = 'umap', group.by = 'sample_type')
	print(plt)
	plt <- FeaturePlot(sc, features = c('YPEL3'))
	print(plt)
	plt <- FeaturePlot(sc, features = c('FOXO3'))
	print(plt)
	dev.off()
}

if(FALSE)
{
	library(Seurat)
	sc = readRDS(file = './working_data/persister_comb.rds')
	#print(sc)
	#print(head(sc@assays$RNA@counts))
	#print(sc@assays)
	#print(Idents(object = sc))
	#sc <- RegroupIdents(sc, metadata = 'cell_line')
	#s
	#Idents(sc) <- 'cell_line'
	#levels(sc)
	#print(table(Idents(sc)))
	#sc <- subset(sc, subset = cell_line == 'A375_SKIN')
	#print(sc)
	Idents(sc) <- 'sample_type'
	#print(table(sc$sample_type))
	#a375.markers <- FindMarkers( sc, ident.1 = 'naive', ident.2 = 'non-cycling')
	#a375.markers <- FindMarkers( sc, ident.1 = 'non-cycling', ident.2 = 'naive', min.pct = 0.01)
	for( i in c('A375_SKIN','PC9_LUNG','HT29_LARGE_INTESTINE','COLO829_SKIN'))
	{
		print(i)
		tmp <- subset(sc, subset = cell_line == i)
		markers <- FindMarkers( tmp, ident.1 = 'non-cycling', ident.2 = 'naive', min.pct = 0, test.use = 'MAST')
		write.table( markers, file = paste0('working_data/',i,'_naive.csv'), quote = F, row.names = T, col.names = T, sep = ',')
		markers <- FindMarkers( tmp, ident.1 = 'non-cycling', ident.2 = 'cycling', min.pct = 0, test.use = 'MAST')
		write.table( markers, file = paste0('working_data/',i,'_cycling.csv'), quote = F, row.names = T, col.names = T, sep = ',')
	}
	for( i in c('HCC827_LUNG'))
	{
		print(i)
		tmp <- subset(sc, subset = cell_line == i)
		markers <- FindMarkers( tmp, ident.1 = 'non-cycling', ident.2 = 'cycling', min.pct = 0, test.use = 'MAST')
		write.table( markers, file = paste0('working_data/',i,'_cycling.csv'), quote = F, row.names = T, col.names = T, sep = ',')
	}
}

if(FALSE)
{
	library(Seurat)
	sc = readRDS(file = './working_data/persister_comb.rds')
	#print(sc)
	#print(head(sc@assays$RNA@counts))
	#print(sc@assays)
	#print(Idents(object = sc))
	#sc <- RegroupIdents(sc, metadata = 'cell_line')
	#s
	#Idents(sc) <- 'cell_line'
	#levels(sc)
	#print(table(Idents(sc)))
	#sc <- subset(sc, subset = cell_line == 'A375_SKIN')
	#sc <- subset(sc, subset = cell_line == 'BT474_BREAST')
	#print(sc)
	Idents(sc) <- 'cell_line'
	pdf("figure/YPEL3.pdf", onefile = TRUE)
	plt <- VlnPlot(sc, features = c('YPEL3'))
	print(plt)
	dev.off()
}

#plotting PCA with only a few genes
if(FALSE)
{
	library(Seurat)
	sc = readRDS(file = './working_data/persister_comb.rds')
	all.genes <- rownames(sc)
	sc <- ScaleData(sc, features = all.genes)
    saveRDS( sc, file = './working_data/persister_scale.rds')
}

if(FALSE)
{
    mat <- read.table( paste('working_data/top_cor_gene.txt'),sep='\t', header = TRUE);
    print(head(mat))
}

if(FALSE)
{
	library(Seurat)
	sc <- readRDS(file = './working_data/persister_comb.rds')
    mat <- read.table( paste('working_data/top_cor_gene.txt'),sep='\t', header = TRUE);
    print(head(mat))
    set.seed(20000)
	#print(sc)
	#print(head(sc@assays$RNA@counts))
	#print(sc@assays)
	#print(Idents(object = sc))
	#sc <- RegroupIdents(sc, metadata = 'cell_line')
	#s
	#Idents(sc) <- 'cell_line'
	#levels(sc)
	#print(table(Idents(sc)))
	sc <- subset(sc, subset = cell_line == 'A375_SKIN')
	#print(sc)
	Idents(sc) <- 'sample_type'
    #sc <- subset(sc, subset = cell_line == 'HT29_LARGE_INTESTINE')
    #sc <- subset(sc, subset = cell_line == 'PC9_LUNG')
    genes <- mat[,'cell']
	#genes <- rownames(sc)[sample(1:length(rownames(sc)), 209)]
    sc <- ScaleData(sc, feature = genes)
	#print(table(sc$sample_type))
	#a375.markers <- FindMarkers( sc, ident.1 = 'naive', ident.2 = 'non-cycling')
	#a375.markers <- FindMarkers( sc, ident.1 = 'non-cycling', ident.2 = 'naive', min.pct = 0.01)
    #sc <- subset(sc, subset = cell_line == 'PC9_LUNG')
    print('success')
    sc <- RunPCA( sc, features = genes)
    sc <- RunUMAP( sc, dims = 1:10)
	pdf("figure/scrna_test_A375_multiple_reg.pdf", onefile = TRUE)
	plt <- DimPlot(sc, reduction = 'umap',  cols = c("darkgrey","coral", "cornflowerblue"))
	print(plt)
	plt <- PCAPlot(sc, cols = c("darkgrey","coral", "cornflowerblue"))
	print(plt)
	dev.off()

    if(FALSE)
    {
	for( i in c('A375_SKIN','PC9_LUNG','HT29_LARGE_INTESTINE','COLO829_SKIN'))
	{
		print(i)
		tmp <- subset(sc, subset = cell_line == i)
		markers <- FindMarkers( tmp, ident.1 = 'non-cycling', ident.2 = 'naive', min.pct = 0, test.use = 'MAST')
		write.table( markers, file = paste0('working_data/',i,'_naive.csv'), quote = F, row.names = T, col.names = T, sep = ',')
		markers <- FindMarkers( tmp, ident.1 = 'non-cycling', ident.2 = 'cycling', min.pct = 0, test.use = 'MAST')
		write.table( markers, file = paste0('working_data/',i,'_cycling.csv'), quote = F, row.names = T, col.names = T, sep = ',')
	}
	for( i in c('HCC827_LUNG'))
	{
		print(i)
		tmp <- subset(sc, subset = cell_line == i)
		markers <- FindMarkers( tmp, ident.1 = 'non-cycling', ident.2 = 'cycling', min.pct = 0, test.use = 'MAST')
		write.table( markers, file = paste0('working_data/',i,'_cycling.csv'), quote = F, row.names = T, col.names = T, sep = ',')
	}
    }
}

if(FALSE)
{
    dep <- readRDS('working_data/dep.rds')
    mat <- read.table( paste('working_data/dep_index.csv'),sep=',', header = TRUE);
    rownames(dep) <- mat[,'X']
    saveRDS(dep, file = './working_data/dep_all.rds')
    
    dep <- readRDS('working_data/pred.rds')
    mat <- read.table( paste('working_data/pred_index.csv'),sep=',', header = TRUE);
    rownames(dep) <- mat[,'X']
    saveRDS(dep, file = './working_data/pred_all.rds')
    
    dep <- readRDS('working_data/exp.rds')
    mat <- read.table( paste('working_data/exp_index.csv'),sep=',', header = TRUE);
    rownames(dep) <- mat[,'X']
    saveRDS(dep, file = './working_data/exp_all.rds')
}

if(FALSE)
{
    dep <- readRDS('working_data/dep_all.rds')
    pred <- readRDS('working_data/pred_all.rds')
    print(head(pred))
    pred <- pred[rownames(dep),]
    print(dim(pred))
    print(head(dep[,c('TP53','MDM2')]))
}

if(FALSE)
{
    library(glmnet)
    dep <- readRDS('working_data/exp_all.rds')
    pred <- readRDS('working_data/pred_all.rds')
    pred <- pred[rownames(dep),]
    mfit <- glmnet( dep, pred, family = 'mgaussian')
    print(mfit)
    plot(mfit)
}

if(FALSE)
{
    library(glmnet)
    dep <- readRDS('working_data/exp_all.rds')
    pred <- readRDS('working_data/pred_all.rds')
    pred <- pred[rownames(dep),]
    mfit <- glmnet( dep, pred, family = 'mgaussian')
    #print(mfit)
    #plot(mfit)
    tmp <- coef(mfit, s = 0.08)
    #pro <- as.matrix( tmp@'proliferate'])
    print(str(tmp))
    gsh <- tmp$GSH
    print(gsh)
    print(gsh@Dimnames[[1]][gsh@i])
    #tmp <- as.data.frame(summary(tmp))
    #print(tmp)
    #print(head(pro))
    #write.table(coef(mfit, s = 0.08), file = 'working_data/elastic_test.csv', sep = ',', col.names = T, row.names = T, quote = F)
}

if(FALSE)
{
    library(glmnet)
    set.seed(20000)
    dep <- readRDS('working_data/exp_all.rds')
    pred <- readRDS('working_data/pred_all.rds')
    pred <- pred[rownames(dep),]
    mfit <- glmnet( dep, pred, family = 'mgaussian')
    #print(mfit)
    #plot(mfit)
    tmp <- coef(mfit, s = 0.08, alpha = 0.5)
    print(str(tmp))
    gsh <- tmp$GSH
    pro <- tmp$proliferate
    #print(gsh@Dimnames[[1]][gsh@i])
    #print(pro@Dimnames[[1]][pro@i])

    library(Seurat)
	sc <- readRDS(file = './working_data/persister_comb.rds')
    Idents(sc) <- 'sample_type'
    sc <- subset(sc, subset = cell_line == 'PC9_LUNG')
    genes <- pro@Dimnames[[1]][pro@i]
	#genes <- rownames(sc)[sample(1:length(rownames(sc)), 209)]
    sc <- ScaleData(sc, feature = genes)
    print('success')
    sc <- RunPCA( sc, features = genes)
    sc <- RunUMAP( sc, dims = 1:10)
	pdf("figure/scrna_test_PC9_elastic_2.pdf", onefile = TRUE)
	plt <- DimPlot(sc, reduction = 'umap',  cols = c("darkgrey","coral", "cornflowerblue"))
	print(plt)
	plt <- PCAPlot(sc, cols = c("darkgrey","coral", "cornflowerblue"))
	print(plt)
	dev.off()
}

if(FALSE)
{
    library(glmnet)
    set.seed(20000)
    dep <- readRDS('working_data/exp_all.rds')
    pred <- readRDS('working_data/pred_all.rds')
    dep <- as.matrix(dep)
    pred <- as.matrix(pred)
    pred <- pred[rownames(dep),]
    #dep <- dep[,c('YPEL3','TKT','GCLC','TP53','MDM2','CCND1','PGD','SCD','GPX4')]
    print(dim(dep))
    mfit <- glmnet( dep, pred, family = 'mgaussian')
    #print(mfit)
    #ass <- predict(mfit, newx = dep[,], type = "response", s = 0.281)
    tmp <- coef(mfit, s = 0.07)
    print(str(tmp))
    #gsh <- tmp$GSH
    #pro <- tmp$proliferate
    gsh <- data.frame(tmp$GSH@x)
    rownames(gsh) <- c('intercept',tmp$GSH@Dimnames[[1]][tmp$GSH@i])
    colnames(gsh) <- c('GSH')
    pro <- data.frame(tmp$proliferate@x)
    rownames(pro) <- c('intercept',tmp$proliferate@Dimnames[[1]][tmp$proliferate@i])
    colnames(pro) <- c('pro')
    result <- cbind(gsh, pro)
    #saveRDS(result, file = 'working_data/elastic_coef_topy.rds')
    #result <- t(result)
    print(head(result))
    #print(gsh@Dimnames[[1]][gsh@i])
    #print(pro@Dimnames[[1]][pro@i])
    ass <- predict(mfit, newx = dep[,], type = "response", s = 0.07)
    ass <- ass[,,1]
    diff <- ass - pred
    #print(head(ass))
    #print(head(pred))
    #print(head(diff))
    diff <- diff ** 2
    #print(dim(pred)[1])
    diff <- sum(colSums(diff)) / 2 / dim(pred)[1]
    print(diff)
    ass <- assess.glmnet(mfit, newx = dep[,], newy = pred, s = 0.07)
    print(ass)
    ass <- assess.glmnet(mfit, newx = dep[,], newy = pred)
    print(ass)
}

if(FALSE)
{
	library(Seurat)
	sc = readRDS(file = './working_data/persister_comb.rds')
    sc <- subset(sc, subset = cell_line == 'PC9_LUNG')
	sc <- FindVariableFeatures( sc, selection.method = 'vst', nFeatures = 200)
    genes <- VariableFeatures(object = sc)
	Idents(sc) <- 'sample_type'
    sc <- ScaleData(sc, feature = genes)
    print('success')
    sc <- RunPCA( sc, features = genes)
    sc <- RunUMAP( sc, dims = 1:10)
	pdf("figure/scrna_test_PC9_positive_control.pdf", onefile = TRUE)
	plt <- DimPlot(sc, reduction = 'umap',  cols = c("darkgrey","coral", "cornflowerblue"))
	print(plt)
	plt <- PCAPlot(sc, cols = c("darkgrey","coral", "cornflowerblue"))
	print(plt)
	dev.off()

    if(FALSE)
    {
        for( i in c('A375_SKIN','PC9_LUNG','HT29_LARGE_INTESTINE','COLO829_SKIN'))
        {
        }
    }
}


if(FALSE)
{
    library(glmnet)
    dep <- readRDS('working_data/exp_all.rds')
    pred <- readRDS('working_data/pred_all.rds')
    pred <- pred[rownames(dep),]
    dep <- dep[rownames(dep),]
    pred <- as.matrix(pred)
    dep <- as.matrix(dep)
    print(dim(pred))
    print(dim(dep))
    print(head(pred))

    #data(QuickStartExample)
    #dep <- QuickStartExample$x
    #pred <- cbind(QuickStartExample$y, QuickStartExample$y)
    #print(dim(pred))
    #print(dim(dep))
    #print(head(pred))
    #print(head(dep))
    #print(head(pred))
    set.seed(20000)
    mfit <- cv.glmnet( x = dep, y = pred, family = 'mgaussian', alpha = 0.05, trace.it = 1)
    #print(mfit)
    pdf("figure/elastic_net_005.pdf", onefile = TRUE)
    plt <- plot(mfit)
	print(plt)
    print(mfit$lambda.min)
    print(mfit$lambda.1se)
	dev.off()
}

if(FALSE)
{
    library(glmnet)
    dep <- readRDS('working_data/exp_all.rds')
    pred <- readRDS('working_data/pred_all.rds')
    pred <- pred[rownames(dep),]
    pred <- as.matrix(pred)
    dep <- as.matrix(dep)
    print(dim(pred))
    print(dim(dep))

    set.seed(20000)
    #mfit <- glmnet( x = dep, y = pred, family = 'mgaussian')
    mfit <- cv.glmnet( x = dep, y = pred, family = 'mgaussian')
    print(mfit)
    saveRDS(mfit, file = 'working_data/elastic_net_alpha_1_test.rds')
    #s <- 0.1199737
    #result <- assess.glmnet(mfit, newx = dep, newy = pred, s = 0.1199737)
    #result <- predict(mfit, newx = dep, s = 0.1199737)
    #print(result)
    #print(str(result))
    #print(mfit)
    #pdf("figure/elastic_net_test.pdf", onefile = TRUE)
    #plt <- plot(mfit)
	#print(plt)
    #print(mfit$lambda.min)
    #print(mfit$lambda.1se)
	#dev.off()
}

#simulation of alpha and lambda selection
if(FALSE)
{
    library(glmnet)
    data(QuickStartExample)
    dep <- QuickStartExample$x
    pred <- cbind(QuickStartExample$y, QuickStartExample$y)
    print(dim(pred))
    print(dim(dep))
    set.seed(20000)
    alpha_grid <- c(0,0.1,0.25,0.5,0.75,0.9,1)

    nfolds <- 10
    train_idx <- sample.int(n = nrow(dep), size = floor(.75*nrow(dep)), replace = F)
    x_train <- dep[train_idx,]
    y_train <- pred[train_idx,]
    x_test <- dep[-train_idx,]
    y_test <- pred[-train_idx,]
    print('finish train test split')

    foldid <- sample(rep(seq(nfolds), length = nrow(x_train)))
    print(foldid)
    result <- as.data.frame(matrix(nrow = length(alpha_grid), ncol = nfolds))
    colnames(result) <- seq(nfolds)
    rowname <- c()
    for( i in alpha_grid)
    {
        rowname <- c(rowname, toString(i))
    }
    rownames(result) <- rowname
    print(result)
    print('finish ID')
    library(doMC)
    para <- 1
    registerDoMC(10)
    if( para == 1 )
    {
        outlist <- foreach( i = seq(nfolds), .packages = c('glmnet')) %dopar%
        {
            which <- foldid == i
            x_sub <- x_train[!which,]
            y_sub <- y_train[!which,]
            x_sub_test <- x_train[which,]
            y_sub_test <- y_train[which,]
            #print('Finish External split')
            result <- c()
            for( j in alpha_grid)
            {
                #print(dim(x_sub))
                #print(dim(y_sub))
                #print(j)
                mfit <- cv.glmnet( x = x_sub, y = y_sub, family = 'mgaussian', alpha =j)
                #print('finish fit')
                assess <- assess.glmnet(mfit, newx = x_sub_test, newy = y_sub_test, s = 'lambda.min')
                #print('finish assess')
                result <- c(result,assess$mse[[1]])
            }
            return(result)
        }
        outlist <- as.data.frame(do.call(cbind, outlist))
        #print(outlist)
        result <- outlist
        rownames(result) <- rowname
        colnames(result) <- seq(nfolds)
    }
    else
    {
        for( i in seq(nfolds) )
        {
            print(i)
            which <- foldid == i
            x_sub <- x_train[!which,]
            y_sub <- y_train[!which,]
            x_sub_test <- x_train[which,]
            y_sub_test <- y_train[which,]
            #print('Finish External split')
            for( j in alpha_grid)
            {
                #print(dim(x_sub))
                #print(dim(y_sub))
                print(j)
                mfit <- cv.glmnet( x = x_sub, y = y_sub, family = 'mgaussian', alpha =j)
                #print('finish fit')
                assess <- assess.glmnet(mfit, newx = x_sub_test, newy = y_sub_test, s = 'lambda.min')
                #print('finish assess')
                result[toString(j),i] <- assess$mse[[1]]
            }
        }
    }
    print(result)
    saveRDS(result, file = 'working_data/glmnet_test.rds')
}


#simulation of alpha and lambda selection
if(FALSE)
{
    library(glmnet)
    dep <- readRDS('working_data/exp_all.rds')
    pred <- readRDS('working_data/pred_all.rds')
    pred <- pred[rownames(dep),]
    pred <- as.matrix(pred)
    dep <- as.matrix(dep)
    print(dim(pred))
    print(dim(dep))
    set.seed(20000)
    #alpha_grid <- c(0,0.1,0.25,0.5,0.75,0.9,1)
    alpha_grid <- c(0,0.01,0.025,0.05,0.075,0.09,0.1)

    nfolds <- 10
    train_idx <- sample.int(n = nrow(dep), size = floor(.75*nrow(dep)), replace = F)
    x_train <- dep[train_idx,]
    y_train <- pred[train_idx,]
    x_test <- dep[-train_idx,]
    y_test <- pred[-train_idx,]
    print('finish train test split')

    foldid <- sample(rep(seq(nfolds), length = nrow(x_train)))
    print(foldid)
    result <- as.data.frame(matrix(nrow = length(alpha_grid), ncol = nfolds))
    colnames(result) <- seq(nfolds)
    rowname <- c()
    for( i in alpha_grid)
    {
        rowname <- c(rowname, toString(i))
    }
    rownames(result) <- rowname
    print(result)
    print('finish ID')
    library(doMC)
    para <- 1
    registerDoMC(10)
    if( para == 1 )
    {
        outlist <- foreach( i = seq(nfolds), .packages = c('glmnet')) %dopar%
        {
            which <- foldid == i
            x_sub <- x_train[!which,]
            y_sub <- y_train[!which,]
            x_sub_test <- x_train[which,]
            y_sub_test <- y_train[which,]
            #print('Finish External split')
            result <- c()
            for( j in alpha_grid)
            {
                #print(dim(x_sub))
                #print(dim(y_sub))
                #print(j)
                mfit <- cv.glmnet( x = x_sub, y = y_sub, family = 'mgaussian', alpha =j)
                #print('finish fit')
                assess <- assess.glmnet(mfit, newx = x_sub_test, newy = y_sub_test, s = 'lambda.min')
                #print('finish assess')
                result <- c(result,assess$mse[[1]])
            }
            return(result)
        }
        outlist <- as.data.frame(do.call(cbind, outlist))
        #print(outlist)
        result <- outlist
        rownames(result) <- rowname
        colnames(result) <- seq(nfolds)
    }
    else
    {
        for( i in seq(nfolds) )
        {
            print(i)
            which <- foldid == i
            x_sub <- x_train[!which,]
            y_sub <- y_train[!which,]
            x_sub_test <- x_train[which,]
            y_sub_test <- y_train[which,]
            #print('Finish External split')
            for( j in alpha_grid)
            {
                #print(dim(x_sub))
                #print(dim(y_sub))
                print(j)
                mfit <- cv.glmnet( x = x_sub, y = y_sub, family = 'mgaussian', alpha =j)
                #print('finish fit')
                assess <- assess.glmnet(mfit, newx = x_sub_test, newy = y_sub_test, s = 'lambda.min')
                #print('finish assess')
                result[toString(j),i] <- assess$mse[[1]]
            }
        }
    }
    print(result)
    saveRDS(result, file = 'working_data/glmnet_alpha_search_finer_0_0.1.rds')
}

#alpha search 0 - 1
if(FALSE)
{
    library(glmnet)
    dep <- readRDS('working_data/exp_all.rds')
    pred <- readRDS('working_data/pred_all.rds')
    pred <- pred[rownames(dep),]
    pred <- as.matrix(pred)
    dep <- as.matrix(dep)
    #data(QuickStartExample)
    #dep <- QuickStartExample$x
    #pred <- cbind(QuickStartExample$y, QuickStartExample$y)
    print(dim(pred))
    print(dim(dep))
    set.seed(20000)
    alpha_grid <- c(0,0.1,0.25,0.5,0.75,0.9,1)
    #alpha_grid <- c(0,0.01,0.025,0.05,0.075,0.09,0.1)

    nfolds <- 10
    train_idx <- sample.int(n = nrow(dep), size = floor(.75*nrow(dep)), replace = F)
    x_train <- dep[train_idx,]
    y_train <- pred[train_idx,]
    x_test <- dep[-train_idx,]
    y_test <- pred[-train_idx,]
    print('finish train test split')

    rowname <- c()
    for( i in alpha_grid)
    {
        rowname <- c(rowname, toString(i))
    }
    
    library(doMC)
    para <- 1
    registerDoMC(10)
    if( para == 1 )
    {
        outlist <- foreach( i = alpha_grid, .packages = c('glmnet')) %dopar%
        {
            result <- c()
            mfit <- cv.glmnet( x = x_train, y = y_train, family = 'mgaussian', alpha = i)
            assess <- assess.glmnet(mfit, newx = x_test, newy = y_test, s = 'lambda.min')
                #print('finish assess')
            result <- c(result,assess$mse[[1]])
            return(result)
        }
        outlist <- as.data.frame(do.call(rbind, outlist))
        #print(outlist)
        result <- outlist
        rownames(result) <- rowname
        colnames(result) <- c('test')
        print(result)
        saveRDS(result, file = 'working_data/glmnet_alpha_0_1_test_result.rds')
    }
}

#ploting CV result
if(FALSE)
{
    summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                          conf.interval=.95, .drop=TRUE) {
        library(plyr)

        # New version of length which can handle NA's: if na.rm==T, don't count them
        length2 <- function (x, na.rm=FALSE) {
            if (na.rm) sum(!is.na(x))
            else       length(x)
        }

        # This does the summary. For each group's data frame, return a vector with
        # N, mean, and sd
        datac <- ddply(data, groupvars, .drop=.drop,
          .fun = function(xx, col) {
            c(N    = length2(xx[[col]], na.rm=na.rm),
              mean = mean   (xx[[col]], na.rm=na.rm),
              sd   = sd     (xx[[col]], na.rm=na.rm)
            )
          },
          measurevar
        )

        # Rename the "mean" column    
        datac <- rename(datac, c("mean" = measurevar))

        datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

        # Confidence interval multiplier for standard error
        # Calculate t-statistic for confidence interval: 
        # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
        ciMult <- qt(conf.interval/2 + .5, datac$N-1)
        datac$ci <- datac$se * ciMult

        return(datac)
    }
    library(reshape2)
    result <- readRDS( file = 'working_data/glmnet_alpha_search.rds')
    test <- readRDS( file = 'working_data/glmnet_alpha_0_1_test_result.rds')
    print(result)
    #result$alpha <- rownames(result)
    result <- melt(as.matrix(result))
    result <- result[,colnames(result) != 'Var2']
    colnames(result) <- c('alpha','MSE')
    result$error <- 'train'
    print(head(result))
    colnames(test) <- 'MSE'
    test$alpha <- rownames(test)
    test$error <- 'test'
    result <- rbind(test,result)
    result$alpha <- as.double(result$alpha)
    #print(result)
    result <- summarySE( data = result, measurevar = 'MSE', groupvars = c('alpha','error'))
    print(result)
    library(ggplot2)
    pdf("figure/elastic_net_result.pdf", onefile = TRUE, width = 10, height = 9)
    graph <- ggplot(result, aes(y = MSE, x = alpha, colour = error)) +
        geom_errorbar( aes(ymin = MSE - se, ymax = MSE + se), width = 0.05)+
        geom_line()+
        geom_point(size=3, shape=21, fill="white")+
        labs(y = 'Mean squared error +- 1 SE', colour = 'data')
    graph <- graph + theme_classic(base_size = 18)+
            theme(legend.position = c(0.1,0.9))
    print(graph)
    dev.off()
}

if(FALSE)
{
    library(glmnet)
    dep <- readRDS('working_data/exp_all.rds')
    pred <- readRDS('working_data/pred_all.rds')
    pred <- pred[rownames(dep),]
    pred <- as.matrix(pred)
    dep <- as.matrix(dep)
    print(dim(pred))
    print(dim(dep))
    set.seed(20000)
    result <- readRDS(file = 'working_data/glmnet_alpha_search.rds')
    print(rowMeans(result))

    train_idx <- sample.int(n = nrow(dep), size = floor(.75*nrow(dep)), replace = F)
    x_train <- dep[train_idx,]
    y_train <- pred[train_idx,]
    x_test <- dep[-train_idx,]
    y_test <- pred[-train_idx,]
    #print('finish train test split')

    library(doMC)
    para <- 1
    registerDoMC(10)
    if( para == 1 )
    {
        mfit <- cv.glmnet( x = x_train, y = y_train, family = 'mgaussian', alpha = 0, parallel = TRUE)
        #assess <- assess.glmnet(mfit, newx = x_test, newy = y_test, s = 'lambda.min')
        assess <- assess.glmnet(mfit, newx = x_test, newy = y_test)
        print(assess$mse)
        #print(assess$mse[[1]])
        #print(mfit)
        pdf("figure/elastic_net_0_trained.pdf", onefile = TRUE)
        plt <- plot(mfit)
        print(plt)
        print(mfit$lambda.min)
        print(mfit$lambda.1se)
        dev.off()
    }
    else
    {
        print('no parallel')
    }
}

if(FALSE)
{
    dep <- readRDS('working_data/exp_all.rds')
    pred <- readRDS('working_data/pred_all.rds')
    pred <- pred[rownames(dep),]
    pred <- as.matrix(pred)
    dep <- as.matrix(dep)
    print(dim(pred))
    print(dim(dep))
    set.seed(20000)
    result <- readRDS(file = 'working_data/glmnet_alpha_search.rds')
    print(rowMeans(result))

    train_idx <- sample.int(n = nrow(dep), size = floor(.75*nrow(dep)), replace = F)
    x_train <- dep[train_idx,]
    y_train <- pred[train_idx,]
    x_test <- dep[-train_idx,]
    y_test <- pred[-train_idx,]
    #print('finish train test split')
    print(head(y_test))
    print(colMeans(y_test))
    diff <- y_test - colMeans(y_test)
    print(head(diff))
    diff <- diff**2
    print(head(diff))
    diff <- sum(colSums(diff)) / 2 / dim(y_test)[1]
    print(diff)
}

if( FALSE )
{
    library('randomForestSRC')
    #data(nutrigenomic)
    dep <- readRDS('working_data/exp_all.rds')
    pred <- readRDS('working_data/pred_all.rds')
    pred <- pred[rownames(dep),]
    #pred <- nutrigenomic$lipids
    #dep <- data.frame(nutrigenomic$genes, diet = nutrigenomic$diet,genotype = nutrigenomic$genotype)
    pred <- as.matrix(pred)
    dep <- as.matrix(dep)
    print(dim(pred))
    print(dim(dep))
    set.seed(20000)
    train_idx <- sample.int(n = nrow(dep), size = floor(.75*nrow(dep)), replace = F)
    x_train <- dep[train_idx,]
    y_train <- pred[train_idx,]
    x_test <- dep[-train_idx,]
    y_test <- pred[-train_idx,]
    mv.obj <- rfsrc(get.mv.formula(colnames(y_train)), data.frame(y_train, x_train), importance = TRUE, ntree = 100)
    print(mv.obj)
    print('finish training')
    yhat <- predict(mv.obj, newdata = data.frame(y_test,x_test))
    print(yhat)
    saveRDS(mv.obj, file = 'working_data/rf_100_train_model_importance.rds')
}

if( FALSE )
{
    library('randomForestSRC')
    library('ggplot2')
    #data(nutrigenomic)
    dep <- readRDS('working_data/exp_all.rds')
    pred <- readRDS('working_data/pred_all.rds')
    pred <- pred[rownames(dep),]
    #pred <- nutrigenomic$lipids
    #dep <- data.frame(nutrigenomic$genes, diet = nutrigenomic$diet,genotype = nutrigenomic$genotype)
    pred <- as.matrix(pred)
    dep <- as.matrix(dep)
    print(dim(pred))
    print(dim(dep))
    set.seed(20000)
    train_idx <- sample.int(n = nrow(dep), size = floor(.75*nrow(dep)), replace = F)
    x_train <- dep[train_idx,]
    y_train <- pred[train_idx,]
    x_test <- dep[-train_idx,]
    y_test <- pred[-train_idx,]
    mv.obj <- readRDS(file = 'working_data/rf_1000_train_model_importance.rds')
    #print(mv.obj)
    #yhat <- predict(mv.obj, newdata = data.frame(pred,dep))
    yhat <- predict(mv.obj, newdata = data.frame(x_test))
    #print(str(yhat))
    print(yhat)
}






if( FALSE )
{
    library('randomForestSRC')
    library('ggplot2')
    #data(nutrigenomic)
    dep <- readRDS('working_data/exp_all.rds')
    pred <- readRDS('working_data/pred_all.rds')
    pred <- pred[rownames(dep),]
    #pred <- nutrigenomic$lipids
    #dep <- data.frame(nutrigenomic$genes, diet = nutrigenomic$diet,genotype = nutrigenomic$genotype)
    pred <- as.matrix(pred)
    dep <- as.matrix(dep)
    print(dim(pred))
    print(dim(dep))
    set.seed(20000)
    train_idx <- sample.int(n = nrow(dep), size = floor(.75*nrow(dep)), replace = F)
    x_train <- dep[train_idx,]
    y_train <- pred[train_idx,]
    x_test <- dep[-train_idx,]
    y_test <- pred[-train_idx,]
    mv.obj <- readRDS(file = 'working_data/rf_1000_train_model_importance.rds')
    #print(mv.obj)
    #yhat <- predict(mv.obj, newdata = data.frame(pred,dep))
    yhat <- predict(mv.obj, newdata = data.frame(y_test,x_test))
    #print(str(yhat))
    #print(yhat)
    print(cor(yhat$yvar$GSH,yhat$regrOutput$GSH$predicted))
    print(cor(yhat$yvar$proliferate,yhat$regrOutput$proliferate$predicted))
    df <- data.frame(matrix(ncol = 2, nrow = length(yhat$yvar$GSH)))
    colnames(df) <- c('pred','true')
    df[,'pred'] <- yhat$regrOutput$GSH$predicted
    df[,'true'] <- yhat$yvar$GSH
    pdf("figure/rf_1000_train_model_scatter_plot.pdf", onefile = TRUE, width = 9, height = 9)
    graph <- ggplot(df, aes(y = true, x = pred)) +
        geom_smooth(method='lm')+
        geom_point(alpha = 0.6, size = 4) +
        labs(title = 'GSH')
    graph <- graph + theme_classic(base_size = 18)
    print(graph)
    df[,'pred'] <- yhat$regrOutput$proliferate$predicted
    df[,'true'] <- yhat$yvar$proliferate
    graph <- ggplot(df, aes(y = true, x = pred)) +
        geom_smooth(method='lm')+
        geom_point(alpha = 0.6, size = 4) +
        labs(title = 'proliferate')
    graph <- graph + theme_classic(base_size = 18)
    print(graph)
    dev.off()
}

if( FALSE )
{
    library('randomForestSRC')
    #data(nutrigenomic)
    options(rf.cores = 10)
    dep <- readRDS('working_data/exp_all.rds')
    pred <- readRDS('working_data/pred_all.rds')
    pred <- pred[rownames(dep),]
    #pred <- nutrigenomic$lipids
    #dep <- data.frame(nutrigenomic$genes, diet = nutrigenomic$diet,genotype = nutrigenomic$genotype)
    pred <- as.matrix(pred)
    dep <- as.matrix(dep)
    print(dim(pred))
    print(dim(dep))
    set.seed(20000)
    x_train <- dep[,]
    y_train <- pred[,]
    mv.obj <- rfsrc(get.mv.formula(colnames(y_train)), data.frame(y_train, x_train), importance = TRUE, ntree = 1000)
    print(mv.obj)
}

if( FALSE )
{
    library('randomForestSRC')
    dep <- readRDS('working_data/exp_all.rds')
    pred <- readRDS('working_data/pred_all.rds')
    pred <- pred[rownames(dep),]
    pred <- as.matrix(pred)
    dep <- as.matrix(dep)
    print(dim(pred))
    print(dim(dep))
    set.seed(20000)
    train_idx <- sample.int(n = nrow(dep), size = floor(.75*nrow(dep)), replace = F)
    x_train <- dep[train_idx,]
    y_train <- pred[train_idx,]
    x_test <- dep[-train_idx,]
    y_test <- pred[-train_idx,]
    #print(dim(x_train))
    #print(head(x_train))
    o <- tune(get.mv.formula(colnames(y_train)), data.frame(y_train, x_train), ntreeTry = 500, trace = TRUE)
    print(str(o))
    saveRDS( o, file = 'working_data/randomforest_tune_test.rds')
}

if( FALSE )
{
    library('randomForestSRC')
    o <- readRDS( file = 'working_data/randomforest_tune.rds')
    print(o)
    print(str(o))
}



if( FALSE )
{
    library('randomForestSRC')
    o <- readRDS( file = 'working_data/randomforest_tune_test.rds')
    #print(o)
    pdf("figure/randomforest_tune.pdf", onefile = TRUE)
    if (library("akima", logical.return = TRUE)) {

      ## nice little wrapper for plotting results
      #plot.tune <- function(o, linear = TRUE) {
      plot.tune <- function(o, linear = TRUE) {
          print('in function')
        x <- o$results[,1]
        print(x)
        y <- o$results[,2]
        print(y)
        z <- o$results[,3]
        print(z)
        so <- interp(x=x, y=y, z=z, linear = linear)
        idx <- which.min(z)
        x0 <- x[idx]
        y0 <- y[idx]
        filled.contour(x = so$x,
                       y = so$y,
                       z = so$z,
                       xlim = range(so$x, finite = TRUE) + c(-2, 2),
                       ylim = range(so$y, finite = TRUE) + c(-2, 2),
                       color.palette =
                         colorRampPalette(c("yellow", "red")),
                       xlab = "nodesize",
                       ylab = "mtry",
                       main = "error rate for nodesize and mtry",
                       key.title = title(main = "OOB error", cex.main = 1),
                       plot.axes = {axis(1);axis(2);points(x0,y0,pch="x",cex=1,font=2);
                                points(x,y,pch=16,cex=.25)})
    }  

    ## plot the surface
    plot.tune(o)

    }
    dev.off()
}


if( FALSE )
{
    library('randomForestSRC')
    library('dplyr')
    exp <- readRDS('working_data/exp_all.rds')
    dep <- readRDS('working_data/dep_all.rds')
    pred <- readRDS('working_data/pred_all.rds')
    dep_col <- c()
    for( i in colnames(dep) )
    {
        dep_col <- c( dep_col, paste0(i,'_dep'))
    }
    colnames(dep) <- dep_col
    exp_col <- c()
    for( i in colnames(exp) )
    {
        exp_col <- c( exp_col, paste0(i,'_exp'))
    }
    colnames(exp) <- exp_col
    pred_col <- c()
    for( i in rownames(pred) )
    {
        pred_col <- c( pred_col, gsub('-','.',i))
    }
    rownames(pred) <- pred_col

    dep <- t(dep)
    exp <- t(exp)
    #print(head(dep))
    #print(head(exp))
    dep <- data.frame(dep)
    exp <- data.frame(exp)

    dep <- bind_rows( dep, exp)
    dep <- t(dep)
    dep <- data.frame(dep)
    dep <- dep[complete.cases(dep),]
    print(head(pred))
    
    exp_col <- grep('_exp', colnames(dep))

    pred <- pred[rownames(dep),]
    pred <- pred[complete.cases(pred),]
    dep <- dep[rownames(pred),exp_col]
    pred <- as.matrix(pred)
    dep <- as.matrix(dep)
    print(dim(pred))
    print(dim(dep))
    
    set.seed(20000)
    x_train <- dep[,]
    y_train <- pred[,]
    mv.obj <- rfsrc(get.mv.formula(colnames(y_train)), data.frame(y_train, x_train), importance = TRUE, ntree = 1000)
    print(mv.obj)
}


if(FALSE)
{
    library(glmnet)
    dep <- readRDS('working_data/exp_all.rds')
    pred <- readRDS('working_data/pred_all.rds')
    pred <- pred[rownames(dep),]
    pred <- as.matrix(pred)
    dep <- as.matrix(dep)
    print(dim(pred))
    print(dim(dep))

    set.seed(20000)
    #mfit <- glmnet( x = dep, y = pred, family = 'mgaussian')
    mfit <- cv.glmnet( x = dep, y = pred, family = 'mgaussian', alpha = 0.05)
    print(mfit)
    saveRDS(mfit, file = 'working_data/elastic_net_alpha_0.05_test.rds')
    #s <- 0.1199737
    #result <- assess.glmnet(mfit, newx = dep, newy = pred, s = 0.1199737)
    #result <- predict(mfit, newx = dep, s = 0.1199737)
    #print(result)
    #print(str(result))
    #print(mfit)
    #pdf("figure/elastic_net_test.pdf", onefile = TRUE)
    #plt <- plot(mfit)
	#print(plt)
    #print(mfit$lambda.min)
    #print(mfit$lambda.1se)
	#dev.off()
}

if(FALSE)
{
    sc <- readRDS(file = 'working_data/persister_rna/all_CL_features.rds')
    print(str(sc))
}

if(FALSE)
{
    drug <- c('trametinib', 'dabrafenib', 'brd3379', 'navitoclax', 'everolimus',
       'jq1', 'azd5591', 'afatinib', 'prexasertib', 'taselisib', 'gemcitabine',
       'bortezomib', 'idasanutlin')
    for( i in drug)
    {
        dep <- readRDS(paste0('working_data/mix_seq_result/',i,'persister_rna.rds'))
        mat <- read.table( paste0('working_data/mix_seq_result/',i,'persister_rna_index.csv'),sep=',', header = TRUE);
        rownames(dep) <- mat[,1]
        saveRDS(dep, file = paste0('./working_data/mix_seq_result/',i,'persister_rna_all.rds'))
        
        dep <- readRDS(paste0('working_data/mix_seq_result/',i,'persister_rna_untreated.rds'))
        mat <- read.table( paste0('working_data/mix_seq_result/',i,'persister_rna_untreated_index.csv'),sep=',', header = TRUE);
        rownames(dep) <- mat[,1]
        saveRDS(dep, file = paste0('./working_data/mix_seq_result/',i,'persister_rna_untreated_all.rds'))
        
        dep <- readRDS(paste0('working_data/mix_seq_result/',i,'persister_rna_untreated_ans.rds'))
        mat <- read.table( paste0('working_data/mix_seq_result/',i,'persister_rna_untreated_ans_index.csv'),sep=',', header = TRUE);
        rownames(dep) <- mat[,1]
        saveRDS(dep, file = paste0('./working_data/mix_seq_result/',i,'persister_rna_untreated_ans_all.rds'))
    }
}

if( FALSE )
{
    library('glmnet')
    library('ggplot2')
    library('ggrepel')
    pred <- readRDS('working_data/pred_all.rds')
    mfit <- readRDS(file = 'working_data/elastic_net_alpha_1_test.rds')

    drug <- c('trametinib', 'dabrafenib', 'brd3379', 'navitoclax', 'everolimus',
       'jq1', 'azd5591', 'afatinib', 'prexasertib', 'taselisib', 'gemcitabine',
       'bortezomib', 'idasanutlin')
    for( i in drug)
    {
        dep <- readRDS(file = paste0('./working_data/mix_seq_result/',i,'persister_rna_all.rds'))
        print(dim(dep))
        dep <- as.matrix(dep)
        yhat <- predict(mfit, newx = dep, s = 'lambda.min')
        yhat <- yhat[,,1]
        df_persister <- data.frame(matrix(ncol = 2, nrow = length(rownames(dep))))
        colnames(df_persister) <- c('GSH','proliferate')
        df_persister[,'GSH'] <- yhat[,1]
        df_persister[,'proliferate'] <- yhat[,2]
        df_persister[,'type'] <- 'persister'
        df_persister[,'name'] <- rownames(dep)
        rownames(df_persister) <- rownames(dep)

        dep <- readRDS(file = paste0('./working_data/mix_seq_result/',i,'persister_rna_untreated_all.rds'))
        dep <- as.matrix(dep)
        print(dim(dep))
        yhat <- predict(mfit, newx = dep, s = 'lambda.min')
        print('success')
        yhat <- yhat[,,1]
        df_estimate <- data.frame(matrix(ncol = 2, nrow = length(rownames(dep))))
        colnames(df_estimate) <- c('GSH','proliferate')
        df_estimate[,'GSH'] <- yhat[,1]
        df_estimate[,'proliferate'] <- yhat[,2]
        df_estimate[,'type'] <- 'estimate_true'
        df_estimate[,'name'] <- rownames(dep)
        rownames(df_estimate) <- rownames(dep)
        
        dep <- readRDS(file = paste0('./working_data/mix_seq_result/',i,'persister_rna_untreated_ans_all.rds'))
        dep[,'name'] <- rownames(dep)
        dep[,'type'] <- 'true_value'

        pred[,'name'] <- ''
        pred[,'type'] <- 'depmap'
        print(head(df_persister))
        print(head(df_estimate))
        print(head(dep))
        colnames(dep) <- c('GSH','proliferate','name','type')

        top <- rbind(df_estimate, dep, df_persister)
        top[,'alpha'] <- 0.9
        pred[,'alpha'] <- 0.9
        tmp_pred <- pred
        print(head(pred))
        print(head(top))
        pred <- rbind(pred,top)

        pdf(paste0("figure/mix_seq_result/",i,"_elastic_model_predict_test_2_d.pdf"), onefile = TRUE, width = 9, height = 9)
        graph <- ggplot( pred, aes(y = proliferate, x = GSH, color = type))+
            geom_point( aes(alpha = alpha), size = 3)+
            scale_color_manual( values = c('depmap' = 'azure2', 'true_value' = 'darkgoldenrod2','estimate_true' = 'cornflowerblue', 'persister' = 'coral' ) )+
            scale_x_continuous(limits = c(-2.5,2.5))+
            scale_y_continuous(limits = c(-2.1,4.1))+
            guides( alpha='none')+
            geom_text_repel( aes(label = name), max.overlaps = 100)+
            labs(title = i)
        graph <- graph + theme_classic(base_size = 18)+
                theme(legend.position = c(0.2,0.9))
        print(graph)

        #create line graph
        line <- readRDS(file = paste0('working_data/mix_seq_result/',i,'persister_line.rds'))
        rownames(line) <- line[,'sample_name']
        print(head(line))
        df_line <- data.frame(matrix(ncol = 4, nrow = length(rownames(df_persister))))
        colnames(df_line) <- c('x1','y1','x2','y2')
        rownames(df_line) <- rownames(line)
        for( j in rownames(line) )
        {
            df_line[j,'x2'] <- df_persister[line[j,'sample_name'],'GSH']
            df_line[j,'y2'] <- df_persister[line[j,'sample_name'],'proliferate']
            df_line[j,'x1'] <- df_estimate[line[j,'cell'],'GSH']
            df_line[j,'y1'] <- df_estimate[line[j,'cell'],'proliferate']
        }

        top <- rbind(df_estimate, df_persister)
        top[,'alpha'] <- 0.9
        pred <- rbind(tmp_pred,top)
        graph <- ggplot( pred, aes(y = proliferate, x = GSH))+
            geom_point( aes(alpha = alpha, color = type), size = 3)+
            scale_color_manual( values = c('depmap' = 'azure2', 'true_value' = 'darkgoldenrod2','estimate_true' = 'cornflowerblue', 'persister' = 'coral' ) )+
            scale_x_continuous(limits = c(-2.5,2.5))+
            scale_y_continuous(limits = c(-2.1,4.1))+
            guides( alpha='none')+
            labs(title = i)
        graph <- graph + geom_segment( aes( x = x1, y = y1, xend = x2, yend = y2),data = df_line, linetype = 'dashed')+
                theme_classic(base_size = 18)+
                theme(legend.position = c(0.2,0.9))
        print(graph)
        graph <- ggplot( pred, aes(y = proliferate, x = GSH))+
            geom_point( aes(alpha = alpha, color = type), size = 3)+
            scale_color_manual( values = c('depmap' = 'azure2', 'true_value' ='azure2','estimate_true' = 'azure2', 'persister' = 'azure2' ) )+
            scale_x_continuous(limits = c(-2.5,2.5))+
            scale_y_continuous(limits = c(-2.1,4.1))+
            guides( alpha='none')+
            labs(title = i)
        graph <- graph + geom_segment( aes( x = x1, y = y1, xend = x2, yend = y2),data = df_line, arrow = arrow(length = unit(0.1,"cm")))+
                theme_classic(base_size = 18)+
                theme(legend.position = c(0.2,0.9))
        print(graph)
        dev.off()
    }
}

if( FALSE )
{
    library('randomForestSRC')
    library('ggplot2')
    library('ggrepel')
    pred <- readRDS('working_data/pred_all.rds')
    mv.obj <- readRDS(file = 'working_data/rf_1000_train_model_importance.rds')

    drug <- c('trametinib', 'dabrafenib', 'brd3379', 'navitoclax', 'everolimus',
       'jq1', 'azd5591', 'afatinib', 'prexasertib', 'taselisib', 'gemcitabine',
       'bortezomib', 'idasanutlin')
    for( i in drug)
    {
        dep <- readRDS(file = paste0('./working_data/mix_seq_result/',i,'persister_rna_all.rds'))
        yhat <- predict(mv.obj, newdata = data.frame(dep))
        df_persister <- data.frame(matrix(ncol = 2, nrow = length(rownames(dep))))
        colnames(df_persister) <- c('GSH','proliferate')
        df_persister[,'GSH'] <- yhat$regrOutput$GSH$predicted
        df_persister[,'proliferate'] <- yhat$regrOutput$proliferate$predicted
        df_persister[,'type'] <- 'persister'
        df_persister[,'name'] <- rownames(dep)
        rownames(df_persister) <- rownames(dep)

        dep <- readRDS(file = paste0('./working_data/mix_seq_result/',i,'persister_rna_untreated_all.rds'))
        yhat <- predict(mv.obj, newdata = data.frame(dep))
        df_estimate <- data.frame(matrix(ncol = 2, nrow = length(rownames(dep))))
        colnames(df_estimate) <- c('GSH','proliferate')
        df_estimate[,'GSH'] <- yhat$regrOutput$GSH$predicted
        df_estimate[,'proliferate'] <- yhat$regrOutput$proliferate$predicted
        df_estimate[,'type'] <- 'estimate_true'
        df_estimate[,'name'] <- rownames(dep)
        rownames(df_estimate) <- rownames(dep)
        
        dep <- readRDS(file = paste0('./working_data/mix_seq_result/',i,'persister_rna_untreated_ans_all.rds'))
        dep[,'name'] <- rownames(dep)
        dep[,'type'] <- 'true_value'

        pred[,'name'] <- ''
        pred[,'type'] <- 'depmap'
        colnames(dep) <- c('GSH','proliferate','name','type')
        print(head(pred))
        print(head(dep))
        print(head(df_persister))
        print(head(df_estimate))

        top <- rbind(df_estimate, dep, df_persister)
        top[,'alpha'] <- 0.9
        pred[,'alpha'] <- 0.9
        tmp_pred <- pred
        pred <- rbind(pred,top)
        
        pdf(paste0("figure/mix_seq_result/",i,"_rf_1000_test_2_d.pdf"), onefile = TRUE, width = 9, height = 9)
        graph <- ggplot( pred, aes(y = proliferate, x = GSH, color = type))+
            geom_point( aes(alpha = alpha), size = 3)+
            scale_color_manual( values = c('depmap' = 'azure2', 'true_value' = 'darkgoldenrod2','estimate_true' = 'cornflowerblue', 'persister' = 'coral' ) )+
            scale_x_continuous(limits = c(-2.5,2.5))+
            scale_y_continuous(limits = c(-2.1,4.1))+
            guides( alpha='none')+
            geom_text_repel( aes(label = name), max.overlaps = 100)+
            labs(title = i)
        graph <- graph + theme_classic(base_size = 18)+
                theme(legend.position = c(0.2,0.9))
        print(graph)

        #create line graph
        line <- readRDS(file = paste0('working_data/mix_seq_result/',i,'persister_line.rds'))
        rownames(line) <- line[,'sample_name']
        print(head(line))
        df_line <- data.frame(matrix(ncol = 4, nrow = length(rownames(df_persister))))
        colnames(df_line) <- c('x1','y1','x2','y2')
        rownames(df_line) <- rownames(line)
        for( j in rownames(line) )
        {
            df_line[j,'x2'] <- df_persister[line[j,'sample_name'],'GSH']
            df_line[j,'y2'] <- df_persister[line[j,'sample_name'],'proliferate']
            df_line[j,'x1'] <- df_estimate[line[j,'cell'],'GSH']
            df_line[j,'y1'] <- df_estimate[line[j,'cell'],'proliferate']
        }

        top <- rbind(df_estimate, df_persister)
        top[,'alpha'] <- 0.9
        pred <- rbind(tmp_pred,top)
        graph <- ggplot( pred, aes(y = proliferate, x = GSH))+
            geom_point( aes(alpha = alpha, color = type), size = 3)+
            scale_color_manual( values = c('depmap' = 'azure2', 'true_value' = 'darkgoldenrod2','estimate_true' = 'cornflowerblue', 'persister' = 'coral' ) )+
            scale_x_continuous(limits = c(-2.5,2.5))+
            scale_y_continuous(limits = c(-2.1,4.1))+
            guides( alpha='none')+
            labs(title = i)
        graph <- graph + geom_segment( aes( x = x1, y = y1, xend = x2, yend = y2),data = df_line, linetype = 'dashed')+
                theme_classic(base_size = 18)+
                theme(legend.position = c(0.2,0.9))
        print(graph)
        graph <- ggplot( pred, aes(y = proliferate, x = GSH))+
            geom_point( aes(alpha = alpha, color = type), size = 3)+
            scale_color_manual( values = c('depmap' = 'azure2', 'true_value' ='azure2','estimate_true' = 'azure2', 'persister' = 'azure2' ) )+
            scale_x_continuous(limits = c(-2.5,2.5))+
            scale_y_continuous(limits = c(-2.1,4.1))+
            guides( alpha='none')+
            labs(title = i)
        graph <- graph + geom_segment( aes( x = x1, y = y1, xend = x2, yend = y2),data = df_line, arrow = arrow(length = unit(0.1,"cm")))+
                theme_classic(base_size = 18)+
                theme(legend.position = c(0.2,0.9))
        print(graph)
        dev.off()
    }
}

if( FALSE )
{
    dep <- readRDS('working_data/persister_rna.rds')
    mat <- read.table( paste('working_data/persister_rna_index.csv'),sep=',', header = TRUE);
    rownames(dep) <- mat[,1]
    saveRDS(dep, file = './working_data/persister_rna_all.rds')
    
    dep <- readRDS('working_data/persister_rna_untreated.rds')
    mat <- read.table( paste('working_data/persister_rna_untreated_index.csv'),sep=',', header = TRUE);
    rownames(dep) <- mat[,1]
    saveRDS(dep, file = './working_data/persister_rna_untreated_all.rds')
    
    dep <- readRDS('working_data/persister_rna_untreated_ans.rds')
    mat <- read.table( paste('working_data/persister_rna_untreated_ans_index.csv'),sep=',', header = TRUE);
    rownames(dep) <- mat[,1]
    saveRDS(dep, file = './working_data/persister_rna_untreated_ans_all.rds')
}
if( FALSE )
{
    library('glmnet')
    library('ggplot2')
    library('ggrepel')
    pred <- readRDS('working_data/pred_all.rds')
    mfit <- readRDS(file = 'working_data/elastic_net_alpha_1_test.rds')

    dep <- readRDS(file = './working_data/persister_rna_all.rds')
    dep <- as.matrix(dep)
    yhat <- predict(mfit, newx = dep, s = 'lambda.min')
    yhat <- yhat[,,1]
    df_persister <- data.frame(matrix(ncol = 2, nrow = length(rownames(dep))))
    colnames(df_persister) <- c('GSH','proliferate')
    df_persister[,'GSH'] <- yhat[,1]
    df_persister[,'proliferate'] <- yhat[,2]
    df_persister[,'type'] <- 'persister'
    df_persister[,'name'] <- rownames(dep)
    rownames(df_persister) <- rownames(dep)

    dep <- readRDS(file = './working_data/persister_rna_untreated_all.rds')
    dep <- as.matrix(dep)
    yhat <- predict(mfit, newx = dep, s = 'lambda.min')
    yhat <- yhat[,,1]
    df_estimate <- data.frame(matrix(ncol = 2, nrow = length(rownames(dep))))
    colnames(df_estimate) <- c('GSH','proliferate')
    df_estimate[,'GSH'] <- yhat[,1]
    df_estimate[,'proliferate'] <- yhat[,2]
    df_estimate[,'type'] <- 'estimate_true'
    df_estimate[,'name'] <- rownames(dep)
    rownames(df_estimate) <- rownames(dep)
    
    dep <- readRDS(file = './working_data/persister_rna_untreated_ans_all.rds')
    dep[,'name'] <- rownames(dep)
    dep[,'type'] <- 'true_value'

    pred[,'name'] <- ''
    pred[,'type'] <- 'depmap'
    #print(head(pred))
    #print(head(dep))
    print(df_persister)
    #print(head(df_estimate))

    top <- rbind(df_estimate, dep, df_persister)
    top[,'alpha'] <- 0.9
    pred[,'alpha'] <- 0.9
    tmp_pred <- pred
    pred <- rbind(pred,top)

    pdf("figure/elastic_model_predict_test_2.pdf", onefile = TRUE, width = 9, height = 9)
    graph <- ggplot( pred, aes(y = proliferate, x = GSH, color = type))+
        geom_point( aes(alpha = alpha), size = 3)+
        scale_color_manual( values = c('depmap' = 'azure2', 'true_value' = 'darkgoldenrod2','estimate_true' = 'cornflowerblue', 'persister' = 'coral' ) )+
        scale_x_continuous(limits = c(-2.5,2.5))+
        scale_y_continuous(limits = c(-2.1,4.1))+
        guides( alpha='none')+
        geom_text_repel( aes(label = name), max.overlaps = 100)+
        labs()
    graph <- graph + theme_classic(base_size = 18)+
            theme(legend.position = c(0.2,0.9))
    print(graph)

    #create line graph
    line <- readRDS(file = 'working_data/persister_line.rds')
    rownames(line) <- line[,'sample_name']
    print(head(line))
    df_line <- data.frame(matrix(ncol = 4, nrow = length(rownames(df_persister))))
    colnames(df_line) <- c('x1','y1','x2','y2')
    rownames(df_line) <- rownames(line)
    for( i in rownames(line) )
    {
        df_line[i,'x2'] <- df_persister[line[i,'sample_name'],'GSH']
        df_line[i,'y2'] <- df_persister[line[i,'sample_name'],'proliferate']
        df_line[i,'x1'] <- df_estimate[line[i,'cell'],'GSH']
        df_line[i,'y1'] <- df_estimate[line[i,'cell'],'proliferate']
    }

    top <- rbind(df_estimate, df_persister)
    top[,'alpha'] <- 0.9
    pred <- rbind(tmp_pred,top)
    graph <- ggplot( pred, aes(y = proliferate, x = GSH))+
        geom_point( aes(alpha = alpha, color = type), size = 3)+
        scale_color_manual( values = c('depmap' = 'azure2', 'true_value' = 'darkgoldenrod2','estimate_true' = 'cornflowerblue', 'persister' = 'coral' ) )+
        scale_x_continuous(limits = c(-2.5,2.5))+
        scale_y_continuous(limits = c(-2.1,4.1))+
        guides( alpha='none')+
        labs()
    graph <- graph + geom_segment( aes( x = x1, y = y1, xend = x2, yend = y2),data = df_line, linetype = 'dashed')+
            theme_classic(base_size = 18)+
            theme(legend.position = c(0.2,0.9))
    print(graph)
    dev.off()
}

if( FALSE )
{
    library('randomForestSRC')
    library('ggplot2')
    library('ggrepel')
    pred <- readRDS('working_data/pred_all.rds')
    mv.obj <- readRDS(file = 'working_data/rf_1000_train_model_importance.rds')

    dep <- readRDS(file = './working_data/persister_rna_all.rds')
    yhat <- predict(mv.obj, newdata = data.frame(dep))
    df_persister <- data.frame(matrix(ncol = 2, nrow = length(rownames(dep))))
    colnames(df_persister) <- c('GSH','proliferate')
    df_persister[,'GSH'] <- yhat$regrOutput$GSH$predicted
    df_persister[,'proliferate'] <- yhat$regrOutput$proliferate$predicted
    df_persister[,'type'] <- 'persister'
    df_persister[,'name'] <- rownames(dep)
    rownames(df_persister) <- rownames(dep)

    dep <- readRDS(file = './working_data/persister_rna_untreated_all.rds')
    yhat <- predict(mv.obj, newdata = data.frame(dep))
    df_estimate <- data.frame(matrix(ncol = 2, nrow = length(rownames(dep))))
    colnames(df_estimate) <- c('GSH','proliferate')
    df_estimate[,'GSH'] <- yhat$regrOutput$GSH$predicted
    df_estimate[,'proliferate'] <- yhat$regrOutput$proliferate$predicted
    df_estimate[,'type'] <- 'estimate_true'
    df_estimate[,'name'] <- rownames(dep)
    rownames(df_estimate) <- rownames(dep)
    
    dep <- readRDS(file = './working_data/persister_rna_untreated_ans_all.rds')
    dep[,'name'] <- rownames(dep)
    dep[,'type'] <- 'true_value'

    pred[,'name'] <- ''
    pred[,'type'] <- 'depmap'
    print(head(pred))
    print(head(dep))
    print(head(df_persister))
    print(head(df_estimate))

    top <- rbind(df_estimate, dep, df_persister)
    top[,'alpha'] <- 0.9
    pred[,'alpha'] <- 0.9
    tmp_pred <- pred
    pred <- rbind(pred,top)

    pdf("figure/rf_1000_train_model_predict.pdf", onefile = TRUE, width = 9, height = 9)
    graph <- ggplot( pred, aes(y = proliferate, x = GSH, color = type))+
        geom_point( aes(alpha = alpha), size = 3)+
        scale_color_manual( values = c('depmap' = 'azure2', 'true_value' = 'darkgoldenrod2','estimate_true' = 'cornflowerblue', 'persister' = 'coral' ) )+
        #scale_x_continuous(limits = c(4.5,7))+
        #scale_y_continuous(limits = c(0,0.07))+
        guides( alpha='none')+
        scale_x_continuous(limits = c(-2.5,2.5))+
        scale_y_continuous(limits = c(-2.1,4.1))+
        geom_text_repel( aes(label = name), max.overlaps = 1000)+
        labs()
    graph <- graph + theme_classic(base_size = 18)+
            theme(legend.position = c(0.2,0.9))
    print(graph)
    #create line graph
    line <- readRDS(file = 'working_data/persister_line.rds')
    rownames(line) <- line[,'sample_name']
    print(head(line))
    df_line <- data.frame(matrix(ncol = 4, nrow = length(rownames(df_persister))))
    colnames(df_line) <- c('x1','y1','x2','y2')
    rownames(df_line) <- rownames(line)
    for( i in rownames(line) )
    {
        df_line[i,'x2'] <- df_persister[line[i,'sample_name'],'GSH']
        df_line[i,'y2'] <- df_persister[line[i,'sample_name'],'proliferate']
        df_line[i,'x1'] <- df_estimate[line[i,'cell'],'GSH']
        df_line[i,'y1'] <- df_estimate[line[i,'cell'],'proliferate']
    }

    top <- rbind(df_estimate, df_persister)
    top[,'alpha'] <- 0.9
    pred <- rbind(tmp_pred,top)
    graph <- ggplot( pred, aes(y = proliferate, x = GSH))+
        geom_point( aes(alpha = alpha, color = type), size = 3)+
        scale_color_manual( values = c('depmap' = 'azure2', 'true_value' = 'darkgoldenrod2','estimate_true' = 'cornflowerblue', 'persister' = 'coral' ) )+
        scale_x_continuous(limits = c(-2.5,2.5))+
        scale_y_continuous(limits = c(-2.1,4.1))+
        guides( alpha='none')+
        labs()
    graph <- graph + geom_segment( aes( x = x1, y = y1, xend = x2, yend = y2),data = df_line, linetype = 'dashed')+
            theme_classic(base_size = 18)+
            theme(legend.position = c(0.2,0.9))
    print(graph)
    dev.off()
}

#study 3 noiseq
if(FALSE)
{
    library(NOISeq)
    cell <- c('HCC4006','NCIH1650','NCIH1975','PC9')
    control <- c('DMSO1day_0','DMSO7day_0','DMSO7day_0','DMSO2day_0')
    treat <- c('AZD1day_0','AZD7day_0','AZD7day_0','AZD2day_0')
    for( i in seq(1,4,by = 1) )
    {
        print(cell[i])
        mat_all <- read.csv( file = paste0('data/persister_rna/GSE165019_',cell[i],'.txt'), sep = '\t');
        mat_all <- mat_all[!duplicated(mat_all[,1]),]
        rownames(mat_all) <- mat_all[,1]
        print(head(mat_all))
        mat <- mat_all[,c(control[i],treat[i])]
        print(head(mat))
        mat_fact = data.frame( condition = c(control[i],treat[i]))
        mydata <- readData( data = mat, factors = mat_fact)
        myresult <- noiseq(mydata, factor = 'condition', k = NULL, norm = 'n', pnr = 0.2, nss = 5, v = 0.02, lc = 1, replicates = 'no')
        saveRDS(myresult, file = paste('working_data/persister_rna/study_3_',cell[i],'.rds',sep=''))
        write.table( myresult@results[[1]],file = paste('working_data/persister_rna/study_3_',cell[i],'.csv',sep=''),quote = F,row.names = T,col.names = T,sep = ',')
	}
}

if(FALSE)
{
    library(NOISeq)
    cell <- c('HCC4006','NCIH1650','NCIH1975','PC9')
    #control <- c('DMSO1day_0','DMSO7day_0','DMSO7day_0','DMSO2day_0')
    #treat <- c('AZD1day_0','AZD7day_0','AZD7day_0','AZD2day_0')
    for( i in seq(1,4,by = 1) )
    {
        print(cell[i])
        mat_all <- read.csv( file = paste0('data/persister_rna/GSE165019_',cell[i],'.txt'), sep = '\t');
        mat_all <- mat_all[!duplicated(mat_all[,1]),]
        rownames(mat_all) <- mat_all[,1]
        #print(head(mat_all))
        cond2 <- colnames(mat_all)[grepl('AZD', colnames(mat_all))]
        for( j in cond2)
        {
            print(j)
            mat <- mat_all[,c('DMSO4hr_0',j)]
            #print(head(mat))
            mat_fact = data.frame( condition = c('DMSO4hr_0',j))
            mydata <- readData( data = mat, factors = mat_fact)
            myresult <- noiseq( mydata, factor = 'condition', k = NULL, norm = 'n', pnr = 0.2, nss = 5, v = 0.02, lc = 1, replicates = 'no')
            saveRDS(myresult, file = paste('working_data/persister_rna/study3/study_3_',cell[i],'_',j,'.rds',sep=''))
            write.table( myresult@results[[1]],file = paste('working_data/persister_rna/study3/study_3_',cell[i],'_',j,'.csv',sep=''),quote = F,row.names = T,col.names = T,sep = ',')
        }
	}
}

#study 5 noiseq
if(FALSE)
{
    mat_all <- readRDS( file = 'working_data/persister_rna/GSE189625.rds')
    print(head(mat_all))
    rownames(mat_all) <- mat_all[,'gene']
    mat_all <- mat_all[,colnames(mat_all) != 'gene']
    print(head(mat_all))
    myfactors = data.frame( cell = c('PC9_WT', 'PC9_WT', 'PC9_WT', 'PC9_P', 'PC9_P', 'PC9_P','PC9_E', 'PC9_E', 'PC9_E', 'PC9_AS', 'PC9_AS', 'PC9_AS','HT29_WT', 'HT29_WT', 'HT29_WT', 'HT29_I', 'HT29_I','HT29_I', 'HT29_V', 'HT29_V', 'HT29_V', 'HT29_AS','HT29_AS', 'HT29_AS') )
    print(myfactors)
    library(NOISeq)
    mydata <- readData( data = mat_all, factors = myfactors)
    if( FALSE )
    {
    cell <- c('PC9_P','PC9_E','PC9_AS')
    for( i in cell)
    {
        print(i)
        myresult <- noiseq(mydata, factor = 'cell', conditions = c('PC9_WT',i), norm = 'n')
        saveRDS(myresult, file = paste('working_data/persister_rna/study_5_',i,'.rds',sep=''))
        write.table( myresult@results[[1]],file = paste('working_data/persister_rna/study_5_',i,'.csv',sep=''),quote = F,row.names = T,col.names = T,sep = ',')
    } 
    }
    cell <- c('HT29_I','HT29_V','HT29_AS')
    for( i in cell)
    {
        print(i)
        myresult <- noiseq(mydata, factor = 'cell', conditions = c('HT29_WT',i), norm = 'n')
        saveRDS(myresult, file = paste('working_data/persister_rna/study_5_',i,'.rds',sep=''))
        write.table( myresult@results[[1]],file = paste('working_data/persister_rna/study_5_',i,'.csv',sep=''),quote = F,row.names = T,col.names = T,sep = ',')
    }
}


#study 8 noiseq
if(FALSE)
{
    mat_all <- readRDS( file = 'working_data/persister_rna/GSE114647.rds')
    print(head(mat_all))
    rownames(mat_all) <- mat_all[,'gene']
    mat_all <- mat_all[,colnames(mat_all) != 'gene']
    print(head(mat_all))
    myfactors = data.frame( cell = c('HCC827_par', 'HCC827_par', 'HCC827_par', 'HCC827_TKI',
       'HCC827_TKI', 'HCC827_TKI', 'HCC4006_par', 'HCC4006_par',
       'HCC4006_par', 'HCC4006_TKI', 'HCC4006_TKI', 'HCC4006_TKI',
       'H1975_par', 'H1975_par', 'H1975_par', 'H1975_TKI',
       'H1975_TKI', 'H1975_TKI', 'PC9_par', 'PC9_par', 'PC9_TKI',
       'PC9_TKI') )
    print(myfactors)
    library(NOISeq)
    mydata <- readData( data = mat_all, factors = myfactors)
    cell <- c('HCC827','HCC4006','H1975','PC9')
    for( i in cell)
    {
        print(i)
        myresult <- noiseq(mydata, factor = 'cell', conditions = c(paste0(i,'_par'),paste0(i,'_TKI')), norm = 'n')
        saveRDS(myresult, file = paste('working_data/persister_rna/study_8_',i,'.rds',sep=''))
        write.table( myresult@results[[1]],file = paste('working_data/persister_rna/study_8_',i,'.csv',sep=''),quote = F,row.names = T,col.names = T,sep = ',')
    }
}

#study 2 edgeR
if(FALSE)
{
	library(Seurat)
    library(edgeR)
    library(magrittr)
    library(plyr)
    library(tidyverse)
    library(Matrix)
	sc <- readRDS(file = './working_data/persister_comb.rds')
    sc <- subset(sc, subset = sample_type != 'cycling')
    set.seed(20000)
	    
    if(TRUE)
    {
        for( i in c('A375_SKIN','PC9_LUNG','HT29_LARGE_INTESTINE','COLO829_SKIN'))
        {
            print(i)
            tmp <- subset(sc, subset = cell_line == i)
            used_genes <- which(Matrix::rowSums(GetAssayData(tmp, 'counts')[, rownames(tmp@meta.data)] > 0) > nrow(tmp@meta.data)*0.01)
            print(length(used_genes))
            dge <- DGEList(GetAssayData(tmp, 'counts')[used_genes,rownames(tmp@meta.data)])
            dge <- calcNormFactors(dge, method = 'TMMwzp')
            df <- tmp@meta.data
            design <- model.matrix(~0 + sample_type, data = df)
            cond2 <- 'non-cycling'
            cond1 <- 'naive'
            treat_group <- as.numeric(grepl(cond2, colnames(design))) %>% magrittr::divide_by(., sum(.))
            control_group <- as.numeric(grepl(cond1, colnames(design))) %>% magrittr::divide_by(., sum(.))
            #get differential treatment effect
            contr_vec <- treat_group-control_group
            dge <- estimateDisp(dge, design = design)
            fit <- glmQLFit(dge, design = design, prior.count = 0.125)
            qlf <- glmQLFTest(fit, contrast = contr_vec)
            tt_res <- topTags(qlf, n = Inf)$table %>% rownames_to_column(var = 'Gene')
            print('success')
            print(head(tt_res))
            saveRDS( tt_res, file = paste0('working_data/persister_rna/study_2_',i,'.rds'))
        }
    }
}

#study 7 edgeR
if(FALSE)
{
    library(edgeR)
    library(magrittr)
    library(plyr)
    library(tidyverse)
    library(Matrix)
    
    mat_all <- readRDS( file = 'working_data/persister_rna/GSE162285.rds')
    print(head(mat_all))
    rownames(mat_all) <- mat_all[,'gene']
    mat_all <- mat_all[,colnames(mat_all) != 'gene']
    #print(head(mat_all))
    
    myfactors = data.frame( cell = c('MDAMB231_WT','MDAMB231_WT','MDAMB231_WT','MDAMB231_D','MDAMB231_D', 'MDAMB231_D') )
    rownames(myfactors) <- colnames(mat_all)
    myfactors$what <- 'what'
    print(head(myfactors))
    
    cell <- c('MDAMB231_D')
    for( i in cell)
    {
        print(i)
        tmp <- mat_all[, rownames(myfactors)[myfactors[,'cell'] %in% c(i,'MDAMB231_WT')] ]
        print(head(tmp))
        #used_genes <- which(Matrix::rowSums(GetAssayData(tmp, 'counts')[, rownames(tmp@meta.data)] > 0) > nrow(tmp@meta.data)*0.05)
        #print(length(used_genes))
        dge <- DGEList(tmp)
        dge <- calcNormFactors(dge, method = 'TMMwzp')
        df <- myfactors[colnames(tmp),]
        design <- model.matrix(~0 + cell, data = df)
        cond2 <- i
        cond1 <- 'MDAMB231_WT'
        treat_group <- as.numeric(grepl(cond2, colnames(design))) %>% magrittr::divide_by(., sum(.))
        control_group <- as.numeric(grepl(cond1, colnames(design))) %>% magrittr::divide_by(., sum(.))
        #get differential treatment effect
        contr_vec <- treat_group-control_group
        dge <- estimateDisp(dge, design = design)
        fit <- glmQLFit(dge, design = design, prior.count = 0.125)
        qlf <- glmQLFTest(fit, contrast = contr_vec)
        tt_res <- topTags(qlf, n = Inf)$table %>% rownames_to_column(var = 'Gene')
        print('success')
        print(head(tt_res))

        saveRDS(tt_res, file = paste('working_data/persister_rna/study_7_edge_',i,'.rds',sep=''))
    } 
}



#study 5 edgeR
if(FALSE)
{
    library(edgeR)
    library(magrittr)
    library(plyr)
    library(tidyverse)
    library(Matrix)
    
    mat_all <- readRDS( file = 'working_data/persister_rna/GSE189625.rds')
    print(head(mat_all))
    rownames(mat_all) <- mat_all[,'gene']
    mat_all <- mat_all[,colnames(mat_all) != 'gene']
    #print(head(mat_all))
    
    myfactors = data.frame( cell = c('PC9_WT', 'PC9_WT', 'PC9_WT', 'PC9_P', 'PC9_P', 'PC9_P','PC9_E', 'PC9_E', 'PC9_E', 'PC9_AS', 'PC9_AS', 'PC9_AS','HT29_WT', 'HT29_WT', 'HT29_WT', 'HT29_I', 'HT29_I','HT29_I', 'HT29_V', 'HT29_V', 'HT29_V', 'HT29_AS','HT29_AS', 'HT29_AS') )
    rownames(myfactors) <- colnames(mat_all)
    myfactors$what <- 'what'
    print(head(myfactors))
    
    cell <- c('PC9_P','PC9_E','PC9_AS')
    for( i in cell)
    {
        print(i)
        tmp <- mat_all[, rownames(myfactors)[myfactors[,'cell'] %in% c(i,'PC9_WT')] ]
        print(head(tmp))
        #used_genes <- which(Matrix::rowSums(GetAssayData(tmp, 'counts')[, rownames(tmp@meta.data)] > 0) > nrow(tmp@meta.data)*0.05)
        #print(length(used_genes))
        dge <- DGEList(tmp)
        dge <- calcNormFactors(dge, method = 'TMMwzp')
        df <- myfactors[colnames(tmp),]
        design <- model.matrix(~0 + cell, data = df)
        cond2 <- i
        cond1 <- 'PC9_WT'
        treat_group <- as.numeric(grepl(cond2, colnames(design))) %>% magrittr::divide_by(., sum(.))
        control_group <- as.numeric(grepl(cond1, colnames(design))) %>% magrittr::divide_by(., sum(.))
        #get differential treatment effect
        contr_vec <- treat_group-control_group
        dge <- estimateDisp(dge, design = design)
        fit <- glmQLFit(dge, design = design, prior.count = 0.125)
        qlf <- glmQLFTest(fit, contrast = contr_vec)
        tt_res <- topTags(qlf, n = Inf)$table %>% rownames_to_column(var = 'Gene')
        print('success')
        print(head(tt_res))

        saveRDS(tt_res, file = paste('working_data/persister_rna/study_5_edge_',i,'.rds',sep=''))
    } 

    cell <- c('HT29_I','HT29_V','HT29_AS')
    for( i in cell)
    {
        print(i)
        tmp <- mat_all[, rownames(myfactors)[myfactors[,'cell'] %in% c(i,'HT29_WT')] ]
        print(head(tmp))
        #used_genes <- which(Matrix::rowSums(GetAssayData(tmp, 'counts')[, rownames(tmp@meta.data)] > 0) > nrow(tmp@meta.data)*0.05)
        #print(length(used_genes))
        dge <- DGEList(tmp)
        dge <- calcNormFactors(dge, method = 'TMMwzp')
        df <- myfactors[colnames(tmp),]
        design <- model.matrix(~0 + cell, data = df)
        cond2 <- i
        cond1 <- 'HT29_WT'
        treat_group <- as.numeric(grepl(cond2, colnames(design))) %>% magrittr::divide_by(., sum(.))
        control_group <- as.numeric(grepl(cond1, colnames(design))) %>% magrittr::divide_by(., sum(.))
        #get differential treatment effect
        contr_vec <- treat_group-control_group
        dge <- estimateDisp(dge, design = design)
        fit <- glmQLFit(dge, design = design, prior.count = 0.125)
        qlf <- glmQLFTest(fit, contrast = contr_vec)
        tt_res <- topTags(qlf, n = Inf)$table %>% rownames_to_column(var = 'Gene')
        print('success')
        print(head(tt_res))

        saveRDS(tt_res, file = paste('working_data/persister_rna/study_5_edge_',i,'.rds',sep=''))
    }
}

#study 8 edgeR
if(FALSE)
{
    library(edgeR)
    library(magrittr)
    library(plyr)
    library(tidyverse)
    library(Matrix)
    
    mat_all <- readRDS( file = 'working_data/persister_rna/GSE189625.rds')
    print(head(mat_all))
    rownames(mat_all) <- mat_all[,'gene']
    mat_all <- mat_all[,colnames(mat_all) != 'gene']
    #print(head(mat_all))
 
    mat_all <- readRDS( file = 'working_data/persister_rna/GSE114647.rds')
    print(head(mat_all))
    rownames(mat_all) <- mat_all[,'gene']
    mat_all <- mat_all[,colnames(mat_all) != 'gene']
    print(head(mat_all))
    myfactors = data.frame( cell = c('HCC827_par', 'HCC827_par', 'HCC827_par', 'HCC827_TKI',
       'HCC827_TKI', 'HCC827_TKI', 'HCC4006_par', 'HCC4006_par',
       'HCC4006_par', 'HCC4006_TKI', 'HCC4006_TKI', 'HCC4006_TKI',
       'H1975_par', 'H1975_par', 'H1975_par', 'H1975_TKI',
       'H1975_TKI', 'H1975_TKI', 'PC9_par', 'PC9_par', 'PC9_TKI',
       'PC9_TKI') )
    print(myfactors)
    rownames(myfactors) <- colnames(mat_all)
    myfactors$what <- 'what'
    
    cell <- c('HCC827','HCC4006','H1975','PC9')
    for( i in cell)
    {
        print(i)
        tmp <- mat_all[, rownames(myfactors)[myfactors[,'cell'] %in% c(paste0(i,'_par'),paste0(i,'_TKI'))] ]
        print(head(tmp))
        dge <- DGEList(tmp)
        dge <- calcNormFactors(dge, method = 'TMMwzp')
        df <- myfactors[colnames(tmp),]
        design <- model.matrix(~0 + cell, data = df)
        cond2 <- paste0(i,'_TKI')
        cond1 <- paste0(i,'_par')
        treat_group <- as.numeric(grepl(cond2, colnames(design))) %>% magrittr::divide_by(., sum(.))
        control_group <- as.numeric(grepl(cond1, colnames(design))) %>% magrittr::divide_by(., sum(.))
        #get differential treatment effect
        contr_vec <- treat_group-control_group
        dge <- estimateDisp(dge, design = design)
        fit <- glmQLFit(dge, design = design, prior.count = 0.125)
        qlf <- glmQLFTest(fit, contrast = contr_vec)
        tt_res <- topTags(qlf, n = Inf)$table %>% rownames_to_column(var = 'Gene')
        print('success')
        print(head(tt_res))

        saveRDS(tt_res, file = paste('working_data/persister_rna/study_8_edge_',i,'.rds',sep=''))
    } 
}

#study 10 edgeR
if(FALSE)
{
    library(Seurat)
    sc_pt <- Read10X_h5('working_data/GSE189638/PC9_PT.h5')
    sc_1d <- Read10X_h5('working_data/GSE189638/PC9_1d.h5')
    sc_3d <- Read10X_h5('working_data/GSE189638/PC9_3d.h5')
    sc_1w <- Read10X_h5('working_data/GSE189638/PC9_1w.h5')
    sc_pt <- CreateSeuratObject(sc_pt, project = 'persister_1')
    sc_1d <- CreateSeuratObject(sc_1d, project = 'persister_2')
    sc_3d <- CreateSeuratObject(sc_3d, project = 'persister_3')
    sc_1w <- CreateSeuratObject(sc_1w, project = 'persister_4')
    sc <- merge( sc_pt, y = c(sc_1d, sc_3d, sc_1w), add.cell.ids = c('PT','1d','3d','1w'), project = 'persister')
    print(sc)
    saveRDS(sc,'working_data/GSE189638/PC9.rds')
}

if(FALSE)
{
    library(Seurat)
    library(edgeR)
    library(magrittr)
    library(plyr)
    library(tidyverse)
    library(Matrix)
    sc <- readRDS('working_data/GSE189638/PC9.rds')
    print(sc)
    sample <- sapply(X = strsplit(colnames(sc), split = "_"), FUN = "[", 1)
    meta.data <- data.frame( sample)
    rownames(meta.data) <- colnames(sc)
    colnames(meta.data) <- 'sample_ID'
    sc@meta.data <- meta.data
    print(head(sc@meta.data))
    print(sc)
    for( i in c('1d','3d','1w'))
    {
        print(i)
        tmp <- subset(sc, subset = (sample_ID == i)| (sample_ID == 'PT'))
        print(tmp)
        used_genes <- which(Matrix::rowSums(GetAssayData(tmp, 'counts')[, rownames(tmp@meta.data)] > 0) > nrow(tmp@meta.data)*0.05)
        print(length(used_genes))
        dge <- DGEList(GetAssayData(tmp, 'counts')[used_genes,rownames(tmp@meta.data)])
        dge <- calcNormFactors(dge, method = 'TMMwzp')
        df <- tmp@meta.data
        design <- model.matrix(~0 + sample_ID, data = df)
        cond2 <- i
        cond1 <- 'PT'
        treat_group <- as.numeric(grepl(cond2, colnames(design))) %>% magrittr::divide_by(., sum(.))
        control_group <- as.numeric(grepl(cond1, colnames(design))) %>% magrittr::divide_by(., sum(.))
        #get differential treatment effect
        contr_vec <- treat_group-control_group
        dge <- estimateDisp(dge, design = design)
        fit <- glmQLFit(dge, design = design, prior.count = 0.125)
        qlf <- glmQLFTest(fit, contrast = contr_vec)
        tt_res <- topTags(qlf, n = Inf)$table %>% rownames_to_column(var = 'Gene')
        print('success')
        print(head(tt_res))
        saveRDS( tt_res, file = paste0('working_data/persister_rna/study_5_',i,'_single_cell.rds'))
    }
}

#study 2 edgeR
if(FALSE)
{
	library(Seurat)
    library(edgeR)
    library(magrittr)
    library(plyr)
    library(tidyverse)
    library(Matrix)
	sc <- readRDS(file = './working_data/persister_comb.rds')
    sc <- subset(sc, subset = sample_type != 'cycling')
    set.seed(20000)
	    
    if(TRUE)
    {
        #for( i in c('A375_SKIN','PC9_LUNG','HT29_LARGE_INTESTINE','COLO829_SKIN'))
        for( i in c('BT474_BREAST'))
        {
            print(i)
            tmp <- subset(sc, subset = cell_line == i)
            used_genes <- which(Matrix::rowSums(GetAssayData(tmp, 'counts')[, rownames(tmp@meta.data)] > 0) > nrow(tmp@meta.data)*0.01)
            print(length(used_genes))
            dge <- DGEList(GetAssayData(tmp, 'counts')[used_genes,rownames(tmp@meta.data)])
            dge <- calcNormFactors(dge, method = 'TMMwzp')
            df <- tmp@meta.data
            design <- model.matrix(~0 + sample_type, data = df)
            cond2 <- 'non-cycling'
            cond1 <- 'naive'
            treat_group <- as.numeric(grepl(cond2, colnames(design))) %>% magrittr::divide_by(., sum(.))
            control_group <- as.numeric(grepl(cond1, colnames(design))) %>% magrittr::divide_by(., sum(.))
            #get differential treatment effect
            contr_vec <- treat_group-control_group
            dge <- estimateDisp(dge, design = design)
            fit <- glmQLFit(dge, design = design, prior.count = 0.125)
            qlf <- glmQLFTest(fit, contrast = contr_vec)
            tt_res <- topTags(qlf, n = Inf)$table %>% rownames_to_column(var = 'Gene')
            print('success')
            print(head(tt_res))
            saveRDS( tt_res, file = paste0('working_data/persister_rna/study_2_',i,'.rds'))
        }
    }
}

if(FALSE)
{
    for( i in c('A375_SKIN','PC9_LUNG','HT29_LARGE_INTESTINE','COLO829_SKIN','BT474_BREAST'))
    {
        print(i)
        res <- readRDS( file = paste0('working_data/persister_rna/study_2_',i,'.rds'))
        print(head(res))
        write.table( res ,file = paste0('working_data/persister_rna/study_2_',i,'.csv'),quote = F,row.names = T,col.names = T,sep = '\t')
    }
}

if(FALSE)
{
    library(ComplexHeatmap)
    library(circlize)
    set.seed(20000)
	cor <- read.csv( paste0('working_data/erastin_rna/cp_cor.csv'), header = TRUE, row.names = 1, sep = ',')
    #print(head(cor))
    #print(head(sample))
    colnames(cor) <- rownames(cor)
    cor <- as.matrix(cor)
    #column_ha = HeatmapAnnotation( study = sample[,'study'], seq = sample[,'seq'])
    col_fun = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
    pdf("figure/erastin_rna_seq_cor_2.pdf", onefile = TRUE, width = 13, height = 10)
    graph <- Heatmap(cor, col = col_fun, column_title = 'cor')
    print(graph)
    dev.off()
}



if(FALSE)
{
    library(ComplexHeatmap)
    library(circlize)
    set.seed(20000)
	cor <- read.csv( paste0('working_data/rna_seq_cor.csv'), header = TRUE, row.names = 1, sep = ',')
	sample <- read.csv( paste0('data/persister_rna/sample'), header = TRUE, row.names = 3, sep = ',')
    #print(head(cor))
    #print(head(sample))
    colnames(cor) <- rownames(cor)
    sample <- sample[rownames(cor),]
    cor <- as.matrix(cor)
    sample[,'study'] <- as.character(sample[,'study'])
    #column_ha = HeatmapAnnotation( study = sample[,'study'], seq = sample[,'seq'])
    column_ha = HeatmapAnnotation( seq = sample[,'seq'])
    row_ha = rowAnnotation(drug = sample[,'drug'])
    col_fun = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
    pdf("figure/rna_seq_cor.pdf", onefile = TRUE, width = 13, height = 10)
    graph <- Heatmap(cor, col = col_fun, column_title = 'cor',top_annotation = column_ha, left_annotation = row_ha)
    print(graph)
    dev.off()
}

if(FALSE)
{
    library(ComplexHeatmap)
    library(circlize)
	cor <- read.csv( paste0('working_data/rna_seq_cor_9.csv'), header = TRUE, row.names = 1, sep = ',')
	sample <- read.csv( paste0('working_data/curate_9.txt'), header = TRUE, sep = '\t')
    #print(head(cor))
    #print(head(sample))
    rownames(sample) <- sample[,'sample_name']
    colnames(cor) <- rownames(cor)
    sample <- sample[rownames(cor),]
    cor <- as.matrix(cor)
    sample[,'study'] <- as.character(sample[,'study'])
    #column_ha = HeatmapAnnotation( study = sample[,'study'], seq = sample[,'seq'])
    column_ha = HeatmapAnnotation( seq = sample[,'seq'])
    row_ha = rowAnnotation(drug = sample[,'drug'])
    col_fun = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
    pdf("figure/rna_seq_cor_9.pdf", onefile = TRUE, width = 13, height = 10)
    graph <- Heatmap(cor, show_row_names = FALSE, show_column_names = FALSE, col = col_fun, column_title = 'cor',top_annotation = column_ha, left_annotation = row_ha)
    print(graph)
    dev.off()
}

#LCD cell lines new GSH & media
if(FALSE)
{
    library(ggplot2)
    library(ggrepel)
    mat <- read.csv( paste('working_data/LCD_ferro_drug.csv'),sep=',');
    print(head(mat))
    #mat[,'LCD'] <- factor(mat[,'LCD'])
    #mat <- mat_all[(mat_all$c >= i),]
    pdf("figure/LCD_ferro_drug.pdf", onefile = TRUE, width = 9, height = 9)

    graph <- ggplot(mat, aes_string(y = 'RSL.3', x = 'erastin')) +
        geom_smooth(method='lm')+
        geom_point(alpha = 0.6, size = 4) +
        #geom_point( aes(fill = culture), colour = 'black', pch = 21, alpha = 0.6, size = 4) +
        #scale_fill_brewer(palette = 'Dark2')+
        #scale_x_continuous(limits = c(4.5,7))+
        #scale_y_continuous(limits = c(0,0.05))+
        #guides( fill=FALSE)+
        geom_text_repel( aes(label = cell), max.overlaps = 50)+
        labs()
    graph <- graph + theme_classic(base_size = 18)
    #graph <- graph + theme_classic(base_size = 18)+
    #        theme(legend.position = c(0.9,0.9))
    print(graph)
    graph <- ggplot(mat, aes_string(y = 'RSL.3', x = 'ML210')) +
        geom_smooth(method='lm')+
        geom_point(alpha = 0.6, size = 4) +
        #geom_point( aes(fill = culture), colour = 'black', pch = 21, alpha = 0.6, size = 4) +
        #scale_fill_brewer(palette = 'Dark2')+
        #scale_x_continuous(limits = c(4.5,7))+
        #scale_y_continuous(limits = c(0,0.05))+
        #guides( fill=FALSE)+
        geom_text_repel( aes(label = cell), max.overlaps = 50)+
        labs()
    graph <- graph + theme_classic(base_size = 18)
    #graph <- graph + theme_classic(base_size = 18)+
    #        theme(legend.position = c(0.9,0.9))
    print(graph)
    graph <- ggplot(mat, aes_string(y = 'RSL.3', x = 'proliferation')) +
        geom_smooth(method='lm')+
        geom_point(alpha = 0.6, size = 4) +
        #geom_point( aes(fill = culture), colour = 'black', pch = 21, alpha = 0.6, size = 4) +
        #scale_fill_brewer(palette = 'Dark2')+
        #scale_x_continuous(limits = c(4.5,7))+
        #scale_y_continuous(limits = c(0,0.05))+
        #guides( fill=FALSE)+
        geom_text_repel( aes(label = cell), max.overlaps = 50)+
        labs()
    graph <- graph + theme_classic(base_size = 18)
    #graph <- graph + theme_classic(base_size = 18)+
    #        theme(legend.position = c(0.9,0.9))
    print(graph)
    graph <- ggplot(mat, aes_string(y = 'erastin', x = 'proliferation')) +
        geom_smooth(method='lm')+
        geom_point(alpha = 0.6, size = 4) +
        #geom_point( aes(fill = culture), colour = 'black', pch = 21, alpha = 0.6, size = 4) +
        #scale_fill_brewer(palette = 'Dark2')+
        #scale_x_continuous(limits = c(4.5,7))+
        #scale_y_continuous(limits = c(0,0.05))+
        #guides( fill=FALSE)+
        geom_text_repel( aes(label = cell), max.overlaps = 50)+
        labs()
    graph <- graph + theme_classic(base_size = 18)
    #graph <- graph + theme_classic(base_size = 18)+
    #        theme(legend.position = c(0.9,0.9))
    print(graph)
    col <- unique(mat[,'culture'])
    for( i in col )
    {
        if( i != '')
        {
        mat_2 <- mat
        mat_2[mat[,'culture'] != i,'culture'] <- ''
        mat_2[mat[,'culture'] != i,'cell'] <- ''
        graph <- ggplot(mat_2, aes(y = RSL.3, x = erastin)) +
            geom_point( aes(fill = culture), colour = 'black', pch = 21, alpha = 0.6, size = 3) +
            scale_fill_manual(values = c("black","red"))+
            geom_text_repel( aes(label = cell))+
            labs()
        graph <- graph + theme_classic(base_size = 18)+
            theme(legend.position = c(0.2,0.9))
        print(graph)
        }
    }
    dev.off()
}


#LCD cell lines new GSH & media
if(FALSE)
{
    library(ggplot2)
    library(ggrepel)
    mat <- read.csv( paste('table/growth_metabolites_LCD.csv'),sep=',',row.names = 1);
    print(head(mat))
    pdf("figure/LCD_GSH_NADP.pdf", onefile = TRUE, width = 9, height = 9)
    graph <- ggplot(mat, aes(y = NADP, x = GSH)) +
        #geom_smooth(method='lm')+
        #geom_point(alpha = 0.6, size = 4) +
        geom_point( aes(color = redox_sel), alpha = 0.6, size = 4) +
        #scale_color_manual( values = c('False' = 'azure2', 'True' = 'cornflowerblue') )+
        scale_color_manual( values = c('other' = 'azure2', 'Low' = 'cornflowerblue', 'Medium' = 'coral', 'High' = 'red') )+
        #scale_fill_brewer(palette = 'Dark2')+
        #scale_x_continuous(limits = c(4.5,7))+
        #scale_y_continuous(limits = c(0,0.05))+
        #guides( color=FALSE)+
        geom_text_repel( aes(label = redox), max.overlaps = 50)+
        labs()
    graph <- graph + theme_classic(base_size = 18)
            #theme(legend.position = c(0.9,0.9))
    print(graph)
    graph <- ggplot(mat, aes(y = NADP, x = GSH)) +
        #geom_smooth(method='lm')+
        #geom_point(alpha = 0.6, size = 4) +
        #geom_point( aes(fill = sel), colour = 'black', pch = 21, alpha = 0.6, size = 4) +
        geom_point( aes(color = redox_sel), alpha = 0.6, size = 4) +
        #scale_color_manual( values = c('False' = 'azure2', 'True' = 'cornflowerblue') )+
        scale_color_manual( values = c('other' = 'azure2', 'Low' = 'cornflowerblue', 'Medium' = 'coral', 'High' = 'red') )+
        #scale_fill_brewer(palette = 'Dark2')+
        #scale_x_continuous(limits = c(4.5,7))+
        #scale_y_continuous(limits = c(0,0.05))+
        #guides( color=FALSE)+
        geom_label_repel( aes(label = redox), max.overlaps = 50)+
        labs()
    graph <- graph + theme_classic(base_size = 18)
            #theme(legend.position = c(0.9,0.9))
    print(graph)
    graph <- ggplot(mat, aes(y = NADP, x = GSH)) +
        #geom_smooth(method='lm')+
        #geom_point(alpha = 0.6, size = 4) +
        geom_point( aes(color = sel), alpha = 0.6, size = 4) +
        #scale_color_manual( values = c('False' = 'azure2', 'Low' = 'cornflowerblue', 'Medium' = 'coral', 'High' = 'red') )+
        scale_color_manual( values = c('False' = 'azure2', 'True' = 'cornflowerblue') )+
        #scale_fill_brewer(palette = 'Dark2')+
        #scale_x_continuous(limits = c(4.5,7))+
        #scale_y_continuous(limits = c(0,0.05))+
        #guides( fill=FALSE)+
        geom_text_repel( aes(label = LCD), max.overlaps = 50)+
        labs()
    graph <- graph + theme_classic(base_size = 18)+
            theme(legend.position = c(0.9,0.9))
    print(graph)
    graph <- ggplot(mat, aes(y = NADP, x = GSH)) +
        #geom_smooth(method='lm')+
        #geom_point(alpha = 0.6, size = 4) +
        #geom_point( aes(fill = sel), colour = 'black', pch = 21, alpha = 0.6, size = 4) +
        geom_point( aes(color = sel), alpha = 0.6, size = 4) +
        scale_color_manual( values = c('False' = 'azure2', 'True' = 'cornflowerblue') )+
        #scale_fill_brewer(palette = 'Dark2')+
        #scale_x_continuous(limits = c(4.5,7))+
        #scale_y_continuous(limits = c(0,0.05))+
        #guides( fill=FALSE)+
        geom_label_repel( aes(label = LCD), max.overlaps = 50)+
        labs()
    graph <- graph + theme_classic(base_size = 18)+
            theme(legend.position = c(0.9,0.9))
    print(graph)
    dev.off()
}

if(FALSE)
{
    library(ComplexHeatmap)
    library(circlize)
    set.seed(20000)
	cor <- read.csv( paste0('space_data/keap1_xpr_cor_wide.csv'), header = TRUE, row.names = 1, sep = ',')
    #print(head(cor))
    #print(head(sample))
    colnames(cor) <- rownames(cor)
    cor <- as.matrix(cor)
    #column_ha = HeatmapAnnotation( study = sample[,'study'], seq = sample[,'seq'])
    col_fun = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
    pdf("figure/keap1_xpr_cor_wide.pdf", onefile = TRUE, width = 20, height = 20)
    graph <- Heatmap(cor, col = col_fun, column_title = 'cor')
    print(graph)
    dev.off()
}

if(FALSE)
{
    library(ComplexHeatmap)
    library(circlize)
    set.seed(20000)
    pdf("figure/keap1_xpr_cor_all_individual.pdf", onefile = TRUE, width = 15, height = 15)
    dep <- read.csv( paste0('space_data/L1000_gene_cor/keap1_dep'), header = TRUE, row.names = 1, sep = ',')
    for( i in colnames(dep))
    {
        cor <- read.csv( paste0('space_data/L1000_gene_cor/',i), header = TRUE, row.names = 1, sep = ',')
        tmpdep <- dep[rownames(cor),]
        colnames(cor) <- rownames(cor)
        cor <- as.matrix(cor)
        drug_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
        column_ha = HeatmapAnnotation( Dep = tmpdep[,i], col = list(Dep = drug_fun))
        col_fun = colorRamp2(c(-0.7, 0, 0.7), c("blue", "white", "red"))
        graph <- Heatmap(cor, col = col_fun, top_annotation = column_ha, column_title = i)
        print(graph)
    }

    dep <- read.csv( paste0('space_data/L1000_gene_cor/keap1_wide_dep'), header = TRUE, row.names = 1, sep = ',')
    for( i in colnames(dep))
    {
        cor <- read.csv( paste0('space_data/L1000_gene_cor/',i), header = TRUE, row.names = 1, sep = ',')
        tmpdep <- dep[rownames(cor),]
        colnames(cor) <- rownames(cor)
        cor <- as.matrix(cor)
        drug_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
        column_ha = HeatmapAnnotation( Dep = tmpdep[,i], col = list(Dep = drug_fun))
        col_fun = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
        graph <- Heatmap(cor, col = col_fun, top_annotation = column_ha, column_title = i)
        print(graph)
    }
    dev.off()
}

