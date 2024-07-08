#all breast cancer
if(FALSE)
{
    library(ggplot2)
    library(ggrepel)
    library(cowplot)
    pdf("figure/Fig_1_lower.pdf", onefile = TRUE, width = 20, height = 10)
    #all breast gene priority demonstration
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

#subtype
#Breast cancer subtype percent-median gene effect plot
if(TRUE)
{
    library(ggplot2)
    library(ggrepel)
    library(cowplot)
    print('Fig_2')
    pdf("figure/Fig_2.pdf", onefile = TRUE, width = 17, height = 17)
    mat <- read.csv( paste('R_data/all_gene.csv'),sep=',',row.names = 1);
    print(head(mat))
    mat$type <- factor(mat$assignment)
    mat$molecular <- factor(mat$molecular)
    graph <- ggplot(mat, aes(y = -median_eff, x = perc)) +
	facet_wrap(~molecular, nrow = 1)+
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
    p1 <- graph + theme_classic(base_size = 18)+theme(legend.position = c(0.1,0.9))
    #graph <- graph + theme_classic(base_size = 18)+theme(legend.position = c(0.1,0.9))

    mat <- read.csv( paste('R_data/breast_priority_open_target.csv'),sep=',',row.names = 1);
    #mat <- read.csv( paste('R_data/all_gene_PIK3CA_2.csv'),sep=',',row.names = 1);
    print(head(mat))
    mat$molecular <- factor(mat$molecular)
    mat$sel <- factor(mat$sel)
    graph <- ggplot(mat, aes(y = -median_eff, x = perc)) +
	facet_wrap(~molecular, nrow = 1)+
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
        labs(x = 'fraction of cell line', y = 'median genetic effect score')
    p2 <- graph + theme_classic(base_size = 18)+theme(legend.position = c(0.1,0.9))
    graph <- plot_grid(p1,p2, labels = c('B','C'), label_size = 30,ncol = 1)
    print(graph)
    dev.off()
}

#enrichment and heatmap
if(FALSE)
{
	library(clusterProfiler)
	library(org.Hs.eg.db)
	library(enrichplot)
	library(ggplot2)
	mat <- read.csv( paste('R_data/all_gene_PIK3CA.csv'),sep=',',row.names = 1);
    	#pdf("figure/Fig_3.pdf", onefile = TRUE, width = 25, height = 20)
    	pdf("figure/Fig_3.pdf", onefile = TRUE, width = 17, height = 25)
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
    col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
    graph_h <- Heatmap(cor, show_row_names = TRUE, show_column_names = TRUE, col =col_fun,heatmap_legend_param = list(title = 'dep', at = c(-2,0,2), labels = c(2,0,-2)), name = 'priority_heatmap', column_title = 'HIG_subtype',top_annotation = column_ha, left_annotation = row_ha)
    gb <- grid.grabExpr(draw(graph_h))
    #gb <- cowplot::plot_grid(gb,NULL, labels = c('D',''), label_size = 30, ncol = 2, rel_widths = c(2,1))
    #gb <- cowplot::plot_grid(gb,NULL, labels = c('D',''), label_size = 30, ncol = 2,)
   graph <- cowplot::plot_grid(graph,graph_2,gb, labels = c('','','D'), label_size = 30, ncol = 1, rel_heights = c(0.8,0.8,1))
	print(graph)
	dev.off()
	#saveRDS(kk,file = 'working_data/BT474_6hr_oraKEGG_down_more.rds') 
}

#mutation association
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
