library(ComplexHeatmap)
library(circlize)
library(ggfortify)
library(scales)
library(ggrepel)
library(dplyr)
library(lattice)
library(DESeq2)
library(data.table)

# Step 1: prepare_experiment, returns basic experiment object
# Step 2: add_condition_to_experiment, called on the returned experiment object
# Step 3: prepare_irmatrix, called on the complete experiment object
# Step 4: call glm with parameters being the experiment object, and irfinder_tables from prepare_irmatrix
# Step 5: with experiment object, glm object, irmatrix where required, call desired plotting functions

# Start here: create an experiment definition (just sample names)
prepare_experiment <- function(sample_names, sample_name_column) {
	experiment = data.frame(factor(sample_names))
	names(experiment) = c(sample_name_column)
	return(experiment)
}

# Add conditions to an experiment object
add_condition_to_experiment <- function(experiment, sample_conditions, condition_name) {
	experiment = cbind(experiment, factor(sample_conditions))
	names(experiment)[[length(names(experiment))]] = condition_name
	return(experiment)
}

# Prepare IR data associated with experiment object by reading in IRFinder output files
# Returns (irmatrix, irfinder_tables) (irmatrix is a matrix of IR ratios from each sample only, irfinder_tables is a list of full IRFinder outputs for each sample)
prepare_irmatrix <- function(experiment, name_by, irfinder_file_names) {
    irfinder_tables_list = read_irfinder_multiple(irfinder_file_names)
    prepared = matrix(irfinder_tables_list[[1]]$irratio)
    rownames(prepared) = get_unique_transcripts(irfinder_tables_list[[1]])
    for (i in 2:length(irfinder_tables_list)) {
        prepared = cbind(prepared, irfinder_tables_list[[i]]$irratio)
    }
    colnames(prepared) = as.character(experiment[,name_by])
    return(list(prepared, irfinder_tables_list))
}

# Reads IRFinder output from a specified file into a frame
read_irfinder <- function(filename) {
    print(paste0("reading ", filename))
    frame <- read.delim(file=filename, header=FALSE, stringsAsFactors=FALSE)
    names(frame) = c(    
        "chromosome",
        "start",
        "end",
        "gene",
        "empty1",
        "strand",
        "exclbases",
        "coverage",
        "introndepth",
        "intron25",
        "intron50",
        "intron75",
        "readsleft",
        "readsright",
        "depthfirst50",
        "depthlast50",
        "splicesleft",
        "splicesright",
        "splicesexact",
        "irratio",
        "warning"
    )
    return(frame)
}

# Reads IRFinder output from specified files into a frame
read_irfinder_multiple <- function(irfinder_names) {
    return(lapply(irfinder_names, read_irfinder))
}

# Performs DESeq2 based GLM differential IR analysis. Parameters are an experiment definition, how the samples are named, the variable to contrast on, and a list of IRFinder output tables (generated from read_irfinder_multiple)
# Modified version of how to perform GLM with Replicates code, available here: https://github.com/williamritchie/IRFinder/wiki/Generalized-Linear-Model-with-Replicates
glm <- function(raw_experiment, name_by, contrast_by, contrast, irfinder_tables_list) {
    experiment = data.frame(raw_experiment[,name_by], raw_experiment[,contrast_by])
    names(experiment) = c(name_by, contrast_by)
    metaList = DESeqDataSetFromIRFinderTable(irfinder_tables_list, experiment)
    dds = metaList$DESeq2Object
    design(dds) = as.formula(paste0("~", contrast_by, " + ", contrast_by, ":IRFinder"))

    dds = DESeq(dds, parallel = FALSE)

    return(results(dds, contrast=list(paste0(contrast_by, contrast[[1]], ".IRFinderIR"), paste0(contrast_by, contrast[[2]], ".IRFinderIR"))))
}

# Helper function to obtain more readable IR names
get_unique_transcripts <- function(frame) {
    return(unname(apply(as.matrix(frame),1,FUN=function(x){return(paste0(x[4],"/",x[1],":",x[2],"-",x[3],":",x[6]))})))
}

# Internal function assisting in loading data into DESeq2
# Modified version of part of IRFinder source code, available here: https://github.com/williamritchie/IRFinder/blob/master/bin/DESeq2Constructor.R
DESeqDataSetFromIRFinderTable <- function(irfinder_tables, designMatrix, designFormula=~1) {
    res = c()
    libsz = c()
    irnames = get_unique_transcripts(irfinder_tables[[1]])
    n = 1
    for (irtab in irfinder_tables){
        tmp1 = as.numeric(as.vector(irtab[,"introndepth"]))
        tmp2 = as.numeric(as.vector(irtab[,"splicesexact"]))
        res = cbind(res,tmp1)
        libsz = cbind(libsz,tmp2)
        n = n+1
    }
    res.rd=round(res)
    libsz.rd=round(libsz)
    colnames(res.rd)=paste("intronDepth",as.vector(designMatrix[,1]),sep=".")
    rownames(res.rd)=irnames
    colnames(libsz.rd)=paste("totalSplice",as.vector(designMatrix[,1]),sep=".")
    rownames(libsz.rd)=irnames
    
    ir=c(rep("IR",dim(designMatrix)[1]),rep("Splice",dim(designMatrix)[1]))
    group=rbind(designMatrix,designMatrix)
    group$IRFinder=ir
    group$IRFinder=factor(group$IRFinder,levels=c("Splice","IR"))

    counts.IRFinder=cbind(res.rd,libsz.rd)

    dd = DESeqDataSetFromMatrix(countData = counts.IRFinder, colData = group, design = designFormula)
    sizeFactors(dd)=rep(1,dim(group)[1])
    rownames(dd)=irnames
    sp=libsz
    final=list(dd,res,sp)
    names(final)=c("DESeq2Object","IntronDepth","SpliceDepth")
    return(final)
}

# Filters differential IR results for significance
glm_filter_results <- function(original_results) {
    results = original_results[complete.cases(original_results),]
    results = results[results$padj < 0.05,]
    results = results[order(results$padj),]
    return(results)
}

# Generate a HTML report based on differential IR results
generate_report <- function(glm_results, n_dir_results=10, irmatrix, experiment, name_by, contrast_by, contrast, show_results=TRUE, show_volcano=TRUE, show_pca=TRUE, show_exp_table=TRUE) {
	tmp = head(glm_filter_results(glm_results), n=n_dir_results)
	results_table = cbind(tmp$log2FoldChange, -log10(tmp$padj))
	colnames(results_table) = c("log2 Fold Change", "-log10 adjusted p-value")
	rownames(results_table) = rownames(tmp)

    # Configure this with the path to the Rmd template for rendering
	rmarkdown::render("/home/williamw/iranalysis/R/render.Rmd")
}

# Helper function for formatting intron names
shorten_ir_name_irfinder <- function(long_name) {
	s = tstrsplit(long_name, split="/")
	t = tstrsplit(long_name, split=":")
	return(paste0(s[[1]]," ",gsub(" ","",t[[2]],fixed=TRUE)))
}

# Helper function for formatting intron names
shorten_ir_name <- function(long_name) {
	s = tstrsplit(long_name, split="/")
	return(paste(
		s[[1]],
		ifelse(s[[4]]=="IR","Int","Exit"),
		ifelse(s[[6]]=="+",s[[2]],s[[3]]),
		sep="-"))
}

# Obtain the gene name from an intron name
gene_ir_name <- function(long_name) {
	s = tstrsplit(long_name, split="/")
	return(s[[5]])
}

# Generate IR scatter plot. Parameters are glm object, irmatrix, experiment definition. Optionally supports hexbins
ir_scatter_plot <- function(dir_results, irmatrix, experiment, name_by, contrast_by, contrast, hex=FALSE, hex_bins=50, scale_log=FALSE) {
	if (length(contrast) != 2) {
		print("contrast must be == 2")
		return()
	}
	exp_entries1 = experiment[experiment[,contrast_by] == contrast[1],]
	exp_entries2 = experiment[experiment[,contrast_by] == contrast[2],]
	samples1 = as.character(exp_entries1[,name_by])
	samples2 = as.character(exp_entries2[,name_by])
	irmeans = cbind(rowMeans(irmatrix[,samples1]), rowMeans(irmatrix[,samples2]))
	colnames(irmeans) = contrast
	scatter = merge(x = dir_results, y = irmeans, by = "row.names")
	rownames(scatter) = scatter$Row.names
	scatter = scatter[complete.cases(scatter),]
	scatter = scatter[scatter[[contrast[1]]] > 0,]
	scatter = scatter[scatter[[contrast[2]]] > 0,]
	if (hex) {
		if (scale_log) {
			g = ggplot(scatter, aes_string(contrast[1], contrast[2], colour="padj")) + geom_hex(aes(fill=cut(scatter$padj, c(-Inf, 0.001, 0.01, 0.05, Inf))), bins=hex_bins) + scale_x_log10(limits=c(0.01,1), breaks = c(0.01, 0.1, 0.4, 0.7, 1)) + scale_y_log10(limits=c(0.01, 1), breaks=c(0.01, 0.1, 0.4, 0.7, 1)) + scale_fill_manual(name="Adjusted p-value", values=c("(-Inf,0.001]" = "red", "(0.001,0.01]" = "orange", "(0.01,0.05]" = "yellow", "(0.05, Inf]" = "grey"), labels=c(expression(""<=0.001), expression(""<=0.01), expression(""<=0.05), expression("">0.05))) + xlab(paste0(contrast_by, " ", contrast[[2]], " IR-ratio")) + ylab(paste0(contrast_by, " ", contrast[[1]], " IR-ratio"))
		} else {
			g = ggplot(scatter, aes_string(contrast[1], contrast[2], colour="padj")) + geom_hex(aes(fill=cut(scatter$padj, c(-Inf, 0.001, 0.01, 0.05, Inf))), bins=hex_bins) + scale_x_continuous(limits=c(0,1), breaks=c(0.1, 0.4, 0.7, 1)) + scale_y_continuous(limits=c(0,1), breaks=c(0.1, 0.4, 0.7, 1)) + scale_fill_manual(name="Adjusted p-value", values=c("(-Inf,0.001]" = "red", "(0.001,0.01]" = "orange", "(0.01,0.05]" = "yellow", "(0.05, Inf]" = "grey"), labels=c(expression(""<=0.001), expression(""<=0.01), expression(""<=0.05), expression("">0.05))) + xlab(paste0(contrast_by, " ", contrast[[2]], " IR-ratio")) + ylab(paste0(contrast_by, " ", contrast[[1]], " IR-ratio"))
		}
	} else {
		if (scale_log) {
			g = ggplot(scatter, aes_string(contrast[1], contrast[2])) + geom_point(aes(colour=cut(scatter$padj, c(-Inf,0.001,0.01,0.05,Inf)))) + scale_x_log10(limits=c(0.01,1), breaks = c(0.01, 0.1, 0.4, 0.7, 1)) + scale_y_log10(limits=c(0.01, 1), breaks=c(0.01, 0.1, 0.4, 0.7, 1)) + scale_colour_manual(name="Adjusted p-value", values=c("(-Inf,0.001]" = "red", "(0.001,0.01]" = "orange", "(0.01,0.05]" = "yellow", "(0.05, Inf]" = "grey"), labels=c(expression(""<=0.001), expression(""<=0.01), expression(""<=0.05), expression("">0.05))) + xlab(paste0(contrast_by, " ", contrast[[2]], " IR-ratio")) + ylab(paste0(contrast_by, " ", contrast[[1]], " IR-ratio"))
		} else {
			g = ggplot(scatter, aes_string(contrast[1], contrast[2])) + geom_point(aes(colour=cut(scatter$padj, c(-Inf,0.001,0.01,0.05,Inf)))) + scale_x_continuous(limits=c(0,1), breaks=c(0.1, 0.4, 0.7, 1)) + scale_y_continuous(limits=c(0,1), breaks=c(0.1,0.4,0.7,1)) + scale_colour_manual(name="Adjusted p-value", values=c("(-Inf,0.001]" = "red", "(0.001,0.01]" = "orange", "(0.01,0.05]" = "yellow", "(0.05, Inf]" = "grey"), labels=c(expression(""<=0.001), expression(""<=0.01), expression(""<=0.05), expression("">0.05))) + xlab(paste0(contrast_by, " ", contrast[[2]], " IR-ratio")) + ylab(paste0(contrast_by, " ", contrast[[1]], " IR-ratio"))
		}
	}
	g
}

# Generate gene expression vs IR plot. Parameters are glm object, differential gene expression object, experiment definition, and graphical parameters for x and y axis limits
ir_gene_vs_irratio_processed <- function(dir_results, dir_gene, irmatrix, experiment, name_by, contrast_by, contrast, xlimit = 1, ylimit = 10) {
	exp_entries1 = experiment[experiment[,contrast_by] == contrast[1],]
	exp_entries2 = experiment[experiment[,contrast_by] == contrast[2],]
	samples1 = as.character(exp_entries1[,name_by])
	samples2 = as.character(exp_entries2[,name_by])
	cond_mean1 = rowMeans(irmatrix[,samples1])
	cond_mean2 = rowMeans(irmatrix[,samples2])
	cond_mean_diff = data.frame(matrix(cond_mean1 - cond_mean2))
	rownames(cond_mean_diff) = rownames(irmatrix)
	colnames(cond_mean_diff) = c("delta_ir")
	dir_results = data.frame(dir_results)
	with_ir_padj = merge(x = dir_results, y = cond_mean_diff, by = "row.names")
	with_ir_padj = with_ir_padj %>% select(Row.names, padj, delta_ir)
	colnames(with_ir_padj) = c("Row.names", "irpadj", "delta_ir")
	with_ir_padj$gene_id = lapply(with_ir_padj$Row.names, gene_ir_name)
	processed = merge(x = dir_gene, y = with_ir_padj, by = "gene_id")
	processed = processed[complete.cases(processed),]
	ggplot(processed, aes(x=delta_ir, y=log2FoldChange)) + geom_point(aes(colour=cut(processed$irpadj, c(-Inf,0.001,0.01,0.05,Inf)))) + scale_x_continuous(limits=c(-xlimit, xlimit)) + scale_y_continuous(limits=c(-ylimit, ylimit)) + scale_color_manual(name="Differential IR adjusted p-value", values = c("(-Inf,0.001]" = "red", "(0.001,0.01]" = "orange", "(0.01,0.05]" = "yellow", "(0.05, Inf]" = "grey"), labels=c(expression(""<=0.001), expression(""<=0.01), expression(""<=0.05), expression("">0.05))) + xlab("Delta IR-Ratio") + ylab("log2 Fold Change Gene Expression")
}

# Generate IR volcano plot. Parameters are glm object, irmatrix (if using delta IR), experiment definition, graphical parameters, and whether to use delta IR metric (if not, fold change from glm object is used)
ir_volcano <- function(dir_results, irmatrix="", experiment, name_by, contrast_by, contrast, xlimit = 10, ylimit = 20, x_sig_limit = log2(3), y_sig_limit = log10(20), x_label_limit = log2(10), y_label_limit = 10, use_delta_ir = FALSE) {
	if (use_delta_ir) {
		exp_entries1 = experiment[experiment[,contrast_by] == contrast[1],]
		exp_entries2 = experiment[experiment[,contrast_by] == contrast[2],]
		samples1 = as.character(exp_entries1[,name_by])
		samples2 = as.character(exp_entries2[,name_by])
		cond_mean1 = rowMeans(irmatrix[,samples1])
		cond_mean2 = rowMeans(irmatrix[,samples2])
		cond_mean_diff = matrix(cond_mean1 - cond_mean2) # Check order
		rownames(cond_mean_diff) = rownames(irmatrix)
		colnames(cond_mean_diff) = c("delta_ir")
		processed = merge(x = dir_results, y = cond_mean_diff, by = "row.names")
		final = cbind(processed$delta_ir, -log10(processed$padj))
		rownames(final) = processed$Row.names
	} else {
		final = cbind(dir_results$log2FoldChange, -log10(dir_results$padj))
		rownames(final) = rownames(dir_results)
	}
	xaxis_label = ifelse(use_delta_ir, "Delta IR-Ratio", "log2 Fold Change")
	colnames(final) = c("foldchange", "p")
	final = final[complete.cases(final),]
	rownames(final) = lapply(rownames(final), shorten_ir_name_irfinder)
	final = data.frame(final)
	final$out = ifelse(final$p > ylimit, "yout", ifelse(abs(final$foldchange) > xlimit, "xout", "in"))
	final$sig = ifelse(final$p > y_sig_limit & abs(final$foldchange) > x_sig_limit, "sig", "nosig")
	final$labels = ifelse(final$p > y_label_limit & abs(final$foldchange) > x_label_limit, rownames(final), "")
	g = ggplot(final, aes(foldchange, p, label=labels, colour = sig, fill = sig)) + geom_point(aes(shape=final$out)) + scale_shape_manual(values=c("yout" = 24, "xout" = 23, "in" = 21)) + scale_x_continuous(limits=c(-xlimit,xlimit), oob=squish) + scale_y_continuous(limits=c(0,ylimit), oob=squish) + scale_color_manual(values = c("nosig" = "grey", "sig" = "red")) + scale_fill_manual(values = c("nosig" = "grey", "sig" = "red")) + geom_hline(yintercept = y_sig_limit, colour="black") + geom_vline(xintercept = x_sig_limit, colour="black") + geom_vline(xintercept = -x_sig_limit, colour="black") + xlab(xaxis_label) + ylab("-log10 adjusted p value") + theme(legend.position = "none") + geom_text_repel(max.iter=100) + geom_vline(xintercept = x_label_limit, colour="grey") + geom_vline(xintercept = -x_label_limit, colour="grey") + geom_hline(yintercept = y_label_limit, colour="grey")
	g
}

# Generate pca plot. Parameters are irmatrix, experiment definition, and graphical parameters
ir_pca <- function(irmatrix, experiment, contrast_by, scale=TRUE, label=TRUE) {
	mat = irmatrix[complete.cases(irmatrix),]
	if (scale) {
		mat = mat[apply(mat, 1, var) != 0,]
	}
	mat = t(mat)
	pca = prcomp(mat, center = TRUE, scale. = scale)
	if (label) {
		autoplot(pca, data = experiment, colour = contrast_by, label = TRUE, label.repel = TRUE, label.show.legend = FALSE)
	} else {
		autoplot(pca, data = experiment, colour = contrast_by, label = FALSE)
	}
}

# Generate IR heatmap. Parameters are glm object, irmatrix, experiment definition, graphical parameters, and whether to use delta IR metric and clustering columns
ir_heatmap_plot <- function(dir_results, irmatrix, experiment, name_by, contrast_by, contrast, trim=50, use_delta_ir=FALSE, colors = c("blue","white","red"), fontsize = 7, cluster_columns = FALSE) {
	if (length(contrast) <= 1) {
		print("contrast must be > 2")
		return()
	} else if (use_delta_ir) {
		if (length(contrast) == 2) {
			exp_entries1 = experiment[experiment[,contrast_by] == contrast[1],]
			exp_entries2 = experiment[experiment[,contrast_by] == contrast[2],]
			samples1 = as.character(exp_entries1[,name_by])
			samples2 = as.character(exp_entries2[,name_by])
			samples = c(samples1, samples2)
			conditions = c(as.character(exp_entries1[,contrast_by]), as.character(exp_entries2[,contrast_by]))
			cond_mean1 = rowMeans(irmatrix[,samples1])
			cond_mean2 = rowMeans(irmatrix[,samples2])
			cond_mean_diff = abs(cond_mean1 - cond_mean2)
			vars = cond_mean_diff
		} else {
			print("if use_delta_ir, contrast must be == 2")
			return()
		}
	} else {
		exp_entries = experiment[experiment[,contrast_by] %in% contrast,]
		samples = as.character(exp_entries[,name_by])
		conditions = as.character(exp_entries[,contrast_by])
		vars = apply(irmatrix[,samples], 1, var)
	}
	irmatrix = cbind(irmatrix, vars)
	mats = merge(x = dir_results, y = irmatrix, by = "row.names")
	mats = mats %>% select(Row.names, padj, one_of(samples), vars)
	rownames(mats) = mats$Row.names
	rownames(mats) = lapply(rownames(mats), shorten_ir_name_irfinder)
	mats = mats[complete.cases(mats),]
	mats = mats[order(mats$vars, decreasing=TRUE),]
	mats = mats[mats$padj < 0.05,]
	mats = head(mats, trim)
	mats = mats %>% select(-padj, -vars, -Row.names)
	ha = HeatmapAnnotation(data.frame(Condition = conditions))
	Heatmap(mats, col = colorRamp2(c(0,0.5,1), colors), row_names_gp=gpar(fontsize=fontsize), column_names_gp=gpar(fontsize=fontsize), heatmap_legend_param=list(title="IR-Ratio"), cluster_columns=cluster_columns, top_annotation = ha)
}

# Helper function for creating IR boxplots for Gviz visualisation
ir_anno <- function(start, end, strand, chromosome, identifier, index, GdObject, GdObject.original, data, experiment, name_by) {
	data = data[start(data) == start]
	data = data[end(data) == end]
	data = data[strand(data) == strand]
	if (length(data) == 1) {
		irs = as.vector(as.data.frame(data)[levels(experiment[,name_by])][1,])
		rownames(irs) = c("ir")
		irs = t(irs)
		irs = data.frame(irs)
		irs$x = 1
		print(bwplot(ir~x, data=irs, horiz=FALSE, xlab = NULL, ylab = NULL, scales=list(x=list(draw=FALSE)), ylim=c(0,1)), newpage = FALSE, prefix = "plot")
	}
}

# Helper function for Gviz visualisation
bad_chromosome <- function(bad_name) {
	return(gsub("chr", "", bad_name))
}

# Proof of concept Gviz visualisation function, based on and consuming output of internal lab code
dev_RedrawIR <- function(IntronModel = NULL, DataModel = NULL, genome = "hg38", Name.IntronModel = "Introns", Name.DataTrack = "Intron Data", DataTrack.Type = "boxplot", DataTrack.Anno = NULL, Zmin = 0, Zmax = 0, EvenSpacedDataTrack = FALSE, SpacedDataTrackRatio = 0.6, experiment = NULL, name_by = NULL, bamfile = NULL, force_bad_chr_names = FALSE) {
	if (force_bad_chr_names) {
		if (!is.null(IntronModel)) {
			IntronModel$chromosome = bad_chromosome(IntronModel$chromosome)
		}
		if (!is.null(DataModel)) {
			seqlevels(DataModel) = bad_chromosome(seqlevels(DataModel))
			levels(seqnames(DataModel)) = bad_chromosome(levels(seqnames(DataModel)))
		}
	}
	
	tracks <- list()
	trackcount = 1
	JnChr = 0
	if(!is.null(IntronModel)) {
		Final.GR <- with(IntronModel, GRanges(chromosome,IRanges(intron_start,intron_end), strand=strand))
		if(JnChr == 0) {
			JnChr <- IntronModel$chromosome[!duplicated(IntronModel$chromosome)]
		}
		DataModel.selected <- DataModel[seqnames(DataModel) == JnChr]
		zZmin <- min(IntronModel$intron_start)
		zZmax <- max(IntronModel$intron_end)
		DataModel.selected <- DataModel.selected[start(DataModel.selected) > zZmin]
		DataModel.selected <- DataModel.selected[end(DataModel.selected) < zZmax]
		anno.tr <- AnnotationTrack(range = Final.GR, name = Name.IntronModel, 
			feature = IntronModel$Feature, group = IntronModel$transcript_name, 
			id = IntronModel$ID, genome = genome, JOI = "pink", fun=ir_anno, detailsFunArgs = setNames(list(DataModel.selected, experiment, name_by), c("data", "experiment", "name_by")))
		trackcount = trackcount + 1
		tracks[[trackcount]] = anno.tr
	}

	if (!is.null(bamfile)) {
		coverage = DataTrack(range = bamfile, genome = genome, type = "l", name = "Coverage", window = -1)
		trackcount = trackcount + 1
		tracks[[trackcount]] = coverage
	}

	if(JnChr != 0) {
		axisTrack <- GenomeAxisTrack(name = JnChr)
	} else {
		axisTrack <- GenomeAxisTrack()
	}
	tracks[[1]] = axisTrack
	
	if(Zmax > Zmin) {
		plotTracks(tracks, showFeatureId=TRUE,showID=TRUE,groupAnnotation="group", horizontal = FALSE,
			from = Zmin, to = Zmax,
			fontcolor.feature = 1, cex.feature=0.6, just.group = "below", collapse = FALSE)
	
	} else {
	plotTracks(tracks, showFeatureId=TRUE,showID=TRUE,groupAnnotation="group", horizontal = FALSE,
		fontcolor.feature = 1, cex.feature=0.6, just.group = "below", collapse = FALSE, sizes=c(1,3,2), details.size=0.75)
	}
}
