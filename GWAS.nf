// Default
params.ncpus     = "auto"
params.threshold = "0.05"
params.impute    = true
params.outdir    = "result"
params.plot_type = "jpg"


// Parse the input parameters
genotype  = file(params.genotype)
phenotype = file(params.phenotype)
methods   = Channel.from(params.methods.split(","))
ncpus     = params.ncpus != "auto" ?: "parallel::detectCores() - 1" 
threshold = params.threshold
impute    = params.impute ? "TRUE" : "FALSE"
outdir    = params.outdir
plot_type = params.plot_type


// Define script Path
SCRIPT_DIR = new File(workflow.projectDir.toString(), "bin")
report     = new File(SCRIPT_DIR, "Plot.R")

// Log
log.info """\
=============================================
G W A S - N F v0.1
=============================================
genotype  : $genotype
phenotype : $phenotype
methods   : $methods
ncpus     : $ncpus
threshold : $threshold
outdir    : $outdir
"""


// PART 1: Preprocess & Check ===================

process DataFormatConvert {
    input:
        file geno from genotype
        file pheno from phenotype

    output:
        file "mvp.geno.desc" into mvp_geno_desc
        file "mvp.geno.bin" into mvp_geno_bin
        file "mvp.geno.map" into mvp_geno_map
        file "mvp.phe" into mvp_phe_plot
        file "*.pheno.txt" into mvp_phe mode flatten

    """
    #!/usr/bin/env Rscript
    library(rMVP)
    MVP.Data(filePhe = "$pheno", fileVCF = "$geno", ncpus=$ncpus, SNP.impute = NULL, out = "mvp")
    pheno <- read.table("mvp.phe", header = TRUE)
    for (i in 2:ncol(pheno)) {
        write.table(pheno[, c(1, i)], file = paste0(colnames(pheno)[i], ".pheno.txt"), quote = FALSE, row.names = FALSE)
    }
    if ($impute) {
        MVP.Data.impute(mvp_prefix = "mvp", ncpus = $ncpus)
    } else {
        bigmat <- attach.big.matrix("mvp.geno.desc")
        if (rMVP:::hasNA(bigmat@address))
            stop("Error: NA in genotype")
    }
    """
}

analysis_plain = mvp_phe.combine(methods)


// PART 2: Run analysis ==========================

process RunAnalysis {
    publishDir outdir

    input:
        file geno_desc from mvp_geno_desc
        file geno_bin from mvp_geno_bin
        file geno_map from mvp_geno_map
        set (file(phe), method) from analysis_plain

    output:
        file "*.csv" into pmap_filter, pmap_plot_qq, pmap_plot_manhattan

    """
    #!/usr/bin/env Rscript
    library(rMVP)
    geno  <- attach.big.matrix("$geno_desc")
    pheno <- read.table("$phe", header = TRUE)
    map   <- read.table("$geno_map", header = TRUE)
    imMVP <- MVP(phe=pheno, geno=geno, map=map, method=c("$method"),
        ncpus=$ncpus, file.output=c("pmap"))
    """
}


// PART 3: Postprocess & Plot ====================

process FilterSignalSNP {
    publishDir outdir

    input:
        file pmap from pmap_filter
    
    output:
        file "*_signals.csv" into mvp_signals_pmap

    """
    #!/usr/bin/env Rscript
    pmap     <- read.csv("$pmap", header = TRUE)
    out_snps <- which(pmap[, ncol(pmap)] < ($threshold / nrow(pmap)))
    out_pmap <- pmap[out_snps, ]
    out_file <- paste0(substr("$pmap", 1, nchar("$pmap") - 4), "_signals.csv")
    write.csv(out_pmap, file = out_file, quote = FALSE, row.names = FALSE)
    """
}


process PlotSnpDensity {
    publishDir outdir

    input:
        file geno_map from mvp_geno_map
    
    output:
        file "*.$plot_type"

    """
    #!/usr/bin/env Rscript
    source("$report")
    map <- read.table("$geno_map", header = TRUE, stringsAsFactors = FALSE)
    drawDensityPlot(map, file.type = "$plot_type")
    """
}


process PlotPhenotypeDistribution {
    publishDir outdir

    input:
        file phe from mvp_phe_plot

    output:
        file "*.$plot_type"

    """
    #!/usr/bin/env Rscript
    source("$report")
    phe <- read.table("$phe", header = TRUE)
    for(i in 2:ncol(phe)) {
        drawPhenotypeHistogram(phe[, i], colnames(phe)[i], file.type = "$plot_type")
    }
    """
}


process PlotQQ {
    publishDir outdir

    input:
        file pmap from pmap_plot_qq

    output:
        file "*.$plot_type"

    """
    #!/usr/bin/env Rscript
    source("$report")
    pmap <- read.csv("$pmap", header = TRUE)
    name <- colnames(pmap)[ncol(pmap)]
    pval <- pmap[, ncol(pmap)]
    drawQQPlot(pval, name, file.type = "$plot_type")
    """
}


process PlotManhattan {
    publishDir outdir

    input:
        file pmap from pmap_plot_manhattan

    output:
        file "*.$plot_type"

    """
    #!/usr/bin/env Rscript
    source("$report")
    pmap <- read.csv("$pmap", header = TRUE)
    name <- colnames(pmap)[ncol(pmap)]
    drawManhattanPlot(pmap[, c(2:3, ncol(pmap))], name, density = FALSE, threshold = $threshold, file.type = "$plot_type")
    """
}