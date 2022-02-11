library(affy)
library(affyPLM)
library(sva)
library(AnnotationDbi)
library(hgu133plus2.db)

filepath = '/project/bf528/project_1/data/GSE39582/CEL_files'

celpath <- system.file('celfiles', package='affydata')
fns <- list.celfiles(path=celpath, full.names=TRUE)

# read in files
celbatch <- ReadAffy(celfile.path=filepath)

# normalize files together
rma(celbatch)

# convert AffyBatch to PLMset
pset <- fitPLM(celbatch, normalize = TRUE, background = TRUE)

# relative log expression
rle_stats <- data.frame(t(affyPLM::RLE(pset, type='stats')))

# normalized unscaled standard errors (NUSE)
nuse_stats <- data.frame(t(NUSE(pset, type = 'stats')))
