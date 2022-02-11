library(affy)
library(affyPLM)
library(sva)
library(AnnotationDbi)
library(hgu133plus2.db)
library(ggplot2)

filepath = '/project/bf528/project_1/data/GSE39582/CEL_files'

celpath <- system.file('celfiles', package='affydata')
fns <- list.celfiles(path=celpath, full.names=TRUE)

# read in files
celbatch <- ReadAffy(celfile.path=filepath)

# normalize files together
rma(celbatch)

# convert AffyBatch to PLMset
pset <- fitPLM(celbatch, normalize = TRUE, background = TRUE)


# relative log expression (RLE)
rle_stats <- data.frame(t(affyPLM::RLE(pset, type='stats')))

# plot rle_stats
rle_medians <- ggplot(rle_stats, aes(x=median)) + 
  geom_histogram(bins=50, 
                 color = 'dodgerblue4', 
                 fill = 'white') +
  labs(title = 'RLE Medians') +
  theme_classic()
rle_medians


# normalized unscaled standard errors (NUSE)
nuse_stats <- data.frame(t(NUSE(pset, type = 'stats')))

# plot nuse_stats
nuse_medians <- ggplot(nuse_stats, aes(x=median)) + 
  geom_histogram(bins=50, 
                 color = 'dodgerblue4', 
                 fill = 'white') +
  labs(title = 'NUSE Medians') +
  theme_classic()
nuse_medians


# correction for batch effects
annotation_filepath = '/project/bf528/project_1/doc/proj_metadata.csv'

