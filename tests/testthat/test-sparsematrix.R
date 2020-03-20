context("Test sparse/HDF5 backends")

sampler <- c(0,0,0,0,0,0,1,2,2,2,4)
num_genes <- 50
num_cells <- 20 
gene_names=paste0("gene", 1:num_genes)
cell_names=paste0("cell", 1:num_cells)

condition <- factor(c(rep("Con",10), rep("test", 10)))
raw_counts <- replicate(num_cells, sample(sampler, num_genes, replace=TRUE))
log_counts <- log2( raw_counts + 1)

colnames(raw_counts) <-cell_names
rownames(raw_counts) <-gene_names
colnames(log_counts) <-cell_names
rownames(log_counts) <-gene_names

test_that('Can make SingleCellAssay from dense SCE', {
# With standard dense matrix
sce        <- SingleCellExperiment(assays=SimpleList(counts=raw_counts, logcounts=log_counts), colData=condition)

sca = SceToSingleCellAssay(sce)
zz_dense <<- zlm( ~X, sca=sca) 
})

test_that('Can make SCA from sparse SCE', {
    log_counts.sm <- Matrix::Matrix(log_counts, sparse = TRUE) 
    counts.sm = Matrix::Matrix(raw_counts, sparse = TRUE)
    sce.sm        <- SingleCellExperiment(
        assays  = SimpleList(counts=counts.sm, logcounts=log_counts.sm),
        colData = condition)
    sca.sm = SceToSingleCellAssay(sce.sm)
    ## Should check memory usage here or something?
    zz_sparse <<- zlm( ~X, sca=sca.sm) 
})

test_that('Results are equal', {
    sdense = summary(zz_dense)$datatable
    ssparse = summary(zz_sparse)$datatable
    ## For some reason, the table is sorted by contrast and then `z` statistic
    ## There are ties that end up with arbitrary sort order.
    m = merge(sdense,  ssparse, by = c('primerid', 'component', 'contrast'))
    expect_equal(m$z.x, m$z.y)
})
