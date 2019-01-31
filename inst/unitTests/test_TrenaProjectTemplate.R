library(TrenaProjectTemplate)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tProj")) {
   message(sprintf("--- creating instance of TrenaProjectTemplate"))
   tProj <- TrenaProjectTemplate();
   }
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()
   test_supportedGenes()
   test_variants()
   test_footprintDatabases()
   test_expressionMatrices()
   test_setTargetGene()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   message(sprintf("--- test_constructor"))

   checkTrue(all(c("TrenaProjectTemplate", "TrenaProject") %in% is(tProj)))

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_supportedGenes <- function()
{
   message(sprintf("--- test_supportedGenes"))

   subset.expected <- c("Abca1")
   checkTrue(all(subset.expected %in% getSupportedGenes(tProj)))

} # test_supportedGenes
#------------------------------------------------------------------------------------------------------------------------
test_variants <- function()
{
   message(sprintf("--- test_variants"))

   checkEquals(getVariantDatasetNames(tProj), character(0))

} # test_variants
#------------------------------------------------------------------------------------------------------------------------
test_footprintDatabases <- function()
{
   message(sprintf("--- test_footprintDatabases"))

   expected <- c("extraembryonic_structure_wellington_16", "extraembryonic_structure_wellington_20",
                 "extraembryonic_structure_hint_16", "extraembryonic_structure_hint_20")
   checkTrue(is.na(getFootprintDatabaseNames(tProj)))
   checkTrue(is.na(getFootprintDatabaseHost(tProj)))

} # test_footprintDatabases
#------------------------------------------------------------------------------------------------------------------------
test_expressionMatrices <- function()
{
   expected <- c("thioglycollate-elicited-peritoneal-macrophages")
   checkTrue(all(expected %in% getExpressionMatrixNames(tProj)))

   mtx <- getExpressionMatrix(tProj, expected[1])
   checkEquals(dim(mtx), c(6890, 14))

} # test_expressionMatrices
#------------------------------------------------------------------------------------------------------------------------
test_setTargetGene <- function()
{
   message(sprintf("--- test_setTargetGene"))

   setTargetGene(tProj, "Abca1")
   checkEquals(getTargetGene(tProj), "PGF")

   message(sprintf("    transcripts"))
   tbl.transcripts <- getTranscriptsTable(tProj)
   checkTrue(nrow(tbl.transcripts) >= 9)

   message(sprintf("    enhancers"))
   tbl.enhancers <- getEnhancers(tProj)
   checkEquals(colnames(tbl.enhancers), c("chrom", "start", "end", "type", "combinedScore", "geneSymbol"))
   checkTrue(nrow(tbl.enhancers) >= 28)

   message(sprintf("    encode DHS"))
   tbl.dhs <- getEncodeDHS(tProj)
   checkTrue(nrow(tbl.dhs) > 1000)
   checkEquals(colnames(tbl.dhs), c("chrom", "chromStart", "chromEnd", "count", "score"))

   chromosome <- unique(c(tbl.dhs$chrom, tbl.enhancers$chrom))
   checkEquals(length(chromosome), 1)
   start <- min(tbl.dhs$chromStart)
   end   <- max(tbl.dhs$chromEnd)

   chromosome <- tbl.transcripts$chr[1]
   start <- tbl.transcripts$start[1]
   end <- start + 500

   message(sprintf("    ChIP-seq"))
   tbl.chipSeq <- getChipSeq(tProj, chrom=chromosome, start=start, end=end, tfs=NA)
   checkTrue(nrow(tbl.chipSeq) > 40)
   checkEquals(colnames(tbl.chipSeq), c("chrom", "start", "endpos", "tf", "name", "strand", "peakStart", "peakEnd"))

   tbl.chipSeq.myc <- getChipSeq(tProj, chrom=chromosome, start=start, end=end, tfs="MYC")
   checkTrue(nrow(tbl.chipSeq.myc) >= 4)
   checkTrue(nrow(tbl.chipSeq.myc) < 10)

} # test_setTargetGene
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
