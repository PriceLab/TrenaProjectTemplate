#----------------------------------------------------------------------------------------------------
#' @import methods
#' @import TrenaProject
#' @importFrom AnnotationDbi select
#' @import org.Hs.eg.db
#'
#' @title TrenaProjectTemplate-class
#'
#' @name TrenaProjectTemplate-class
#' @rdname TrenaProjectTemplate-class
#' @aliases TrenaProjectTemplate
#' @exportClass TrenaProjectTemplate
#'

.TrenaProjectTemplate <- setClass("TrenaProjectTemplate",
                                  contains="TrenaProject")

#----------------------------------------------------------------------------------------------------
#' Define an object of class TrenaProjectTemplate
#'
#' @description
#' Expression, variant and covariate data for the genes of interest (perhaps unbounded) for pre-term birth studies
#'
#' @rdname TrenaProjectTemplate-class
#'
#' @export
#'
#' @return An object of the TrenaProjectTemplate class
#'

TrenaProjectTemplate <- function(quiet=TRUE)

{
   genomeName <- "mm10"

   directory <- system.file(package="TrenaProjectTemplate", "extdata", "geneSets")
   geneSet.files <- list.files(directory)
   geneSets <- list()
   for(file in geneSet.files){
      full.path <- file.path(directory, file)
      genes <- scan(full.path, sep="\t", what=character(0), quiet=TRUE)
      geneSet.name <- sub(".txt", "", file)
      geneSets[[geneSet.name]] <- genes
      }

   footprintDatabaseNames <- NA_character_;
   expressionDirectory <- system.file(package="TrenaProjectTemplate", "extdata", "expression")
   variantsDirectory <- system.file(package="TrenaProjectTemplate", "extdata", "variants")
   footprintDatabaseHost <- NA_character_;

   covariatesFile <- NA_character_;

   stopifnot(file.exists(expressionDirectory))
   #stopifnot(file.exists(variantsDirectory))
   # stopifnot(file.exists(covariatesFile))

   .TrenaProjectTemplate(TrenaProject(supportedGenes=geneSets[[1]],
                                      genomeName=genomeName,
                                      footprintDatabaseHost=footprintDatabaseHost,
                                      footprintDatabaseNames=footprintDatabaseNames,
                                      expressionDirectory=expressionDirectory,
                                      variantsDirectory=variantsDirectory,
                                      covariatesFile=covariatesFile,
                                      quiet=quiet
                                      ))

} # TrenaProjectTemplate, the constructor
#----------------------------------------------------------------------------------------------------
