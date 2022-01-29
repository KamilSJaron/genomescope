#!/usr/bin/env Rscript

## GenomeScope: Fast Genome Analysis from Unassembled Short Reads
##
## This is the automated script for computing genome characteristics
## from a kmer histogram file, k-mer size, and readlength

# all used funcations can be found in this file
source('R/regular_genomescope_functions.R')

## Number of rounds before giving up
NUM_ROUNDS=4

## Coverage steps to trim off between rounds
START_SHIFT=5

## Typical cutoff for sequencing error
TYPICAL_ERROR = 15

## Max rounds on NLS
MAX_ITERATIONS=20

## Overrule if two scores are within this percent (0.05 = 5%) but larger difference in het
SCORE_CLOSE = 0.20

## Overrule heterozygosity if there is a large difference in het rate
SCORE_HET_FOLD_DIFFERENCE = 10

## Print out VERBOSEging messages (0/1)
VERBOSE = 0

## Suppress the warnings if the modeling goes crazy, those are in try/catch blocks anyways
options(warn=-1)

## Colors for plots
COLOR_BGCOLOR  = "light grey"
COLOR_HIST     = "#56B4E9"
COLOR_4PEAK    = "black"
COLOR_2PEAK    = "#F0E442"
COLOR_ERRORS   = "#D55E00"
COLOR_KMERPEAK = "black"
COLOR_RESIDUAL = "purple"
COLOR_COVTHRES = "red"

## Main program starts here
###############################################################################

args <- commandArgs(TRUE)

if(length(args) < 4) {
	cat("USAGE: genomescope.R histogram_file k-mer_length read_length output_dir [kmer_max] [verbose]\n")
} else{

    ## Load the arguments from the user
	histfile   <- args[[1]]
	k          <- as.numeric(args[[2]])
	readlength <- as.numeric(args[[3]])
	foldername <- args[[4]]

    maxCovGenomeLen = -1

    if ((length(args) >= 5)) {
        maxCovGenomeLen = as.numeric(args[[5]])
    }

    if ((length(args) == 6) && (as.numeric(args[[6]] == 1))) { VERBOSE = 1 }

    ## values for testing
    #histfile <- "~/build/genomescope/simulation/simulation_results/Arabidopsis_thaliana.TAIR10.26.dna_sm.toplevel.fa_het0.01_br1_rl100_cov100_err0.01_reads.fa21.hist"
    #k <- 21
    #readlength <- 100
    #foldername <- "~/build/genomescope/simulation/simulation_analysis/Arabidopsis_thaliana.TAIR10.26.dna_sm.toplevel.fa_het0.01_br1_rl100_cov100_err0.01_reads.fa21.hist"

    if (k > readlength) { stop("K cannot be greater than readlength") }

    cat(paste("GenomeScope analyzing ", histfile, " k=", k, " readlen=", readlength, " outdir=", foldername, "\n", sep=""))

	dir.create(foldername, showWarnings=FALSE)

	kmer_prof <- read.table(file=histfile, header=FALSE)

    minkmerx = 1;
    if (kmer_prof[1,1] == 0) {
        if (VERBOSE) { cat("Histogram starts with zero, reseting minkmerx\n");  }
        minkmerx = 2;
    }

	kmer_prof <- kmer_prof[c(minkmerx:(length(kmer_prof[,2])-1)),] #get rid of the last position
    kmer_prof_orig <- kmer_prof

    ## Initialize the status
    progressFilename <- paste(foldername,"/progress.txt",sep="")
	cat("starting", file=progressFilename, sep="\n")

    ## try to find the local minimum between errors and the first (heterozygous) peak
    start <- which(kmer_prof[,2]==min(kmer_prof[1:TYPICAL_ERROR,2]))

    maxCovIndex = -1

    ## Figure out which kmers to exclude, if any
    if(maxCovGenomeLen == -1){
        maxCovIndex <- length(kmer_prof[,1])
    }
    else
    {
        ## Figure out the index we should use for this coverage length
        x <- kmer_prof[,1]
        maxCovIndex <- length(x[x<=maxCovGenomeLen])
    }

    if (VERBOSE) { cat(paste("using maxCovGenomeLen:", maxCovGenomeLen, " with index:", maxCovIndex, "trying 4peak model... \n")) }

    ## terminate after NUM_ROUND iterations, store best result so far in container
	round <- 0
	best_container <- list(NULL,0)

	while(round < NUM_ROUNDS)
    {
        cat(paste("round", round, "trimming to", start, "trying 4peak model... "), file=progressFilename, sep="", append=TRUE)
        if (VERBOSE) { cat(paste("round", round, "trimming to", start, "trying 4peak model... \n")) }

        ## Reset the input trimming off low frequency error kmers
        kmer_prof=kmer_prof_orig[1:maxCovIndex,]
        x <- kmer_prof[start:maxCovIndex,1]
        y <- kmer_prof[start:maxCovIndex,2]

        model_4peaks <- estimate_Genome_4peak2(kmer_prof, x, y, k, readlength, round, foldername)

        if (!is.null(model_4peaks[[1]])) {
          cat(paste("converged. score: ", model_4peaks[[2]]$all[[1]]), file=progressFilename, sep="\n", append=TRUE)

          if (VERBOSE)
          {
            mdir = paste(foldername, "/round", round, sep="")
	        dir.create(mdir, showWarnings=FALSE)
	        report_results(kmer_prof,kmer_prof_orig, k, model_4peaks, mdir)
          }
        } else {
          cat(paste("unconverged"), file=progressFilename, sep="\n", append=TRUE)
        }

		#check if this result is better than previous
        if (!is.null(model_4peaks[[1]]))
        {
          if (is.null(best_container[[1]]))
          {
            if (VERBOSE) { cat(paste("no previous best, updating best")) }
            best_container = model_4peaks
          }
          else
          {
            pdiff = abs(model_4peaks[[2]]$all[[1]] - best_container[[2]]$all[[1]]) / max(model_4peaks[[2]]$all[[1]], best_container[[2]]$all[[1]])

            if (pdiff < SCORE_CLOSE)
            {
              hetm = summary(model_4peaks[[1]])$coefficients['r',][[1]]
              hetb = summary(best_container[[1]])$coefficients['r',][[1]]

              if (hetb * SCORE_HET_FOLD_DIFFERENCE < hetm)
              {
                if (VERBOSE) { cat(paste("model has significantly higher heterozygosity but similar score, overruling")) }
                best_container = model_4peaks
              }
              else if (hetm * SCORE_HET_FOLD_DIFFERENCE < hetb)
              {
                if (VERBOSE) { cat(paste("previous best has significantly higher heterozygosity and similar score, keeping")) }
              }
              else if (model_4peaks[[2]]$all[[1]] < best_container[[2]]$all[[1]])
              {
                if (VERBOSE) { cat(paste("score is marginally better but het rate is not extremely different, upating")) }
                best_container = model_4peaks
              }
            }
            else if (model_4peaks[[2]]$all[[1]] < best_container[[2]]$all[[1]])
            {
              if (VERBOSE) { cat(paste("score is significantly better, upating")) }
              best_container = model_4peaks
            }
          }
        }

        ## Ignore a larger number of kmers as errors
        start <- start + START_SHIFT
		round <- round + 1
	}
    ## Report the results, note using the original full profile
	report_results(kmer_prof,kmer_prof_orig, k, best_container, foldername)
}
