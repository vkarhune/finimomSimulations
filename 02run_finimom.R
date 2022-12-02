# rm(list=ls())

args <- commandArgs(trailingOnly = T)

generegion <- args[1]
insampleLD <- as.logical(args[2])
chainlength <- as.integer(args[3])
burninprop <- as.numeric(args[4])
tau <- as.numeric(args[5])
maxsize <- as.integer(args[6])
simnumstart <- as.integer(args[7])
u <- as.numeric(args[8])
stdize <- as.logical(args[9])

niter <- chainlength/(1 - burninprop)

if(0){
 generegion <- "MDGA2"
 chainlength <- 10000
 burninprop <- 0.2
 tau <- 0.0083
 maxsize <- 10
 niter <- chainlength/(1 - burninprop)
 simnumstart <- 501
 insampleLD <- FALSE
 u <- 1.5
 stdize <- TRUE
}

# R-4.1.1
.libPaths(paste0("/projappl/minmanni/project_rpackages_", gsub("\\.", "", getRversion())))

library(finimom) # remotes::install_github("vkarhune/finimom")


 
filename <- list.files("data", pattern = paste0("simdata", simnumstart, "_", generegion, ".*.Rds"), full.names = T)

filenamesplit <- strsplit(filename, "_")

generegion <- filenamesplit[[1]][2]
samplesize <- filenamesplit[[1]][3]
clump <- filenamesplit[[1]][4]


LDfile <- ifelse(insampleLD,
 paste0("data/LDmat_", paste0(c(generegion, samplesize, clump), collapse = "_"), ".Rds"),
 paste0("data/LDmat_NFBC1986_", generegion, ".Rds")
)

LD <- readRDS(LDfile)

LDmat <- LD[[1]]
eafs <- LD[[2]]

N <- as.numeric(gsub("\\D", "", samplesize))



files <- list.files("data", pattern = paste0("simdata.*._", paste0(filenamesplit[[1]][2:6], collapse = "_")), full.names = T)

reslist <- lapply(files, function(file){
# reslist <- lapply(files[1:6], function(file){
# file <- files[1]
# niter <- 12500

 dat <- readRDS(file)

 summarystats <- dat[[1]]

 out <- posterior_samples(beta = summarystats[,1], se = summarystats[,2],
                  eaf = eafs, R = LDmat,
                  maxsize = maxsize, tau0 = tau, r0 = 1, niter = niter,
                  burnin = 0, excl.burnin = F, p = nrow(summarystats),
                  standardize = stdize, clump = T, clump_r2 = 0.99^2,
                  check_ld = !(insampleLD),
                  a0 = 1, b0 = nrow(summarystats)^u,
                  verbose = F)



 outshort <- extract_subchain(out, burnin = chainlength*burninprop, niter = chainlength)

 truecausals <- dat[[2]]
 p <- nrow(summarystats)

 num_signals_short <- names(which.max(table(outshort[[2]])))
 sets_short <- get_credible_sets(samples = outshort, num_signals = num_signals_short, level = 0.95)
 pips_short <- get_pips(samples = outshort)

 num_signals_best <- outshort[[2]][match(names(sort(table(outshort[[3]]), decreasing = T)[1]), outshort[[3]])]
 sets_best <- get_credible_sets(samples = outshort, num_signals = num_signals_best, level = 0.95)

 sets_bf <- list(get_csbf(beta = summarystats[,1]*sqrt(2*eafs*(1-eafs)), se = summarystats[,2]*sqrt(2*eafs*(1-eafs)),
  tau = tau, r = 1, level = 0.95))
   
 d <- data.frame(
  pip = pips_short[,2],
  y = ifelse(seq_len(p) %in% truecausals, 1, 0)
 )  
 
 cat(sprintf("%s\n", file))
 
 return(list(d, sets_short, sets_best, sets_bf, table(outshort[[2]])))
 
})

reffile <- ifelse(insampleLD, "insampleLD", "refLD")

saveRDS(reslist,
 file = paste0(sub(".Rds", "",
         paste0(c("results/finimom", reffile, filenamesplit[[1]][2:6], u, stdize),
          collapse = "_")),
         ".Rds")
)



Sys.Date()

sessionInfo()
