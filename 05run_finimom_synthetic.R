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
applyclumping <- as.logical(args[10])
ala <- as.logical(args[11])

niter <- chainlength/(1 - burninprop)

if(0){
 generegion <- "500variants"
 insampleLD <- TRUE
 chainlength <- 50000
 burninprop <- 0.2
 tau <- 0.000416
 maxsize <- 10
 simnumstart <- 501
 u <- 2
 stdize <- TRUE
 applyclumping <- TRUE
 ala <- TRUE
 
 niter <- chainlength/(1 - burninprop)
}

if(u == 0){ u <- NULL }
if(is.null(u)){ out_u <- 0
 } else {
 out_u <- u
}


# .libPaths("/home/vkarhune/rpackages/")

# R-4.2.1
.libPaths(paste0("/projappl/minmanni/project_rpackages_", gsub("\\.", "", getRversion())))

library(finimom) # remotes::install_github("vkarhune/finimom")


 
filename <- list.files("data", pattern = paste0("simdatasynthetic", simnumstart, "_", generegion, ".*.Rds"), full.names = T)

filenamesplit <- strsplit(filename, "_")

generegion <- filenamesplit[[1]][2]
samplesize <- filenamesplit[[1]][3]
clump <- filenamesplit[[1]][4]


LDfile <- ifelse(insampleLD,
 paste0("data/LDmatsynthetic_", paste0(c(sub("variants", "", generegion), samplesize, "0"), collapse = "_"), ".Rds"),
 paste0("data/refLDmatsynthetic_", paste0(c(sub("variants", "", generegion), samplesize, "0"), collapse = "_"), ".Rds")
)

LD <- readRDS(LDfile)

LDmat <- LD[[1]]
eafs <- LD[[2]]



# fix flipped alleles in reference LD
if(0){
 LD2 <- readRDS(paste0("data/LDmatsynthetic_", paste0(c(sub("variants", "", generegion), samplesize, "0"), collapse = "_"), ".Rds"))
 plot(LDmat[lower.tri(LDmat)],LD2[[1]][lower.tri(LD2[[1]])])
 
 slopes <- sapply(seq_len(nrow(LDmat)), function(i){
  m1 <- lm(LD2[[1]][i,-i] ~ LDmat[i,-i]) # update: ith element deleted
  coef(m1)[2]
 })
 
 which(slopes < 0)

}

if(!(insampleLD) & generegion %in% "2000variants"){
 LDmat[1493,] <- -LDmat[1493,]
 LDmat[,1493] <- -LDmat[,1493]
 
 stopifnot(all(diag(LDmat) == 1))
}

if(!(insampleLD) & generegion %in% "5000variants"){
 LDmat[c(385, 542, 587, 631, 812, 1091, 2353, 2464, 4126, 4806),] <- -LDmat[c(385, 542, 587, 631, 812, 1091, 2353, 2464, 4126, 4806),]
 LDmat[,c(385, 542, 587, 631, 812, 1091, 2353, 2464, 4126, 4806)] <- -LDmat[,c(385, 542, 587, 631, 812, 1091, 2353, 2464, 4126, 4806)]
 
 stopifnot(all(diag(LDmat) == 1))
}

if(!(insampleLD) & generegion %in% "10000variants"){
 swapinds <- c(1504, 2364, 3483, 3559, 3564, 3568, 3569, 3575, 3577, 3592,
  3594, 3599, 3602, 3605, 3606, 3607, 3608, 3609, 3614,
  3616, 3622, 3624, 3627, 3628, 3629, 3638, 4452, 5083,
  5085, 6292, 6732, 7021, 7025, 7027, 7031, 7042, 7056,
  8403, 8995, 8998, 9101)

 LDmat[swapinds,] <- -LDmat[swapinds,]
 LDmat[,swapinds] <- -LDmat[,swapinds]
 
 stopifnot(all(diag(LDmat) == 1))
}



N <- as.numeric(gsub("\\D", "", samplesize))



files <- list.files("data", pattern = paste0("simdata.*._", paste0(filenamesplit[[1]][2:6], collapse = "_")), full.names = T)

reslist <- lapply(files, function(file){
# reslist <- lapply(files[1:6], function(file){
# file <- files[1]
# niter <- 12500

 dat <- readRDS(file)

 summarystats <- dat[[1]]
 
 if(is.null(u)){ u <- log(nrow(summarystats) - 1)/log(nrow(summarystats)) }
 if(u == -999){ u <- log(999)/log(nrow(summarystats)) }
 
prc <- proc.time()
 out <- finimom(beta = summarystats[,1], se = summarystats[,2],
  eaf = eafs, R = LDmat,
  maxsize = maxsize, tau0 = tau, r0 = 1,
  niter = niter, burnin = niter - chainlength,
  excl.burnin = TRUE, # p = nrow(summarystats),
  standardize = stdize, clump = applyclumping, clump_r2 = 0.99^2,
  insampleLD = insampleLD,
  check_ld = !(insampleLD),
  a0 = 1, u = u,
  ala = ala,
  verbose = TRUE,
  cs = TRUE,
  pip = TRUE
 )
elapsed <- (proc.time() - prc)[[3]]

 truecausals <- dat[[2]]
 # p <- ncol(out[[1]]) FIX
 p <- nrow(summarystats)


  
 if(0){
  library(susieR)
  z <- summarystats[,1]/summarystats[,2]
  prc <- proc.time()
  fitted_rss <- susie_rss(z, LDmat, L = maxsize, n = N)
  proc.time() - prc
 }
  
 d <- data.frame(
  pip = out$pip,

  y = ifelse(seq_len(p) %in% truecausals, 1, 0)
 )  
 
 cat(sprintf("%s\n", file))
 
 return(list(d, out$sets, table(out$samples[[2]]), elapsed))
 
})

reffile <- ifelse(insampleLD, "insampleLD", "refLD")

outfilename <- paste0(sub(".Rds", "",
  paste0(c("results/finimomsynthetic", reffile, filenamesplit[[1]][2:6], tau, out_u, stdize),
   collapse = "_")),
 ".Rds")

saveRDS(reslist, file = outfilename)



Sys.Date()

sessionInfo()
