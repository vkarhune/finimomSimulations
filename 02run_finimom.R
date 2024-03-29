# rm(list=ls())

args <- commandArgs(trailingOnly = T)

generegion <- args[1]
insampleLD <- as.logical(args[2])
chainlength <- as.integer(args[3])
burninprop <- as.numeric(args[4])
tau <- as.numeric(args[5])
maxsize <- as.integer(args[6])
simnumstart <- as.integer(args[7])
u <- as.numeric(args[8]) # set to zero for varying b0
stdize <- as.logical(args[9])
applyclumping <- as.logical(args[10])
ala <- as.logical(args[11])

niter <- chainlength/(1 - burninprop)

if(0){
 generegion <- "CSNK1A1L"
 insampleLD <- FALSE
 chainlength <- 50000
 burninprop <- 0.2
 tau <- 0.00385
 maxsize <- 10
 simnumstart <- 501
 u <- 2.25
 stdize <- TRUE
 applyclumping <- TRUE
 ala <- TRUE
 
 niter <- chainlength/(1 - burninprop)
}

if(u == 0){ u <- NULL }# NOTE: this changed later
if(is.null(u)){ out_u <- 0
 } else {
 out_u <- u
}

# R-4.2.1
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

# there were some flipped alleles in the reference LD - fix these
if(0){
 LD2 <- readRDS(paste0("data/LDmat_", paste0(c(generegion, samplesize, clump), collapse = "_"), ".Rds"))
 plot(LDmat[lower.tri(LDmat)],LD2[[1]][lower.tri(LD2[[1]])])
 
 slopes <- sapply(seq_len(nrow(LDmat)), function(i){
  m1 <- lm(LD2[[1]][i,] ~ LDmat[i,])
  coef(m1)[2]
 })
 
 which(slopes < 0)
 # LDmat[i, ] LDmat[i, ]
 #      1628       1631
}

if(!(insampleLD) & generegion %in% "MDGA2"){
 LDmat[c(490, 492, 1792, 1832, 1834),] <- -LDmat[c(490, 492, 1792, 1832, 1834),]
 LDmat[,c(490, 492, 1792, 1832, 1834)] <- -LDmat[,c(490, 492, 1792, 1832, 1834)]
  
 stopifnot(all(diag(LDmat) == 1))
}

if(!(insampleLD) & generegion %in% "GCNT2"){
 LDmat[c(757, 760, 767),] <- -LDmat[c(757, 760, 767),]
 LDmat[,c(757, 760, 767)] <- -LDmat[,c(757, 760, 767)]

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
 p <- nrow(summarystats)

  
  
 d <- data.frame(
  pip = out$pip,

  y = ifelse(seq_len(p) %in% truecausals, 1, 0)
 )  
 
 cat(sprintf("%s\n", file))
 
 return(list(d, out$sets, table(out$samples[[2]]), elapsed))
 
})

reffile <- ifelse(insampleLD, "insampleLD", "refLD")

outfilename <- paste0(sub(".Rds", "",
  paste0(c("results/finimom", reffile, filenamesplit[[1]][2:6], tau, out_u, stdize),
   collapse = "_")),
 ".Rds")

saveRDS(reslist, file = outfilename)



Sys.Date()

sessionInfo()
