library("AnnotationHub")
library("BSgenome.Hsapiens.UCSC.hg19")
library("JASPAR2014")
library("TFBSTools")

# function to turn matches into GRanges
gr_of_matches <- function(gr, matches, genome){
  if (length(matches) > 0){
    GRanges(seqnames = seqnames(gr),
            ranges = IRanges(start = start(matches@views) + start(gr),
                             end = end(matches@views) + start(gr)),
            strand = matches@strand,
            score = matches@score,
            relscore = relScore(matches),
            seqinfo = seqinfo(genome))
  }
} 

# get CTCF motif PWM
opts <- list(name="CTCF", species = "9606") #species = human
pfm <- getMatrixSet(JASPAR2014, opts)[[1]]
pwm <- TFBSTools:::toPWM(pfm)

# get CTCF peaks in IMR90 cells
hub <- AnnotationHub()
ctcf_imr90 <- query(hub, c("ENCODE", "CTCF", "IMR90"))[[1]] # Uniformly processed peaks 

# Find position of best CTCF motif match in each peak
ctcf_imr90_seqs <- getSeq(Hsapiens, ctcf_imr90)
imr90_matches <- lapply(ctcf_imr90_seqs, function(s){
  searchSeq(pwm, s, min.score = "75%")}) # 80% is default cutoff, reduce to increase % of peaks with a match

imr90_gr_list <- Map(function(g,m){gr_of_matches(g,m, Hsapiens)}, 
                     as.list(ctcf_imr90), imr90_matches)

imr90_gr_list <- imr90_gr_list[sapply(imr90_gr_list, length) > 0]
imr90_gr_list <- lapply(imr90_gr_list, function(g){
  g[which.max(g$score)]
})

all_gr <- do.call("c", imr90_gr_list)

# split by strand and export bed files
minus_gr <- all_gr[strand(all_gr)=="-"]
plus_gr <- all_gr[strand(all_gr)=="+"]

export.bed(object = minus_gr, con = "imr90_ctcf_minus.bed")
export.bed(object = plus_gr, con = "imr90_ctcf_plus.bed")

sessionInfo()
