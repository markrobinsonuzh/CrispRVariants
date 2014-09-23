library(Rsamtools)
library(GenomicAlignments)
param <- ScanBamParam(what = c("seq"))
bam_fname <- "/Users/helen/Analyses/Zebrafish/Gagnon/SRR1264585_mg_s.bam"


#all_targets <- read.table("/Users/helen/Desktop/all_target_mappings.bed", sep = "\t")
#target_chr <- gsub(":.*","",gsub("[0-9]+_chr","", all_targets$V4))
#targets <- all_targets[all_targets$V1 == target_chr,]
#target_number <- gsub("_.*", "", targets$V4)
#ambiguous <- which(target_number %in% names(which(table(target_number) > 3)))
targets_f <- "/Users/helen/Analyses/Zebrafish/Gagnon/manual_targets.bed"
#write.table(targets[-ambiguous,], file = targets_f, sep = "\t", row.names = FALSE, 
#            quote = FALSE, col.names = FALSE)

# PCR primers for target 60 are really repetitive

targets <- read.table(targets_f, sep = "\t", stringsAsFactors = FALSE)
# On target = between pcr primers and covers target region
all_ranges <- GRanges(targets$V1, IRanges(targets$V2, targets$V3))

temp <- lapply(split(all_ranges, ceiling(1:nrow(targets)/3)), range)
pcr_ranges <- do.call(c,unlist(temp, use.names = FALSE))
strands <- strand(targets$V6[grep("forward", targets$V4)])
strand(pcr_ranges) <- strands






all_ranges[- grep("target", targets$V4)]

gsub("_.*", "", targets$V4)




# Strategy 1: read everything:
bam <- readGAlignments(bam_fname, param = param, use.names = TRUE) # 2.900

# Get the names, only find chimeras for on target alignments
names_only <- ScanBamParam(what = c("qname"))
system.time(x <- scanBam(bam_fname, param = names_only)) # 1.750

system.time(x <- which(table(names(bam)) > 1)) # 36 sec

system.time(x <- sort(names(bam))) # ~ 39 sec

system.time(x <- order(names(bam))) # ~ 39 sec

system.time(x <- split(1:length(bam), names(bam))) # ~ 35 sec

system.time(x <- which(names(bam) %in% names(bam)[duplicated(names(bam))])) # ~ 0.2 sec


# It's quicker to sort/order a factor than a character vector
system.time(x <- order(names(bam)[chimera_idxs])) # 3.7 sec
system.time(x <- order(as.factor(names(bam)[chimera_idxs]))) # 1.8 sec

system.time(x <- sort(names(bam)[chimera_idxs])) # 3.6 sec
system.time(x <- sort(as.factor(names(bam)[chimera_idxs]))) # 2.4 sec

# Reduce to smaller object then subset better than subset and reduce
system.time(nms <- names(bam[chimera_idxs])) # ~ 0.2 sec
system.time(nms2 <- names(bam)[chimera_idxs]) # ~ 0.005 sec

> all.equal(nms, nms2)
[1] TRUE



system.time(sqs <- seqnames(bam[chimera_idxs])) #0.569
system.time(sqs <- seqnames(bam)[chimera_idxs]) #0.149


system.time(strds <- strand(bam)[chimera_idxs]) # 0.141
system.time(strds <- strand(bam[chimera_idxs])) # 0.598 


system.time(del_lns1 <- start(bam[chimera_idxs[-1]]) - end(bam[chimera_idxs][-length(chimera_idxs)])) # 1.552 
system.time(del_lns2 <- start(bam)[chimera_idxs[-1]] - end(bam)[chimera_idxs][-length(chimera_idxs)]) # 0.067 
system.time(del_lns3 <- start(bam)[chimera_idxs[-1]] - end(bam)[chimera_idxs[-length(chimera_idxs)]]) # 0.062

system.time(cigars <- cigar(bam[chimera_idxs]))
system.time(cigars <- cigar(bam)[chimera_idxs])

system.time(chimeras1 <- bam[chimera_idxs][one_chr & same_strd & has_genome_gap & has_read_gap]) # 0.660
system.time(chimeras2 <- bam[chimera_idxs[one_chr & same_strd & has_genome_gap & has_read_gap]]) # 0.111

system.time(hits_pcr <- findOverlaps(bam, pcr_w_tol, type = "within", ignore.strand = TRUE)) # 2.688
system.time(hits_pcr2 <- subsetByOverlaps(bam, pcr_w_tol, type = "within", ignore.strand = TRUE)) # 3.755
system.time(hits_pcr3 <- findOverlaps(bam, pcr_targets) #2.123
system.time(hits_pcr4 <- findOverlaps(as(bam, "GRanges"), pcr_targets, ignore.strand = TRUE)) # 1.379 
system.time(hits_pcr5 <- findOverlaps(as(bam, "GRanges"), pcr_targets, ignore.strand = TRUE, type = "equal")) # 1.17
system.time(hits_pcr6 <- findOverlaps(as(bam, "GRanges"), pcr_targets, ignore.strand = TRUE, type = "equal", maxgap = 10)) # 1.922


 #nms <- rle(names(bam[rr]))
    #same_name <- rep(nms$lengths > 1, nms$lengths)
    #chimeras <- rr[same_name]
    
       # Merge long gap chimeras

    # Remove the clipped regions
    #unclipped <- width(cigarRangesAlongPairwiseSpace(cigars)) > 0
    


