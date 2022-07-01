library(data.table)
library(dplyr)
library(TFBSTools)
library(ggseqlogo)
library(motifmatchr)

id = "M001_0.6"
motifs <- gsub(".txt", "", list.files("pwms_all_motifs/", pattern = ".txt"))

lapply(motifs, function(id){
  print(id)
  x <- fread(paste0("pwms_all_motifs/",id,".txt"))[,-1] %>%
    data.frame() %>% data.matrix
  colnames(x) <- c("A", "C", "G", "T")
  pwm <- PFMatrix(ID = id, name = id,profileMatrix = t(x)*10000) %>%
    toPWM
  pwm
}) -> list_of_pwms

pfmList <- do.call(PWMatrixList, list_of_pwms)
saveRDS(pfmList, file = "cisbp_RNA.motifs.rds")

# Verify it works with motif matchr
matchMotifs(pfmList, "AAATTTCCAAAGGGAA")
#yay