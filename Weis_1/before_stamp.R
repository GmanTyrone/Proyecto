pathway <- read.table("L2_3_path_abun_unstrat_descrip.tsv", sep = "\t", header = TRUE)
pathway2 <- pathway[pathway[,2] != "not_found",] 
pathway2 <- pathway2[,-1]
write.table(pathway2, file = "pathway_stamp.tsv", row.names=FALSE, sep="\t", quote = FALSE)
