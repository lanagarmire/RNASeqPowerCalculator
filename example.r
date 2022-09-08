#my_directory <- ""
#setwd(my_directory)
source("rs_simulations.r")

read.table("bottomly_count_table.txt", sep = "\t", header=T) -> rawdata
rownames(rawdata) <- rawdata$gene
rawdata <- as.matrix(rawdata[,-1])
head(rawdata)

condition <- c(rep("A", 10), rep("B", 11))
print(condition)

estimate_params(rawdata=rawdata, condition=condition, designtype="one factor") -> params
save(params, file="bottomly_params.Rdata")

RS_simulation(sims=5, params=params, budget=3000, designtype = "one factor", nmax = 20, nmin = 2, program="DESeq2") -> results

plot(rownames(results),rowMeans(results, na.rm=T), main="DESeq simulations on Bottomly Dataset", xlab = "number of replicates", ylab = "Power")
