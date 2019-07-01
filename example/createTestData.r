library(AlphaSimR)
options(scipen=999)
args = commandArgs(trailingOnly=TRUE)


# print(args)
founderPop = runMacs(nInd = 100, nChr = 1, segSites = 1000)

# ---- Set Simulation Parameters ----

SP = SimParam$new(founderPop)
SP$addTraitA(nQtlPerChr=1000, mean=0, var=1)
SP$addSnpChip(nSnpPerChr=1000)
SP$setVarE(h2=.4)


SP$setTrackPed(TRUE)
SP$setTrackRec(TRUE)



nGenerations = 2
generationOutput = vector("list",nGenerations)


generationOutput[[1]] = newPop(founderPop)
for(gen in 2:nGenerations){ 
    sires = selectInd(generationOutput[[gen-1]], nInd = 25, gender = "M")
    dams = selectInd(generationOutput[[gen-1]], nInd = 50, gender = "F")
    generationOutput[[gen]] = randCross2(females = dams, males = sires, nCrosses = 100, nProgeny = 1)
}

output = do.call("c",generationOutput)

write.table2 = function(x, file, ...) {
    write.table(x, file, col.names=F, row.names=F,quote=F,...)
}
pedigree = data.frame(output@id, output@father, output@mother, stringsAsFactors = FALSE)
genotypes = AlphaSimR::pullSnpGeno(output, 1)
haplotypes = AlphaSimR::pullSnpHaplo(output, 1)

write.table2(pedigree, "genotypes/pedigree.txt")
write.table2(cbind(pedigree[,1], genotypes), "genotypes/trueGenotypes.txt")
write.table2(cbind(rep(pedigree[,1], each = 2), haplotypes), "genotypes/trueHaplotypes.txt")


maskedGenotypes = genotypes
maskedHaplotypes = haplotypes
# Let's say 20% of them are missing at 50% of the loci

ldIds = 101:200
ldSnps = seq(1, 1000, length.out = 200)

maskedGenotypes[101:200, -ldSnps] = 9
maskedHaplotypes[201:400, -ldSnps] = 9

write.table2(cbind(pedigree[,1], maskedGenotypes), "genotypes/maskedGenotypes.txt")
write.table2(cbind(rep(pedigree[,1], each = 2), maskedHaplotypes), "genotypes/maskedHaplotypes")
