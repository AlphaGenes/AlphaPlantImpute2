library(data.table)
getFile = function(fileName) {
    mat = as.matrix(fread(fileName))
    mat = mat[order(mat[,1]),]
    return(mat)
}

getSubset = function(mat, ids){
    mat = mat[mat[,1] %in% ids,]
    mat = mat[order(mat[,1]),]
    return(mat)
}

true = getFile("genotypes/trueGenotypes.txt")
masked = getFile("genotypes/maskedGenotypes.txt")

pedigree = getFile("genotypes/pedigree.txt")
pedigree = cbind(pedigree,-1)
getGeneration = function(index) {
    if(pedigree[index,2] == 0) return(0)
    sire = pedigree[index, 2]
    return(pedigree[pedigree[,1] == sire, 4] + 1)
}

for(i in 1:nrow(pedigree)) {
    pedigree[i,4] = getGeneration(i)
}
generations = pedigree[,4]

ids = pedigree[,1]


stratifyByGeneration = function(true, masked, mat) {
    vals = lapply(unique(generations), function(gen) {
        subTrue = true[pedigree[,4] == gen,-1]
        subMat = mat[pedigree[,4] == gen,-1]

        accuracy = sapply(1:nrow(subMat), function(index){
            cor(subMat[index,], subTrue[index,])
        })


        return(c(gen, mean(accuracy)))
    })

    vals = do.call("rbind", vals)
    colnames(vals) = c("gen", "acc")
    return(vals)
}


getAccuracy= function(fileName){
    mat = getFile(fileName)
    mat = getSubset(mat, ids)
    print(fileName)
    print(stratifyByGeneration(true, masked, mat))
}

print("Genotypes")
getAccuracy("out.dosages")
