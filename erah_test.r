# erah test
library(erah)
createdt(path = 'erah_test/')
ex <- newExp(instrumental = 'erah_test/erah_test_inst.csv', phenotype = 'erah_test/erah_test_pheno.csv')
metaData(ex)
phenoData(ex)

ex.dec.par <- setDecPar(min.peak.width = 5, avoid.processing.mz = c(73:75, 147:149))
ex <- deconvolveComp(ex, ex.dec.par)
ex1 <- ex

ex.al.par <-setAlPar(min.spectra.cor = 1, max.time.dist = 3, mz.range = 50:600)
ex <- alignComp(ex, alParameters = ex.al.par)

ex <- identifyComp(ex, mz.range = 50:600)

id.list <- idList(ex)
head(id.list)
