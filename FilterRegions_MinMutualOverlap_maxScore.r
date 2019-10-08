## Fitlering of redundant overlapping regions with minimal mutual overlap and selectino of highest score
##
## This code was tested on r3.5.0
## 
## the function mergeOverlapingDomains takes as input a GRange object of regions <regions> (e.g. HiC regions or TADs) with score(regions) beeing a score for fitlering 
## (e.g. average contact enrichment, or other). 
## <MinMutOverlap> is the minimum mutual overlap or 2 features to be considered to be redundant, as a fraction of domains' length (value between 0 and 1). For each comparison, for each region, the fraction of overlap over region length as to be higher than <MinMutOverlap> for the two domaisn to be considered redundant.
## Two domains overlapping less than <MinMutOverlap> of each other are not redundant. Similarly, A region completely encompassed in another bigger region but smaller than <MinMutOverlap> of the bigger region is not considered redundant.

library(GenomicAlignments)
library(rtracklayer)
library(ggplot2)
library(dplyr)
library(tidyr)

mergeOverlapingDomains <- function(domains, MinMutOverlap, MinWidth){
    domains <- domains[width(domains)>MinWidth]
    names(domains) <- values(domains)[["name"]]
    domainsMerge <- domains
    nbSelectDomain <- 1
    nbDomain <- length(domainsMerge)

    n=0
    while (nbDomain!= nbSelectDomain){
        nbDomain <- length(domainsMerge)
        domainsMergeOlap <- findOverlaps(domainsMerge, domainsMerge)

        qHits=queryHits(domainsMergeOlap)[queryHits(domainsMergeOlap)!=subjectHits(domainsMergeOlap)]
        sHits=subjectHits(domainsMergeOlap)[queryHits(domainsMergeOlap)!=subjectHits(domainsMergeOlap)]

        domainsMergepIntersect <- pintersect(domainsMerge[qHits], domainsMerge[sHits])

        qHits_Olap_MinOlap <- qHits[width(domainsMergepIntersect)/width(domainsMerge[qHits])>=MinMutOverlap & width(domainsMergepIntersect)/width(domainsMerge[sHits])>=MinMutOverlap]
        sHits_Olap_MinOlap <- sHits[width(domainsMergepIntersect)/width(domainsMerge[qHits])>=MinMutOverlap & width(domainsMergepIntersect)/width(domainsMerge[sHits])>=MinMutOverlap]


        qHits_Olap_MinOlap_highestScore <- qHits_Olap_MinOlap[score(domainsMerge[qHits_Olap_MinOlap])>= score(domainsMerge[sHits_Olap_MinOlap])]
        sHits_Olap_MinOlap_highestScore <- sHits_Olap_MinOlap[score(domainsMerge[sHits_Olap_MinOlap])> score(domainsMerge[qHits_Olap_MinOlap])]

        Hits_Olap_MinOlap_highestScore <- unique(c(qHits_Olap_MinOlap_highestScore, sHits_Olap_MinOlap_highestScore))

        Hits_Olap_MinOlap <- unique(c(qHits_Olap_MinOlap, sHits_Olap_MinOlap))
        noOlap_MinOlap <- names(domainsMerge)[!(names(domainsMerge) %in% names(domainsMerge[Hits_Olap_MinOlap]))]

        domainsMerge_select <-  domainsMerge[c(names(domainsMerge)[Hits_Olap_MinOlap_highestScore], noOlap_MinOlap)]
        print(paste0(n, "    |    domainsMerge: ", length(domainsMerge), "    |    domainsMerge_select: ", length(domainsMerge_select)))
        nbSelectDomain <- length(domainsMerge_select)
        domainsMerge <- domainsMerge_select
        n=n+1
    }
    return(sort(domainsMerge))
}
