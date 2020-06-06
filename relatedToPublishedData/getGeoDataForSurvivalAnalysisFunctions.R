#######################################################################
#######################################################################
##this function performs survival analysis for input genes
##based on provided expression and phenotype data
#######################################################################
survivalAnalysis = function(geneOfInterest, expressionData, phenotypeData,
                            survivalPlot = FALSE, survivalData = FALSE, writeDirectory = '', ...){
  ##
  geneOfInterestExprs = expressionData %>%
    filter(grepl(geneOfInterest, symbol))
  geneOfInterestExprs$medExprs = apply(geneOfInterestExprs[,4:ncol(geneOfInterestExprs)], 1, function(x) median(x, na.rm = TRUE))
  geneOfInterestSort = geneOfInterestExprs %>%
    arrange(desc(medExprs))
  ##
  geneSurvivalInput = geneOfInterestSort[1,1:(ncol(geneOfInterestSort) - 1)] %>%
    pivot_longer(cols = colnames(geneOfInterestSort)[4]:colnames(geneOfInterestSort)[(ncol(geneOfInterestSort) - 1)], 
                 names_to = 'geo_accession', values_to = 'rnaExprs') %>%
    right_join(phenotypeData) %>%
    arrange(desc(rnaExprs)) 
  ##
  geneSurvivalInput$geneLevel = ''
  geneSurvivalInput[1:round(nrow(geneSurvivalInput) * 0.25), 'geneLevel'] = 'high'
  geneSurvivalInput[round(nrow(geneSurvivalInput) * 0.75):nrow(geneSurvivalInput), 'geneLevel'] = 'low'
  geneSurvivalInput$geneLevel = factor(geneSurvivalInput$geneLevel, levels = c('low','high'))
  ##
  
  ##
  survivalFit = survfit(Surv(ovs, status) ~ geneLevel, data = geneSurvivalInput)
  survivalDiff = survdiff(Surv(ovs, status) ~ geneLevel, data = geneSurvivalInput)
  survivalPValue = 1 - pchisq(survivalDiff$chisq, length(survivalDiff$n) - 1)
  coxStats = coxph(Surv(ovs, status) ~ geneLevel, data = geneSurvivalInput)
  coxZScore = coef(coxStats)/sqrt(diag(vcov(coxStats)))
  ##
  if (survivalPlot == TRUE){
    if (nchar(writeDirectory) <= 3){
      message(paste('No write directory specified. Plots will be written to', baseRepository))
    }
    ggsurvplot(survivalFit, pval = survivalPValue, conf.int = FALSE,
               risk.table = FALSE, linetype = "strata", ggtheme = theme_classic(), palette = brewer.pal(8,'RdBu')[c(8,1)])
    ggsave(paste(baseRepository, writeDirectory, '/survival_', 
                 geneOfInterest, '.pdf', sep = ''), width = 4, height = 4, useDingbats = FALSE)
  }
  ##
  if (survivalData == TRUE){
    if (nchar(writeDirectory) <= 3){
      message(paste('No write directory specified. Results will be written to', baseRepository))
    }
    survivalOutput = summary(survivalFit)$table
    head(survivalOutput)
    write.table(survivalOutput, paste(baseRepository, writeDirectory, '/survival_', geneOfInterest, '.csv', sep = ''),
                col.names = TRUE, row.names = TRUE, quote = FALSE, sep = ',')
  }
  ##
  return(coxZScore)
}
#######################################################################
#######################################################################