#mediation analysis: microbial impacts on host phenotypes through metabolites

library(lme4)
library(mediation)
model.m=lmer(meta~microbe+age+sex+(1|time),data)
model.y=lmer(pheno~microbe+meta+age+sex+(1|time),data)
summary=summary(mediate(model.m ,model.y,treat = "microbe", mediator = "meta",boot = F,sims = 1000))