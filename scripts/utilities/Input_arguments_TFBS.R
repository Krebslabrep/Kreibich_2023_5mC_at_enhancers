#Important arguments for SM call and analysis

#Arguments TFBS analysis
cO_TFBS       = 10    #coverage cutoff per TFBS
cO_fraction   = 5     #coverage cutoff per fraction
cO_BOUND      = 0.05  #cutoff bound fraction 5%
WDW_SIZE      = 30    #window size of TFBS for nearby CpG
ADD           = "_min2_closed2_NWCGW_2"

COR.i.antago  = 0.5       #CMH test: common odds ratio cutoff for antagonist sites
COR.i.ago     = 2         #CMH test: common odds ratio cutoff for agonist sites
PVAL.i        = 0.05  #pvalue cutoff for antagonist sites
MEDELTA.i     = 0.2   #methylation delta cutoff for antagonist sites (NOT USED ANYMORE WITH CMH test)
AVME.i.1      = 0.1   #lower average methylation cutoff (CMH cannot handle high/low methylation)
AVME.i.2      = 0.8   #upper average methylation cutoff (CMH cannot handle high/low methylation)
COR.i.neut.1  = 0.8   #lower common odds ratio cutoff for neutral sites
COR.i.neut.2  = 1.2   #higher common odds ratio cutoff for neutral sites


LEVELS_state <- c("fully methylated", "antagonist", "agonist", "neutral",  "unmethylated")
LEVELS_state2 <- c("fully methylated", "antagonist", "agonist", "neutral+",  "neutral", "unmethylated")
LEVELS_state_crude = c("fully methylated", paste("COR >", COR.i.antago),  paste("COR <", COR.i.antago), "unmethylated")
