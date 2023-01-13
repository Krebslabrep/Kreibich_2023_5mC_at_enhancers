#Important arguments for SM call and analysis

#Arguments CA analysis
cO_bin          = 10        #coverage cutoff for bins 
cO_CMH          = 30        #coverage cutoff for SM analysis
WDW_SIZE        = 101       #CG window size
ADD2            = "min2"    #number of replicates needed for CMH test
AVME.i.1        = 0.1       #average methylation - lower cutoff
AVME.i.2        = 0.6       #average methylation - upper cutoff
COR.i.antago    = 0.5       #CMH test: common odds ratio cutoff for antagonist sites
COR.i.ago       = 2         #CMH test: common odds ratio cutoff for agonist sites
PVAL.i          = 0.05      #CMH test: pvalue cutoff for antagonist sites
COR.ntrl.i.1    = 0.9000    #CMH test: lower common odds ratio cutoff for neutral+ sites
COR.ntrl.i.2    = 1.1111    #CMH test: upper common odds ratio cutoff for neutral+ sites
PVAL.ntrl.i     = 1         #CMH test: pvalue cutoff for neutral+ sites


level_state0  <- c("neutral", "antagonist", "agonist", "low_5mC", "high_5mC")
level_state   <- c("neutral", "antagonist", "ICR", "agonist", "low_5mC", "high_5mC")
level_state2  <- c("neutral", "neutral+", "antagonist", "ICR", "agonist", "low_5mC", "high_5mC")
