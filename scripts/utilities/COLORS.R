# TITLE: Plot functions and colors
# AUTHOR: Elisa Kreibich


# Plot design -----------------------------------------------------------------------

theme_EK <- function (base_size = 14, base_family = "") {
  theme_bw() %+replace% 
    theme(
      text = element_text(color = "black"),
      axis.text.x = element_text(color = "black", size = 10, margin = margin(t=5)),
      axis.text.y = element_text(color = "black", size = 10, margin = margin(r=5)),
      axis.ticks.length = unit(2, "mm"),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", size = 1, fill = NA),
      axis.ticks = element_line(color = "black")
    )
}

# Colors -----------------------------------------------------------------------

STATES                  <- c("neutral", "neutral+", "antagonist", "ICR", "high_5mC", "low_5mC")
STATES_Ago              <- c("neutral", "neutral+", "antagonist", "ICR", "agonist", "high_5mC", "low_5mC")

STATES_5mC              <- c("neutral", "neutral+", "antagonist", "fully methylated", "unmethylated")
STATES_Ago_5mC          <- c("neutral", "neutral+", "antagonist", "agonist", "fully methylated", "unmethylated")

STATES_crude_5mC        <- c("COR > 0.5", "COR < 0.5","unmethylated")


COLORS_STATE_v1         <- c("#000000", "#000000","#84334A", "#ED9696", "#F59851", "#468FAF")
names(COLORS_STATE_v1)  <- STATES
COLORS_STATE_v1s        <- COLORS_STATE_v1[1:4]

COLORS_STATE_v2         <- c("grey20", "grey20", COLORS_STATE_v1[-c(1:2)]) 
names(COLORS_STATE_v2)  <- STATES
COLORS_STATE_v2s        <- COLORS_STATE_v2[1:4]

COLORS_STATE_v3         <- c("#000000", "#000000","#84334A", "#ED9696", "#33846D", "#F59851", "#468FAF")
names(COLORS_STATE_v3)  <- STATES_Ago
COLORS_STATE_v3s        <- COLORS_STATE_v3[1:5]

COLORS_STATE_v4         <- c("#000000", "#000000","#84334A", "#33846D", "royalblue4", "grey60")
names(COLORS_STATE_v4)  <- STATES_Ago_5mC

COLORS_STATE_v5         <- c("#000000", "#000000","#84334A", "royalblue4", "grey60")
names(COLORS_STATE_v5)  <- STATES_5mC

COLORS_STATE_v6         <- c("#000000", "#84334A", "grey60")
names(COLORS_STATE_v6)  <- STATES_crude_5mC



TYPES                   <- c("ES WT", "DNMT TKO", "TET TKO")
COLORS_TYPES            <- c("#000000","#5e0017","#c16441")
names(COLORS_TYPES)     <- TYPES

COLORS_TYPES_F1         <- c("#000000","#5e0017","#c16441")
names(COLORS_TYPES_F1)  <- TYPES

COLORS_TYPES_SOMATIC    <- c("#000000","#5e0017","#c16441")
names(COLORS_TYPES_SOMATIC) <- TYPES


COLORS_METH             <- c("#03045e","#023e8a","#0077b6","#0096c7","#00b4d8","#48cae4","#90e0ef","#ade8f4","#caf0f8")


COLOR_SMF               <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(9)[c(4, 3, 2)]
COLOR_SMF               <- c(COLOR_SMF[-3],"black")

