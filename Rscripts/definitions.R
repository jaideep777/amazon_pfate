exprlabel = function(a,b,c=NULL){
  a = substitute(a)
  b = substitute(b)
  c = substitute(c)
  
  if (is.null(c)){
    expr(
      atop(
        textstyle(!!a), 
        textstyle(!!b)
      )
    ) 
  } else {
    expr(
      atop(
        atop(
          textstyle(!!a), 
          textstyle(!!b)
        ),
        textstyle(!!c)
      )
    )
  }
}


labels = c(
  GPP   = exprlabel("Gross primary", "productivity", "(kgC m"^"-2"~"yr"^"-1"*")"),
  NPP   = exprlabel("Net primary", "productivity", "(kgC m"^"-2"~"yr"^"-1"*")"),
  GS    = exprlabel("Stomatal", "conductance", "(mol m"^"-2"~"s"^"-1"*")"),
  VCMAX = exprlabel("V"["cmax"], "("*mu*"mol m"^"-2"~"s"^"-1"*")"),
  BA    = exprlabel("Basal area", "(m"^2~"ha"^"-1"*")"),
  AGB   = exprlabel("Aboveground", "biomass", "(kgC ha"^"-1"*")"),
  CFR   = exprlabel("Fine root", "biomass", "(kgC ha"^"-1"*")"),
  LMA   = exprlabel("Leaf mass", "per area", "(g m"^"-2"*")"),
  HMAT  = exprlabel("Max. height", "(m)"),
  WD    = exprlabel("Wood density", "(kg m"^"-3"*")"),
  P50X  = exprlabel("Xylem P"[50],"(MPa)"),
  LAI   = exprlabel("Leaf area", "index", "(-)"),
  RAU   = exprlabel("Autotrophic", "respiration", "(kgC m"^"-2"~"yr"^"-1"*")"),
  Z     = exprlabel("Canopy layer", "heights (m)"),
  TRANS = exprlabel("Transpiration", "(mm day-1)"),
  MORT  = exprlabel("Biomass", "mortality rate", "(kgC m"^"-2"~"yr"^"-1"*")")
)

amz_theme = function(){
  theme(strip.placement = "outside",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
}

