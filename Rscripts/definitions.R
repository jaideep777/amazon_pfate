# exprlabel = function(a,b,c=NULL){
#   a = substitute(a)
#   b = substitute(b)
#   c = substitute(c)
#   
#   if (is.null(c)){
#     expr(
#       atop(
#         textstyle(!!a), 
#         textstyle(!!b)
#       )
#     ) 
#   } else {
#     expr(
#       atop(
#         atop(
#           textstyle(!!a), 
#           textstyle(!!b)
#         ),
#         textstyle(!!c)
#       )
#     )
#   }
# }
# 
# 
# labels = c(
#   GPP   = exprlabel("Gross primary", "productivity", "(kgC m"^"-2"~"yr"^"-1"*")"),
#   NPP   = exprlabel("Net primary", "productivity", "(kgC m"^"-2"~"yr"^"-1"*")"),
#   GS    = exprlabel("Stomatal", "conductance", "(mol m"^"-2"~"s"^"-1"*")"),
#   VCMAX = exprlabel("V"["cmax"], "("*mu*"mol m"^"-2"~"s"^"-1"*")"),
#   BA    = exprlabel("Basal area", "(m"^2~"ha"^"-1"*")"),
#   AGB   = exprlabel("Aboveground", "biomass", "(kgC ha"^"-1"*")"),
#   CFR   = exprlabel("Fine root", "biomass", "(kgC ha"^"-1"*")"),
#   LMA   = exprlabel("Leaf mass", "per area", "(g m"^"-2"*")"),
#   HMAT  = exprlabel("Max. height", "(m)"),
#   WD    = exprlabel("Wood density", "(kg m"^"-3"*")"),
#   P50X  = exprlabel("Xylem P"[50],"(MPa)"),
#   LAI   = exprlabel("Leaf area", "index", "(-)"),
#   RAU   = exprlabel("Autotrophic", "respiration", "(kgC m"^"-2"~"yr"^"-1"*")"),
#   Z     = exprlabel("Canopy layer", "heights (m)"),
#   TRANS = exprlabel("Transpiration", "(mm day-1)"),
#   MORT  = exprlabel("Biomass", "mortality rate", "(kgC m"^"-2"~"yr"^"-1"*")")
# )

labels3 = c(
  GPP   = "Gross<br>productivity<br>(kgC m<sup>&minus;2</sup> yr<sup>&minus;1</sup>)",
  NPP   = "Net<br>productivity<br>(kgC m<sup>&minus;2</sup> yr<sup>&minus;1</sup>)",
  GS    = "Stomatal<br>conductance,<br>*g*<sub>c</sub> (mol m<sup>&minus;2</sup> s<sup>&minus;1</sup>)",
  VCMAX = "Photosynthetic<br>capacity,<br>*V*<sub>cmax,25</sub> (&mu;mol m<sup>&minus;2</sup> s<sup>&minus;1</sup>)",
  BA    = "<br>Basal area<br>(m<sup>2</sup> ha<sup>&minus;1</sup>)",
  AGB   = "Aboveground<br>biomass<br>(kgC ha<sup>&minus;1</sup>)",
  CFR   = "Fine root<br>biomass<br>(kgC ha<sup>&minus;1</sup>)",
  LMA   = "Leaf mass<br>per area<br>(g m<sup>&minus;2</sup>)",
  HMAT  = "<br>Max. height<br>(m)",
  WD    = "<br>Wood density<br>(kg m<sup>&minus;3</sup>)",
  P50X  = "Xylem<br>vulnerability,<br>*&psi;*<sub>50</sub> (MPa)",
  LAI   = "Leaf area<br>index<br>(-)",
  RAU   = "Autotrophic<br>respiration<br>(kgC m<sup>&minus;2</sup> yr<sup>&minus;1</sup>)",
  Z     = "Canopy layer<br>heights,<br>z<sup>*</sup> (m)",
  TRANS = "<br>Transpiration<br>(mm day<sup>&minus;1</sup>)",
  MORT  = "Biomass<br>mortality rate<br>(kgC m<sup>&minus;2</sup> yr<sup>&minus;1</sup>)",
  z2     = "Canopy layer 2<br>height<br>(m)",
  z3     = "Canopy layer 3<br>height<br>(m)",
  SEEDS  = "<br>Seed rain<br>(yr<sup>&minus;1</sup>)"
)

multi_breaks = str_count(labels3, pattern = "<br>")==2
labels2 = labels3
labels2[multi_breaks] = stringr::str_replace(labels3[multi_breaks], "<br>", " ")
names(labels2) = names(labels3)

unitstart = labels3 %>% str_locate_all(pattern = "<br>") %>% lapply(function(x){x[2,"start"]}) %>% unlist() 
labels_nounit = substr(labels3, 1, unitstart-1) %>% str_replace_all(",","")
names(labels_nounit) = names(labels3)

amz_theme = function(){
  theme_bw()+
    theme(strip.placement = "outside",
          strip.text = ggtext::element_markdown(lineheight=1.2),
          axis.title = ggtext::element_markdown(lineheight=1.2),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank())
}

col_ele = "#fa9500"
col_obs = "#0065fa"
col_amb = "grey40"
col_amb_dark = "grey10"
