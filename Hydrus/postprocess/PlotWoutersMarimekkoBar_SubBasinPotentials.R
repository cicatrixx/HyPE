#Nat=natural flow                 = Full_GWh_year
#Cons= 100-(nat - cons)/cons*100  = Cons
#Season=Dry/Wet                   = pottype
#Basin= Basin name                = basin
#scenyear=scenario name for panels= pottype

## try w matfile
# library(R.matlab)
# pathname <- 'G:/SurfDrive/HPmodel/output/Figs_trial/MainScenarios4WouterBar.mat'
# data <- readMat(pathname)
# print(data)

##----
#ColorOrder: Ora#ColorOrder: Orange, Sky blue, BluishGreen, Yellow, Blue, Vermillion, Black, Reddish purple
cmap8_Wong<-c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#000000", "#CC79A7")
myroot<- "G:/SurfDrive/HPmodel"   # "D:/HPmodel"  #
setwd(myroot)
library(readxl) # For reading
library(dplyr)  #for data wrangling
library(ggplot2) #for main plot
library(gridExtra) #for grid tiling
library(grid)
idata<-read_excel(file.path(myroot,"output","Figs_trial","SubBasinData4WouterBar.xlsx"))

idata$basin<-factor(idata$basin, unique(idata$basin))
idata$planttype<-factor(idata$planttype,unique(idata$planttype))
#idata$pottype<-factor(idata$pottype,unique(idata$pottype)[c(1,5,2,3,4)])
idata$pottype<-factor(idata$pottype,unique(idata$pottype))

idata$Cons<-idata$Remain_GWh_year/idata$Full_GWh_year*100; 

subBasinTots<-idata #data.frame()
subBasinTots$Nat<-idata$Full_GWh_year/1000 # Convert to TWh


subBasinTots_RPonly<-subBasinTots %>% filter(pottype!="Theoretical" & pottype!="Visualized") %>% group_by(pottype) %>% filter(planttype=="River Power Plant") %>%  summarise(Nat=sum(Nat,na.rm = TRUE))

#=====
myMarimekkoTheme<-function(){
  theme(panel.grid = element_blank(),  #wouters gray color =#f5f5f2
        plot.background = element_rect(fill = "NA" , colour = "NA"),
        panel.background = element_rect(fill = "NA",color = "NA"),
        strip.background = element_rect(fill = "NA"),
        legend.background = element_rect(fill = "NA"),
        panel.grid.major.y = element_line(color = alpha("gray40", 0.7), size = 0.4),
        panel.grid.minor.y = element_line(color = alpha("gray60", 0.7), size = 0.4),
        #panel.grid.major.x = element_line(color = alpha("gray40", 0.7), size = 0.4),
        #text=element_text(family="Helvetica"),
        axis.ticks = element_line(color="black"), 
        axis.ticks.length=unit(1,"mm"), 
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(color="black", size = 0.5), 
        axis.text.y = element_text(color="black", size=10),
        axis.text.x = element_text(margin=margin(5,0,0,0), color="gray30", size=10), 
        axis.title = element_text(color="black", size=10, face="bold"),
        strip.text = element_text(color="black", size=9, face="bold", hjust = 0.5), 
        panel.spacing = unit(0, "cm"),
        legend.title = element_text(color="black", size=10, face="plain"),
        legend.text = element_text(color="black", size=8, face="plain"),
        legend.position = "bottom", 
        legend.margin=margin(0,2,0,0),
        legend.box.margin=margin(-5,0,0,0), 
        legend.key.width = unit(1.2, "cm"), 
        legend.key.height = unit(0.4, "cm"), 
        legend.key = element_rect(fill = NA),
        legend.spacing.x = unit(0.02, "cm"), 
        legend.justification = "right") 
}

vispotticks<-seq(0,300,50);
vispotlim<-c(-5,320);
theorypotticks<-seq(0,1500,200);
theorypotlim<-c(-5,1600);
prctticks<-seq(0,100,50); #c(0,25,50,75,100)
prctlim<-c(-5,110);

my_alpha_fill<-0.9
my_ylabel<-"Full potential (TWh/yr)"; #xy is flipped here!
my_xlabel<-"Remaining potential as % of full potential"; 

#====
# Plot Tech, Fin, Sust ====
p_TFS<-ggplot(filter(subBasinTots,pottype!="Theoretical" & pottype!="Visualized"), aes(x = Cons/2, y = Nat, width = Cons, group = planttype))+
  #create outer volume boxes stacked horizontally
  geom_col(position=position_stack(reverse = T), stat="identity", aes(x = 100/2, width = 100, fill = basin, col = basin),fill = NA, alpha = 1, lwd = 0.5) +
  #create filled boxes stacked horizontally
  geom_col(position=position_stack(reverse = T), stat="identity", alpha = my_alpha_fill, col = NA, aes(fill = basin, col = basin)) +
  # Add lines and labels for two plant types
  geom_hline(data = subBasinTots_RPonly[subBasinTots_RPonly$pottype == "Technical",], aes(yintercept = Nat), x = c(-3, 130), col ="grey40", linetype = "33", lwd = .04)+
  geom_hline(data = subBasinTots_RPonly[subBasinTots_RPonly$pottype != "Technical",], aes(yintercept = Nat), x = c(-3, 105), col ="grey40", linetype = "33", lwd = .04)+ #sprintf('Diversion canal plant \u2192')
  geom_text(data = subBasinTots_RPonly[subBasinTots_RPonly$pottype == "Technical",], aes(y = Nat - 35), x = 108, col = "grey40", label =  sprintf(' \u2190 River power plant') , size=3, fontface = "bold", inherit.aes = F)+
  geom_text(data = subBasinTots_RPonly[subBasinTots_RPonly$pottype == "Technical",], aes(y = Nat + 36), x = 108, col = "grey40", label = sprintf('Diversion canal plant \u2192'), size=3, fontface = "bold", inherit.aes = F)+ 
  # make plot pretty, add labels and legends
  facet_grid(rows = vars(pottype), )+
  #facet_wrap(rows = vars(pottype), scales="free_y")+
  labs(x = my_xlabel, y = my_ylabel)+
  scale_y_continuous(breaks = vispotticks, limits = vispotlim)+
  # Same format for all figures
  scale_x_continuous(breaks = prctticks, limits = prctlim)+
  coord_flip(expand = F)+
  scale_fill_manual(values=cmap8_Wong, name = "Sub-basins  ", guide = 'none')+
  scale_color_manual(values=cmap8_Wong, name = "Sub-basins ", guide = 'none')+
  #scale_fill_viridis_d(option = "inferno", begin = 0.1111111, end = 0.8888889, name = "Basin  ")+
  #scale_color_viridis_d(option = "inferno", begin = 0.1111111, end = 0.8888889, name = "Total flow (km3)", guide = 'none')+
  myMarimekkoTheme()

# Plot ONLY THEORY====
p_theory<-ggplot(filter(subBasinTots,pottype=="Theoretical"), aes(x = Cons/2, y = Nat, width = Cons, group = planttype))+
  #create outer volume boxes stacked horizontally
  geom_col(position=position_stack(reverse = T), stat="identity", aes(x = 100/2, width = 100, fill = basin, col = basin),fill = NA, alpha = 1, lwd = 0.5) +
  #create filled boxes stacked horizontally
  geom_col(position=position_stack(reverse = T), stat="identity", alpha = my_alpha_fill, col = NA, aes(fill = basin, col = basin)) +
  ##make plot pretty, add labels and legends
  facet_grid(rows = vars(pottype), )+
  #facet_wrap(rows = vars(pottype), scales="free_y")+
  labs(x ="", y = my_ylabel)+
  scale_y_continuous(breaks = theorypotticks, limits = c(-20,1600))+
  # fill w manual colors
  scale_x_continuous(breaks =prctticks, limits = prctlim)+
  coord_flip(expand = F)+
  scale_fill_manual(values=cmap8_Wong, name = "Sub-basins  ", guide = 'none')+
  scale_color_manual(values=cmap8_Wong, name = "Sub-basins ", guide = 'none')+
  myMarimekkoTheme()

# Plot ONLY Visualized_hatched====
p_viz<-ggplot(filter(subBasinTots,pottype=="Visualized"), aes(x = Cons/2, y = Nat, width = Cons, group = planttype))+
  #create outer volume boxes stacked horizontally
  geom_col(position=position_stack(reverse = T), stat="identity", aes(x = 100/2, width = 100, fill = basin, col = basin),fill = NA, alpha = 1, lwd = 0.5) +
  #create filled boxes stacked horizontally
  geom_col(position=position_stack(reverse = T), stat="identity", alpha = my_alpha_fill, col = NA, aes(fill = basin, col = basin))+
  ##make plot pretty, add labels and legends
  facet_grid(rows = vars(pottype))+
  #facet_wrap(rows = vars(pottype), scales="free_y")+
  labs(x ="", y = my_ylabel)+
  scale_y_continuous(breaks = vispotticks, limits = vispotlim)+
  # Same format for all figures
  scale_x_continuous(breaks = prctticks, limits = prctlim)+
  coord_flip(expand = F)+
  # fill w manual colors
  scale_fill_manual(values=cmap8_Wong, name = "Sub-basins ")+ #
  scale_color_manual(values=cmap8_Wong, name = "Sub-basins ", guide = 'none')+
  myMarimekkoTheme()+
  guides(fill = guide_legend(title.position = "left", order = 3, nrow = 1,                             
                             title.vjust = 0.9, label.position = "bottom"))

# Plot ONLY Visualized_hatched====
p_viz_hatched<-ggplot(filter(subBasinTots,pottype=="Visualized"), aes(x = Cons/2, y = Nat, width = Cons, group = planttype))+
  #create outer volume boxes stacked horizontally
  geom_col(position=position_stack(reverse = T), stat="identity", aes(x = 100/2, width = 100, fill = basin, col = basin),fill = NA, alpha = 1, lwd = 0.5) +
  #create filled boxes stacked horizontally
  ggpattern::geom_col_pattern(position=position_stack(reverse = T), stat="identity", alpha = my_alpha_fill, col = NA, aes(fill = basin, col = basin,pattern_colour=basin),
                              pattern = 'stripe',pattern_fill='white',pattern_alpha=my_alpha_fill,pattern_angle=45,  pattern_density = 0.6) +
  ggpattern::scale_pattern_color_manual(values = cmap8_Wong, guide = 'none')+ #https://coolbutuseless.github.io/package/ggpattern/index.html
  ##make plot pretty, add labels and legends
  facet_grid(rows = vars(pottype))+
  #facet_wrap(rows = vars(pottype), scales="free_y")+
  labs(x ="", y = my_ylabel)+
  scale_y_continuous(breaks = vispotticks, limits = vispotlim)+
  # Same format for all figures
  scale_x_continuous(breaks = prctticks, limits = prctlim)+
  coord_flip(expand = F)+
  # fill w manual colors
  scale_fill_manual(values=cmap8_Wong, name = "Sub-basins ")+ #
  scale_color_manual(values=cmap8_Wong, name = "Sub-basins ", guide = 'none')+
  myMarimekkoTheme()+
  guides(fill = guide_legend(title.position = "left", order = 3, nrow = 1,                             
                             title.vjust = 0.9, label.position = "bottom"))


#####
# Plot individually =====
p_TFS
p_theory
p_viz
#save fig
#ggsave(filename = file.path(getwd(),"output","Figs_trial","Sub-basinMarimekkoBar_Onlytheory.pdf"),plot=p_theory,width=7, height=2.5,units=c("in"), device = cairo_pdf,dpi=300)

#interactive
#library(plotly)
#ggplotly(p)

#save fig
# ggsave(filename = file.path(getwd(),"Sub-basinMarimekkoBar.pdf"),plot=p,width =8, height=5)

# Plot tiled ====
#Tile my plots ### 0---did not work
# g1 <- ggplotGrob(p_theory)
# g2 <- ggplotGrob(p_TFS)
#g3 <- ggplotGrob(p_viz)
#g <- cbind(g2, g3, size = "first")

# Tile using grid.arrange
#grid.arrange(g1, g2, g3, ncol=1) #equal size subplots
#gall<-grid.arrange(p_theory, p_TFS, p_viz, ncol=1,heights=c(1,3,1.2)) #https://cran.r-project.org/web/packages/gridExtra/vignettes/arrangeGrob.html
gall<-grid.arrange(p_theory, p_TFS, p_viz, ncol=1,heights=c(6,4*3+2+1,6+1.5))
#ggsave(plot=gall,filename = file.path(getwd(),"output","Figs_trial","cl","F4_Sub-basinMarimekkoBar_All5.pdf"),width=10, height=7,units=c("in"), device = cairo_pdf)

gall_hatched<-grid.arrange(p_theory, p_TFS, p_viz_hatched, ncol=1,heights=c(6,4*3+2+1,6+1.5))
#ggsave(plot=gall_hatched,filename = file.path(getwd(),"output","Figs_trial","cl","F4_Sub-basinMarimekkoBar_All5_hatched.pdf"),width=10, height=7,units=c("in"), device = cairo_pdf)

# stitch them together as in https://stackoverflow.com/a/53471718
#g00 <- grid.arrange(g3, g3, nrow = 2)
# final touch
# g200 <- cowplot::ggdraw(g00) + 
#   theme(plot.background = element_rect(fill="white", color = "white"))

#interactive
library(plotly)
ggplotly(p_theory)
