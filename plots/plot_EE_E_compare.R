####
#important note: HEzGP refers to EEzGP in chapter 4, naming was changed at the last stage and the code is still under HEzGP for namings 
####

rm(list=ls())
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#RMSE
nmethod=11
nsim=30
#all
RMSE_mat2=matrix(NA,nrow=nsim,ncol=nmethod)
#EzGP
RMSE_mat2[,1]=method_EzGP_1[1:nsim]#twin
RMSE_mat2[,2]=method_EzGP_2[1:nsim]#NN
RMSE_mat2[,3]=method_EzGP_3[1:nsim]#la
RMSE_mat2[,4]=method_EzGP_4[1:nsim]#LE-1
RMSE_mat2[,5]=method_EzGP_6[1:nsim] #"Vecchia"

#HEzGP
RMSE_mat2[,6]=method_HEzGP_1[1:nsim]#twin
RMSE_mat2[,7]=method_HEzGP_2[1:nsim]#NN
RMSE_mat2[,8]=method_HEzGP_3[1:nsim]#la
RMSE_mat2[,9]=method_HEzGP_4[1:nsim]#LE-1
RMSE_mat2[,10]=method_HEzGP_6[1:nsim] #"Vecchia"
RMSE_mat2[,11]=method_HEzGP_7[1:nsim] #"SVecchia"




Methods=c(rep(c('Twin',"NN","La","LE","VA"),2),"SVA")
Models=c(rep(c("EzGP","EEzGP"), each= 5*nsim),rep("EEzGP",nsim))

#boxplot
RMSE_boxplot=NULL
for(i in 1:nmethod){
  #for(j in 1:(P)){
  RMSE_boxplot=  c(RMSE_boxplot, (RMSE_mat2[,i]))
  #}
}

df_RMSE_boxplot = data.frame(Methods=rep(Methods,each=nsim), RMSE= RMSE_boxplot, Models=Models)


library(ggplot2)
library(dplyr)
library(latex2exp)
library(ggpubr)
library("gridExtra")

RMSE_median=apply(RMSE_mat2,2,median)

RMSE_plot=ggplot(df_RMSE_boxplot, aes(x=Models, y=log(RMSE), color=Methods)) +
  geom_boxplot(size=1)+ ggtitle("(a) log(RMSE)")+
  theme(legend.position = "bottom", legend.box = "horizontal",legend.text=element_text(size = 12,face="bold"),legend.title = element_text(size=12,face = "bold"),
        axis.title.y = element_text(size = 12,face="bold"),axis.title.x = element_text(size = 12,face="bold"),
        axis.text=element_text(size = 12,face = "bold",colour="black"),plot.title = element_text(hjust = 0.5,face = "bold",size=12),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))+guides(color = guide_legend(nrow = 1))+  
  geom_hline(yintercept = log(min(RMSE_median)), linetype = "dashed", color = "darkgray", size = 1)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#time
time_secs2=round(colMeans(cbind(method_EzGP1_1[1:nsim],method_EzGP2_2[1:nsim],
                               method_EzGP3_3[1:nsim],method_EzGP4_4[1:nsim],
                               method_EzGP6_6[1:nsim],
                               method_HEzGP1_1[1:nsim],method_HEzGP2_2[1:nsim],
                               method_HEzGP3_3[1:nsim],method_HEzGP4_4[1:nsim],
                               method_HEzGP6_6[1:nsim], method_HEzGP7_7[1:nsim])),digits = 2)

Models=c(rep("EzGP",5),rep("EEzGP",6))
time_s <- data.frame(Models=Models,Time = time_secs2, Methods = Methods)
# Subset the data for each model
time_s_EzGP <- subset(time_s, Models == "EzGP")
time_s_HEzGP <- subset(time_s, Models == "EEzGP")

# Create the plot
base_plot <- ggplot(time_s, aes(x = Methods, y = Time, group = Models)) + 
  geom_point(aes(color = Methods), size = 5) +  # Points are colored by Methods
  geom_line(linetype = "dashed", color = "black") +  # Lines are black and dashed
  facet_wrap(~ Models, scales = "free_x") +  # Move facet labels to the bottom
  xlab("Methods") +
  ylab("Time (mins)") +
  ggtitle("(b) Time (mins)") +
  scale_color_discrete(name = "Methods") +
  # theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12, face = "bold", colour = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    strip.background = element_blank(),  # Optional: Remove background of facet labels
    strip.text = element_text(size = 12, face = "bold")  # Optional: Adjust text size and style
  ) +
  guides(color = guide_legend(nrow = 1))

# Create zoomed-in plot for EzGP
zoomed_plot_EzGP <- ggplot(time_s_EzGP, aes(x = Methods, y = Time, group = Models)) +
  geom_point(aes(color = Methods), size = 3) +
  geom_line(linetype = "dashed", color = "black") +
  scale_color_discrete(name = "Methods") +
  coord_cartesian(ylim = c(0.5,4.2)) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 9, face = "bold", colour = "black"),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.position = "none",
    axis.title = element_blank()
  )

# Create zoomed-in plot for HEzGP
zoomed_plot_HEzGP <- ggplot(time_s_HEzGP, aes(x = Methods, y = Time, group = Models)) +
  geom_point(aes(color = Methods), size = 3) +
  geom_line(linetype = "dashed", color = "black") +
  scale_color_discrete(name = "Methods") +
  coord_cartesian(ylim = c(0, 1.7)) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour= "black"),
    axis.text = element_text(size =9 , face = "bold", colour = "black"),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.position = "none",
    axis.title = element_blank()
  )

# Convert the zoomed plots to grobs
zoomed_grob_EzGP <- ggplotGrob(zoomed_plot_EzGP)
zoomed_grob_HEzGP <- ggplotGrob(zoomed_plot_HEzGP)

# Define positions for insets (adjust as needed)
inset_x_start_EzGP <- 3
inset_x_end_EzGP <- 5
inset_y_bottom_EzGP <- max(time_s_EzGP$Time) - 8
inset_y_top_EzGP <- max(time_s_EzGP$Time) + 1

inset_x_start_HEzGP <- 3
inset_x_end_HEzGP <- 5
inset_y_bottom_HEzGP <- max(time_s_HEzGP$Time) - 8
inset_y_top_HEzGP <- max(time_s_HEzGP$Time) + 1

# Display the base plot first
Time_plot <- base_plot

# Add the zoomed-in plots manually using grid viewports
grid.newpage()
print(Time_plot, newpage = FALSE)

# For the EzGP facet (left side)
pushViewport(viewport(x = 0.1, y = 0.85, width = 0.25, height = 0.25, just = c("left", "top")))
grid.draw(zoomed_grob_EzGP)
popViewport()

# For the HEzGP facet (right side)
pushViewport(viewport(x = 0.5, y = 0.87, width = 0.3, height = 0.2, just = c("left", "top")))
grid.draw(zoomed_grob_HEzGP)
popViewport()


filename<- paste("Boxplot_Time_E_HE.pdf",sep="")
pdf(file=filename,width =15, # The width of the plot in inches
    height =10)
ggarrange( RMSE_plot,Time_plot,ncol=1, nrow=2, common.legend = TRUE, legend="bottom")
dev.off()

save.image("SE_E.Rdata")
