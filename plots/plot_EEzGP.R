####
#important note: HEzGP refers to EEzGP in chapter 4, naming was changed at the last stage and the code is still under HEzGP for namings 
####
rm(list=ls())
#compare diff EzGP
nsim=30

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#HEzGP
method_HEzGP_1=c()
method_HEzGP1_1=c()

method_HEzGP_2=c()
method_HEzGP2_2=c()

method_HEzGP_3=c()
method_HEzGP3_3=c()

method_HEzGP_4=c()
method_HEzGP4_4=c()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nsim_start=1
nsim_end=30
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#HEZGP
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#1
model_name="Twin-HEzGP-relabel"
filename= paste0("RMSE","_",model_name,"_",".rds",sep="")
method_HEzGPsep_1=readRDS(filename)
method_HEzGP_1=c(method_HEzGP_1,method_HEzGPsep_1[nsim_start:nsim_end])
filename= paste0("time","_",model_name,"_",".rds",sep="")
method_HEzGPsep1_1=as.numeric(readRDS(filename))
method_HEzGP1_1=c(method_HEzGP1_1,method_HEzGPsep1_1[nsim_start:nsim_end])

#2  
model_name="NN-HEzGP-relabel"
filename= paste0("RMSE","_",model_name,"_",".rds",sep="")
method_HEzGPsep_2=readRDS(filename)
method_HEzGP_2=c(method_HEzGP_2,method_HEzGPsep_2[nsim_start:nsim_end])
filename= paste0("time","_",model_name,"_",".rds",sep="")
method_HEzGPsep2_2=as.numeric(readRDS(filename))
method_HEzGP2_2=c(method_HEzGP2_2,method_HEzGPsep2_2[nsim_start:nsim_end])

#3
model_name="La-alc-HEzGP-relabel"
filename= paste0("RMSE","_",model_name,"_",".rds",sep="")
method_HEzGPsep_3=readRDS(filename)
method_HEzGP_3=c(method_HEzGP_3,method_HEzGPsep_3[nsim_start:nsim_end])
filename= paste0("time","_",model_name,"_",".rds",sep="")
method_HEzGPsep3_3=as.numeric(readRDS(filename))
method_HEzGP3_3=c(method_HEzGP3_3,method_HEzGPsep3_3[nsim_start:nsim_end])

#4
model_name="LHEzGP-ns-1-relabel"
filename= paste0("RMSE","_",model_name,"_",".rds",sep="")
method_HEzGPsep_4=readRDS(filename)
method_HEzGP_4=c(method_HEzGP_4,method_HEzGPsep_4[nsim_start:nsim_end])
filename= paste0("time","_",model_name,"_",".rds",sep="")
method_HEzGPsep4_4=as.numeric(readRDS(filename))
method_HEzGP4_4=c(method_HEzGP4_4,method_HEzGPsep4_4[nsim_start:nsim_end])

#5
model_name="HEz-Vecchia_encode-1"
filename= paste0("RMSE","_",model_name,"_",".rds",sep="")
method_HEzGP_6=readRDS(filename)
filename= paste0("time","_",model_name,"_",".rds",sep="")
method_HEzGP6_6=as.numeric(readRDS(filename))
#6
model_name="HEz-SVecchia_encode-1"
filename= paste0("RMSE","_",model_name,"_",".rds",sep="")
method_HEzGP_7=readRDS(filename)
filename= paste0("time","_",model_name,"_",".rds",sep="")
method_HEzGP7_7=as.numeric(readRDS(filename))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#RMSE
nmethod=6
nsim=30

RMSE_mat=matrix(NA, nrow=nsim,ncol=nmethod)
RMSE_mat[,1]=method_HEzGP_1[1:nsim]#twin
RMSE_mat[,2]=method_HEzGP_2[1:nsim]#NN
RMSE_mat[,3]=method_HEzGP_3[1:nsim]#la
RMSE_mat[,4]=method_HEzGP_4[1:nsim]#LE-1
RMSE_mat[,5]=method_HEzGP_6[1:nsim] #"Vecchia"
RMSE_mat[,6]=method_HEzGP_7[1:nsim] #"SVecchia"




Methods=c('Twin',"NN","La","LE","VA","SVA")
Models=c(rep(c("SEzGP"), each= 6*nsim))

#boxplot
RMSE_boxplot=NULL
for(i in 1:nmethod){
  RMSE_boxplot=  c(RMSE_boxplot, (RMSE_mat[,i]))
}

df_RMSE_boxplot = data.frame(Methods=rep(Methods,each=nsim), RMSE= RMSE_boxplot, Models=Models)


library(ggplot2)
library(dplyr)
library(latex2exp)
library(ggpubr)
library("gridExtra")

library(ggplot2)
library(grid)

RMSE_plot=ggplot(df_RMSE_boxplot, aes(x=Methods, y=log(RMSE),color=Methods)) +
  geom_boxplot(size=1)+ ggtitle("(a) log(RMSE)")+
  theme(legend.position = "bottom", legend.box = "horizontal",legend.text=element_text(size=12,face="bold"),legend.title = element_text(size=12,face = "bold"),
        axis.title.y = element_text(size=12,face="bold"), axis.title.x =element_blank(),
        axis.text=element_text(size=12,face = "bold",colour="black"),plot.title = element_text(hjust = 0.5,face = "bold",size=12),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))+guides(color = guide_legend(nrow = 1))


RMSE_median=apply(RMSE_mat,2,median)

RMSE_plot=ggplot(df_RMSE_boxplot, aes(x=Methods, y=log(RMSE),color=Methods)) +
  geom_boxplot(size=1)+ ggtitle("(a) log(RMSE)")+
  theme(legend.position = "bottom", legend.box = "horizontal",legend.text=element_text(size=12,face="bold"),legend.title = element_text(size=12,face = "bold"),
        axis.title.y = element_text(size=12,face="bold"), axis.title.x =element_blank(),
        axis.text=element_text(size=12,face = "bold",colour="black"),plot.title = element_text(hjust = 0.5,face = "bold",size=12),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))+guides(color = guide_legend(nrow = 1))+
  geom_hline(yintercept = log(min(RMSE_median)), linetype = "dashed", color = "darkgray", size = 1)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#time
time_secs=round(colMeans(cbind(method_HEzGP1_1[1:nsim],method_HEzGP2_2[1:nsim],
                               method_HEzGP3_3[1:nsim],method_HEzGP4_4[1:nsim],
                                method_HEzGP6_6[1:nsim],method_HEzGP7_7[1:nsim])),digits = 2)

Models=c(rep(c("SEzGP"), each= 6))
time_s <- data.frame(Models=Models,Time = time_secs, Methods = Methods)



# Main plot
base_plot <- ggplot(time_s, aes(x = Methods, y = Time, group = Models)) + 
  geom_point(aes(color = Methods), size = 5) +  
  geom_line(linetype = "dashed", color = "black") +  
  facet_wrap(~ Models, scales = "free_x") +  
  ylab("Time (mins)") + 
  ggtitle("(b) Time (mins)") + 
  scale_color_discrete(name = "Methods") + 
  coord_cartesian(ylim = c(0, max(time_s$Time)+1)) +  
  theme(
    legend.position = "bottom",  
    legend.box = "horizontal",
    legend.text = element_text(size=12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size=12, face = "bold"),
    #axis.title.x = element_text(size=12, face = "bold"),
    axis.title.x =element_blank(),
    axis.text = element_text(size=12, face = "bold", colour = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    strip.background = element_blank(),
    strip.text = element_blank()
  ) +
  guides(color = guide_legend(nrow = 1))

# Zoomed-in plot
zoomed_plot <- ggplot(time_s, aes(x = Methods, y = Time, group = Models)) + 
  geom_point(aes(color = Methods), size = 3) +  
  geom_line(linetype = "dashed", color = "black") +  
  facet_wrap(~ Models, scales = "free_x") +  
  scale_color_discrete(name = "Methods") + 
  coord_cartesian(ylim = c(0,2)) +  
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 8, face = "bold", colour = "black"),
        strip.background = element_blank(),  # Remove background of facet labels
        strip.text = element_blank(),        # Hide the facet label text
        legend.position = "none",            # Hide the legend
        axis.title = element_blank()
  )

# Convert zoomed plot to a grob
zoomed_grob <- ggplotGrob(zoomed_plot)

# Define approximate positions for the inset
# For example, let's place it in the top-right corner
inset_x_start <-3  # Adjust as needed
inset_x_end <- 6  # Adjust as needed
inset_y_bottom <- max(time_s$Time) -7# Place inset just below the max y-axis value
inset_y_top <- max(time_s$Time)+1         # Adjust as needed to fit the inset plot

# Add the zoomed-in plot as an inset
Time_plot=base_plot +
  annotation_custom(
    grob = zoomed_grob,
    xmin = inset_x_start,
    xmax = inset_x_end,
    ymin = inset_y_bottom,
    ymax = inset_y_top
  )



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ms
#1
model_name="HEz-Vecchia_encode-ms2"
filename= paste0("RMSE","_",model_name,"_",".rds",sep="")
Vecchia_RMSE=readRDS(filename)
filename= paste0("time","_",model_name,"_",".rds",sep="")
Vecchia_time=as.numeric(readRDS(filename))
#2
model_name="HEz-SVecchia_encode-ms2"
filename= paste0("RMSE","_",model_name,"_",".rds",sep="")
SVecchia_RMSE=readRDS(filename)
filename= paste0("time","_",model_name,"_",".rds",sep="")
SVecchia_time=as.numeric(readRDS(filename))

Methods=c("VA","SVA")

RMSE_vec=cbind(Vecchia_RMSE,SVecchia_RMSE)
Time_vec=cbind(Vecchia_time,SVecchia_time)


RMSE_Vec=NULL
Time_Vec=NULL

for(i in 1:2){
  RMSE_Vec=  c(RMSE_Vec, (RMSE_vec[,i]))
  Time_Vec=  c(Time_Vec, (Time_vec[,i]))
}


df_vec = data.frame(Methods=rep(Methods,each=8), RMSE= RMSE_Vec,Time=Time_Vec,ms=1:8)



# Main plot
ms_RMSE <- ggplot(df_vec, aes(x = ms, y = log(RMSE), group = Methods)) + 
  geom_point(aes(color = Methods), size = 5) +  
  geom_line(linetype = "dashed", color = "black") +  
  ylab("log(RMSE)") + 
  xlab(expression(bold(m[s])))+
  ggtitle("(a) log(RMSE)") + 
  scale_color_discrete(name = "Methods") + 
  theme(
    legend.position = "bottom",  
    legend.box = "horizontal",
    legend.text = element_text(size=12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size=12, face = "bold"),
    axis.text = element_text(size=12, face = "bold", colour = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    strip.background = element_blank(),
    strip.text = element_blank()
  ) +
  guides(color = guide_legend(nrow = 1))

ms_time<- ggplot(df_vec, aes(x = ms, y = Time, group = Methods)) + 
  geom_point(aes(color = Methods), size = 5) +  
  geom_line(linetype = "dashed", color = "black") +  
  ylab("Time (min)") + 
  xlab(expression(bold(m[s])))+
  ggtitle("(b) Time (min)") + 
  scale_color_discrete(name = "Methods") + 
  theme(
    legend.position = "bottom",  
    legend.box = "horizontal",
    legend.text = element_text(size=12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size=12, face = "bold"),
    axis.text = element_text(size=12, face = "bold", colour = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    strip.background = element_blank(),
    strip.text = element_blank()
  ) +
  guides(color = guide_legend(nrow = 1))

filename<- paste("HEzGP_ms.pdf",sep="")
pdf(file=filename,width =10, # The width of the plot in inches
    height =5)
ggarrange( ms_RMSE,ms_time,ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
dev.off()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ms
#1
model_name="NN-HEzGP-mixed"
filename= paste0("RMSE","_",model_name,"_",".rds",sep="")
NN_mixed_RMSE=readRDS(filename)
filename= paste0("time","_",model_name,"_",".rds",sep="")
NN_mixed_time=as.numeric(readRDS(filename))
#2
model_name="NN-HEzGP-relabel"
filename= paste0("RMSE","_",model_name,"_",".rds",sep="")
NN_relabel_RMSE=readRDS(filename)
filename= paste0("time","_",model_name,"_",".rds",sep="")
NN_relabel_time=as.numeric(readRDS(filename))

Methods=c("Gower","one-hot")

RMSE_dis=cbind(NN_mixed_RMSE,NN_relabel_RMSE)
Time_dis=cbind(NN_mixed_time,NN_relabel_time)


RMSE_dis2=NULL
Time_dis2=NULL
for(i in 1:2){
  RMSE_dis2=  c(RMSE_dis2, (RMSE_dis[,i]))
  Time_dis2=  c(Time_dis2, (Time_dis[,i]))
  
}


df_dis = data.frame(Methods=rep(Methods,each=nsim), RMSE= RMSE_dis2,Time=Time_dis2)

RMSE_median_dis=apply(RMSE_dis,2,median)
RMSE_plot_dis=ggplot(df_dis, aes(x=Methods, y=log(RMSE),color=Methods)) +
  geom_boxplot(size=1)+ ggtitle("(a) log(RMSE)")+  ylab("log(RMSE)") + 
  theme(legend.position = "bottom", legend.box = "horizontal",legend.text=element_text(size=12,face="bold"),legend.title = element_text(size=12,face = "bold"),
        axis.title.y = element_text(size = 12,face="bold"), axis.title.x =element_blank(),
        axis.text=element_text(size=12,face = "bold",colour="black"),plot.title = element_text(hjust = 0.5,face = "bold",size=12),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))+guides(color = guide_legend(nrow = 1))+
  geom_hline(yintercept = log(min(RMSE_median_dis)), linetype = "dashed", color = "darkgray", size = 1)

Time_plot_dis<- ggplot(df_dis, aes(x=Methods, y=Time,color=Methods)) +
  geom_boxplot(size=1)+ ggtitle("(b) Time (min)")+  ylab("Time (min)") + 
  theme(legend.position = "bottom", legend.box = "horizontal",legend.text=element_text(size=12,face="bold"),legend.title = element_text(size=12,face = "bold"),
        axis.title.y = element_text(size=12,face="bold"), axis.title.x =element_blank(),
        axis.text=element_text(size=12,face = "bold",colour="black"),plot.title = element_text(hjust = 0.5,face = "bold",size=12),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))+guides(color = guide_legend(nrow = 1))

filename<- paste("Boxplot_Time_SEzGP.pdf",sep="")
pdf(file=filename,width =10, # The width of the plot in inches
    height =5)
ggarrange( RMSE_plot,Time_plot,ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
dev.off()

save.image("HEzGP.Rdata")
