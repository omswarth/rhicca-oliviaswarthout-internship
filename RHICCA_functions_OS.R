
#function to create a histogram
histFun <- function(dat, var, groupVar, log=FALSE, xLim=NULL, xLab = "", title=NA,
                    bin=50, colors=c("grey", "#00AFBB", "#E7B800")){
  
  #list of detection limits for biomarkers, choose proper one
  limits <- c("ICAM1"=250000, "CD163"=1000, "ILb"=2.44140625, "IL6"=2.44140625)
  limit <- limits[var]
  xaes <- dat[,var]
  
  #if variable is to be log transformed, transform data as well as detection limit
  if(log==TRUE){
    xaes <- log(xaes)
    limit <- log(limit)
  }
  
  #set dimensions for shaded boxes that show detection limit  
  low <- ifelse(var=="ICAM1"||var=="CD163", limit, -Inf)
  high <- ifelse(var=="ICAM1"||var=="CD163", Inf, limit)
  
  #main function
  #one geom_hist for the overall data, one split into different colors by groupVar
  #data and aes can't go directly in ggplot() or else shading box won't be drawn
  #geom_vline and geom_rect draw the  detection limit bar and box, alpha controls transparency
  #scale_color and scale_fill manual define colors
  #theme changes appearance 
  ggplot()+
    geom_histogram(data=dat, mapping=aes(x=xaes, color="ALL", fill="All"), bins=bin) +
    geom_histogram(data=dat, aes(x=xaes, color=dat[,groupVar], fill=dat[,groupVar]), 
                   position = "identity", bins=bin, alpha=0.6)+
    geom_vline(aes(xintercept = limit), color="red", size=0.5, alpha=0.4) +
    geom_rect(aes(xmin=low,xmax=high, ymin=0,ymax=Inf),alpha=0.1,fill="red") +
    scale_color_manual(guide=FALSE, values = colors) +
    scale_fill_manual(" ",values = colors)+
    xlab(xLab) +
    scale_x_continuous(labels = comma, limits = xLim)+
    theme(panel.grid = element_line(linetype = "dotted", color = "grey"), 
          panel.background = element_rect(fill = "white"),
          legend.key.size = unit(0.03, "npc"),
          legend.position = "top")
}

#function to create violin plot
violinFun<-function(dat,var,groupVar,log=FALSE,yLim=NULL,yLab=NULL,labels=NULL,title=NA,
                    fontSize=10,colors=c("#E7B800", "#00AFBB")){
  #if labels aren't defined, assign labels to be categories of groupVar 
  if(is.null(labels)){labels<-levels(dat[,groupVar])}
  if(is.null(yLab)){yLab<-var}
  
  #set detection limit   
  limits <- c("ICAM1"=250000, "CD163"=1000, "ILb"=2.44140625, "IL6"=2.44140625)
  limit <- limits[var]
  yaes <- dat[,var]
  
  #log transform data if log = true
  if(log==TRUE){
    yaes <- log(yaes)
    limit <- log(limit)
  }
  
  #main function
  #no shading here so we can put data directly in ggplot()
  #geom_boxplot adds boxplot on top
  #use predefined theme
  ggplot(data=dat,mapping=aes(x=get(groupVar),y=yaes,fill=get(groupVar))) +
    geom_violin(show.legend=F) +
    scale_fill_manual(values=colors) +
    geom_boxplot(width=0.05, fill="white") +
    geom_hline(yintercept=limit, color="red", size=0.5, alpha=0.4)+
    xlab("") +
    ylab(yLab) +
    theme_pubclean(base_size=fontSize) +
    scale_x_discrete(breaks=levels(dat[,groupVar]),labels=labels)+
    scale_y_continuous(labels = comma, limits = yLim, n.breaks=7)
}

#function to print two given plots side by side
#create one viewport anchored at left side of page and one anchored in middle
printPlots <- function(a, b){
  
  vp.Right<-viewport(height=unit(1, "npc"), width=unit(0.5, "npc"),just=c("left","top"),y=1, x=0.5)
  vp.Left<-viewport(height=unit(1, "npc"), width=unit(0.5, "npc"),just=c("left","top"),y=1, x=0)
  
  print(a, vp=vp.Left)
  print(b, vp=vp.Right)
}

#print histogram split by a chosen facet (categorical  variable)
#in my case this was age group 
facetsFun <- function(dat, var, groupVar, xLim, facets, log=TRUE, colors=c("grey", "#00AFBB", "#E7B800"), xLab=""){
  xaes = dat[,var]
  if(log==TRUE){
    xaes <- log(xaes)
  }

  #same as above, except for facet_wrap
  #facet_wrap simply divides the data into different plots depending on variable which was passed
  ggplot()+
    geom_histogram(data=dat, mapping=aes(x=xaes, color="ALL", fill="All")) +
    geom_histogram(data=dat, aes(x=log(dat[,var]), color=dat[,groupVar], fill=dat[,groupVar]), position = "identity",   alpha=0.6)+
    scale_x_continuous( limits = xLim)+
    xlab(xLab) +
    scale_color_manual(guide=FALSE, values = colors) +
    scale_fill_manual(" ",values = colors)+
    theme(panel.grid = element_line(linetype = "dotted", color = "grey"), 
          strip.background = element_rect(fill="white"),
          panel.background = element_rect(fill = "white"),
          legend.key.size = unit(0.03, "npc"),
          legend.position = "top")+
    facet_wrap(~ dat[,facets])
}

