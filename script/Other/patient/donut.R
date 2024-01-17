# Library -----------------------------------------------------------------
library(stringr)
library(plyr)
library(ggplot2)
library(plotly)
library(dplyr)

# Test donut --------------------------------------------------------------
patient="FL140304"
load(file=paste0("/home/boris/Documents/lipinskib/Boris_Manon/flinovo/result/",patient,"/R/VDJ.RData"))
data <- data.frame(clonotype=vloupe$clonotype_id,frequency=vloupe$frequency)
data$fraction = data$frequency / sum(data$frequency) # Compute percentages
data$ymax = cumsum(data$fraction) # Compute the cumulative percentages (top of each rectangle)
data$ymin = c(0, head(data$ymax, n=-1)) # Compute the bottom of each rectangle
data$labelPosition <- (data$ymax + data$ymin) / 2 # Compute label position
data %>% plot_ly(labels = ~clonotype, values = ~frequency) %>% add_pie(hole = 0.6) %>% layout(title = "Donut charts using Plotly",  showlegend = F, xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE), yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

# Final donut -------------------------------------------------------------
ax <- list(title = "",zeroline = FALSE,showline = FALSE,showticklabels =  FALSE,showgrid = FALSE)
ax2 <- list(title = "",zeroline = T,showline = F,showticklabels =  FALSE,showgrid = FALSE)

for (patient in c( "FL08G0431")) {
  ## -- Donut -- ##
  load(file=paste0("/home/boris/Documents/lipinskib/Boris_Manon/flinovo/result/",patient,"/R/VDJ.RData"))

  if(length(Light[[2]])==3){
    Light[[1]]<- c(Light[[1]][1],Light[[1]][3])
    Light[[2]]<- c(Light[[2]][1],Light[[2]][3])
  }

  p1 <- plot_ly(x = Heavy[[1]], y = Heavy[[2]], name = "Heavy", type = "bar", showlegend = FALSE,hovertemplate = paste0("Locus : %{x}\nProportion : ", round(Heavy[[2]]/sum(Heavy[[2]]),3)*100,"%")) %>% layout(xaxis= list(showticklabels = T, tickangle = 45), yaxis = ax2)
  p2 <- plot_ly(x = Light[[1]], y = Light[[2]], name = "Light", type = "bar", showlegend = FALSE,hovertemplate = paste0("Locus : %{x}\nProportion : ", round(Light[[2]]/sum(Light[[2]]),3)*100,"%")) %>% layout(xaxis= list(showticklabels = T, tickangle = 45), yaxis = ax2)#plot_ly(x = D[[1]], y = D[[2]], name = "D", type = "bar", showlegend = T) %>% layout(xaxis= list(showticklabels = FALSE)),
  p3 <- plot_ly(x = V[[1]], y = V[[2]], name = "V", type = "bar", showlegend = F) %>% layout(xaxis= list(showticklabels = FALSE), yaxis = ax2)
  p4 <- plot_ly(x = J[[1]], y = J[[2]], name = "J", type = "bar", showlegend = F) %>% layout(xaxis= list(showticklabels = FALSE), yaxis = ax2)

  e <- plot_ly(hole = 0.9 ,type = "pie",labels = vloupe$clonotype_id, values = vloupe$frequency, showlegend = FALSE,textposition = 'side',textinfo = 'percent', hovertemplate = paste0('Clonotype : ', vloupe$clonotype_id, "\nProportion : ", round(vloupe$proportion,3)*100, "% \nType : ", vloupe$type, "\nIsotype : ", vloupe$igh_c_genes,"\nHeavy : ", vloupe$V_lourde, " / ", vloupe$D_lourde, " / ", vloupe$J_lourde, "\nLight : ", vloupe$V_legere, " / ", vloupe$J_legere), domain = list(x = c(0, 1), y = c(0, 1)) ) %>% layout(title = patient, autosize = T) %>%
        subplot(plot_ly() %>% layout(xaxis = ax, yaxis = ax), p1, p2 ,p3, p4, nrows = 4, shareX = F, shareY = F, margin = 0.04, heights = c(0,0.5,0.25,0.05), widths = c(0.2,0.2))

  ## -- Sunburst -- ##
  load(file = paste0("/home/boris/Bureau/Flinovo/result/analyse_patient/",patient,"/All_All_",patient,".RData"))

  f <-plot_ly(singlet@tools$sunburst, ids = ~ids ,labels = ~labels, parents = ~parents, values = ~values, marker = list(colors = c( "#BEBADA", "#8DD3C7", "#FB8072", "#80B1D3", "#FDB462")),type = 'sunburst', branchvalues = 'total', hoverinfo = "text", hovertext = paste(singlet@tools$sunburst$labels, ":", round((singlet@tools$sunburst$values/singlet@tools$sunburst$values[1])*100,2),"%", "\nTotal : " ,singlet@tools$sunburst$values)) 
  
  save(e, f, file = paste0("/home/boris/Bureau/scShiny/inst/app/www/",patient,"/",patient,"_VDJ.RData"))
}
