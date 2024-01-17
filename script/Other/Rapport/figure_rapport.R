library(Seurat)
library(monocle)
library(escape)
library(dittoSeq)
library(plotly)
library(dplyr)
library(shiny)
library(shinyjs)
library(shinybusy)
library(shinyWidgets)
library(shinydashboard)
library(dashboardthemes)
library(flexdashboard)
library(processx)
library(rmarkdown)
Sys.setenv("PATH" = "/home/boris/Software/miniconda3/bin")

setwd("/")
siege <- c("FL140304","FL12C1888","FL09C1164","FL08G0293","FL02G095","FL05G0330")
#siege <- c("FL1085",  "FL120316",  "FL1214",  "FL1481")
patient <- siege[6]

load(file = paste0("/home/boris/Documents/analyse/singlet_", patient,".RData"))

# Métadonnées : 
fig <- plot_ly(singlet@tools$sunburst, ids = ~ids ,labels = ~labels, parents = ~parents, values = ~values, marker = list(colors = c( "#BEBADA", "#8DD3C7", "#FB8072", "#80B1D3", "#FDB462")), type = 'sunburst', branchvalues = 'total', hoverinfo = "text", hovertext = paste(singlet@tools$sunburst$labels, ":", round((singlet@tools$sunburst$values/singlet@tools$sunburst$values[1])*100,2),"%", "\nTotal : " ,singlet@tools$sunburst$values))
orca(fig, paste0("/home/boris/Bureau/scShiny/www/", patient, "/metadata_", patient,".png"), more_args = c('--disable-gpu'))

# Réduction de dimension :
plots <- PCAPlot(object = singlet, group.by = "Phénotype", split.by = NULL, label.size = 0.0, pt.size = 1) & theme(title = element_text(size=20),
              legend.position = "top",
              legend.title = element_text(size=10),
              legend.text = element_text(size=10)
) & guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 2))) & xlab(label = paste0("PCA 1 : ", round(Stdev(singlet[["pca"]])[1],2), " %")) & ylab(label = paste0("PCA 2 : ", round(Stdev(singlet[["pca"]])[2],2), " %"))
orca(plots, paste0("/home/boris/Bureau/scShiny/www/", patient, "/ACP_", patient,".png"), more_args = c('--disable-gpu'))


# VDJ : 
load(file=paste0("/home/boris/Documents/lipinskib/Boris_Manon/flinovo/result/",patient,"/R/VDJ.RData"))
e <- plot_ly(hole = 0.85 ,type = "pie",labels = vloupe$clonotype_id, values = vloupe$frequency, showlegend = FALSE,
             textposition = 'side',textinfo = 'percent',
             hovertemplate = paste0('Clonotype : ', vloupe$clonotype_id, "\nProportion : ", round(vloupe$proportion,3)*100, "% \nType : ", vloupe$type, "\nIsotype : ", vloupe$igh_c_genes,"\nHeavy : ", vloupe$V_lourde, " / ", vloupe$D_lourde, " / ", vloupe$J_lourde, "\nLight : ", vloupe$V_legere, " / ", vloupe$J_legere),
             domain = list(x = c(0, 1), y = c(0, 1)) ) %>% layout(title = patient, autosize = T) %>%
  subplot(nrows = 3, shareX = F, shareY = F, margin = 0.04, heights = c(0.2,0.2,0.2), widths = c(0.2,0.2),
          plot_ly(x = V[[1]], y = V[[2]], name = "V", type = "bar", showlegend = T) %>% layout(xaxis= list(showticklabels = FALSE)),
          plot_ly(x = Heavy[[1]], y = Heavy[[2]], name = "Heavy", type = "bar", showlegend = FALSE,
                  hovertemplate = paste0("Locus : %{x}\nProportion : ", round(Heavy[[2]]/sum(Heavy[[2]]),3)*100,"%")) %>% layout(xaxis= list(showticklabels = T, tickangle = 45)),
          plot_ly(x = Light[[1]], y = Light[[2]], name = "Light", type = "bar", showlegend = FALSE,
                  hovertemplate = paste0("Locus : %{x}\nProportion : ", round(Light[[2]]/sum(Light[[2]]),3)*100,"%")) %>% layout(xaxis= list(showticklabels = T, tickangle = 45)),
          plot_ly(x = D[[1]], y = D[[2]], name = "D", type = "bar", showlegend = T) %>% layout(xaxis= list(showticklabels = FALSE)),
          plot_ly(x = J[[1]], y = J[[2]], name = "J", type = "bar", showlegend = T) %>% layout(xaxis= list(showticklabels = FALSE))
  )

orca(e, paste0("/home/boris/Bureau/scShiny/www/", patient, "/VDJ_", patient,".png"), more_args = c('--disable-gpu'))




## Clonotype : 
## -- Clonotype -- ##
fig <- plot_ly(x = singlet@tools$Clonotype[[1]], y = singlet@tools$Clonotype[[2]], name = "Info", type = "bar", textposition="inside",
               text = paste0("Proportion : ", round(singlet@tools$vloupe$proportion[1:5],3)*100,
                             "% \nType : ", singlet@tools$vloupe$type[1:5], 
                             "\nIsotype : ", singlet@tools$vloupe$igh_c_genes[1:5], 
                             "\n\nHeavy : \n", singlet@tools$vloupe$V_lourde[1:5], " / \n", singlet@tools$vloupe$D_lourde[1:5], " / \n", singlet@tools$vloupe$J_lourde[1:5], 
                             "\n\nLight : \n", singlet@tools$vloupe$V_legere[1:5], " / \n", singlet@tools$vloupe$J_legere[1:5])) %>% layout(title='Frequencies of the mains clonotypes', yaxis =list(title="Number of cells"))
       
fig
orca(fig, paste0("/home/boris/Bureau/scShiny/www/", patient, "/clonotype_", patient,".png"), more_args = c('--disable-gpu'))


## VDJ : 
## -- VDJ -- ##
fig <- plot_ly(x = singlet@tools$V[[1]], y = singlet@tools$V[[2]], name = "Clonotype", type = "bar") %>% 
  layout(title='Frequencies V genes : Heavy and Lights chains', xaxis = list(tickangle = 45), yaxis =list(title="Number of cells"))
fig
orca(fig, paste0("/home/boris/Bureau/scShiny/www/", patient, "/V_", patient,".png"), more_args = c('--disable-gpu'))


fig <- plot_ly(x = singlet@tools$D[[1]], y = singlet@tools$D[[2]], name = "Clonotype", type = "bar") %>% 
  layout(title='Frequencies D genes : Heavy chain',xaxis = list(tickangle = 45), yaxis =list(title="Number of cells"))
fig
orca(fig, paste0("/home/boris/Bureau/scShiny/www/", patient, "/D_", patient,".png"), more_args = c('--disable-gpu'))


fig <- plot_ly(x = singlet@tools$J[[1]], y = singlet@tools$J[[2]], name = "Clonotype", type = "bar") %>% 
  layout(title='Frequencies J genes : Heavy and Lights chains',xaxis = list(tickangle = 45), yaxis =list(title="Number of cells"))
fig
orca(fig, paste0("/home/boris/Bureau/scShiny/www/", patient, "/J_", patient,".png"), more_args = c('--disable-gpu'))



## Heavy : 
## -- Heavy chain -- ##
fig <- plot_ly(x = singlet@tools$Heavy[[1]], y = singlet@tools$Heavy[[2]], name = "Heavy Chain", type = "bar",
        hovertemplate = paste0("Locus : %{x}\nProportion : ", round(singlet@tools$Heavy[[2]]/sum(singlet@tools$Heavy[[2]]),3)*100,"%")) %>%  
  layout(title='Frequencies Heavy Chain' , xaxis = list(tickangle = 45), yaxis =list(title="Number of cells"))
fig
orca(fig, paste0("/home/boris/Bureau/scShiny/www/", patient, "/Heavy_", patient,".png"), more_args = c('--disable-gpu'))


fig <- plot_ly(x = singlet@tools$Isotype[[1]], y = singlet@tools$Isotype[[2]], name = "Heavy Chain", type = "bar",
        hovertemplate = paste0("Locus : %{x}\nProportion : ", round(singlet@tools$Isotype[[2]]/sum(singlet@tools$Isotype[[2]]),3)*100,"%")) %>%  
  layout(title='Frequencies Heavy Chain with details' , xaxis = list(tickangle = 45), yaxis =list(title="Number of cells"))
fig
orca(fig, paste0("/home/boris/Bureau/scShiny/www/", patient, "/Isotype_", patient,".png"), more_args = c('--disable-gpu'))



## Light : 
## -- Light chain -- ##
fig <- plot_ly(x = singlet@tools$Light[[1]], y = singlet@tools$Light[[2]], name = "Light Chain", type = "bar",
        hovertemplate = paste0("Locus : %{x}\nProportion : ", round(singlet@tools$Light[[2]]/sum(singlet@tools$Light[[2]]),3)*100,"%")) %>% 
  layout(title='Frequencies Light Chain',yaxis =list(title="Number of cells"))
fig
orca(fig, paste0("/home/boris/Bureau/scShiny/www/", patient, "/Light_", patient,".png"), more_args = c('--disable-gpu'))


fig <- plot_ly(x = singlet@tools$Type[[1]], y = singlet@tools$Type[[2]], name = "Clonotype", type = "bar",
        hovertemplate = paste0("Locus : %{x}\nProportion : ", round(singlet@tools$Type[[2]]/sum(singlet@tools$Type[[2]]),3)*100,"%")) %>% 
  layout(title='Frequencies Light Chain with details', yaxis =list(title="Number of cells"))
fig
orca(fig, paste0("/home/boris/Bureau/scShiny/www/", patient, "/Type_", patient,".png"), more_args = c('--disable-gpu'))



for (patient in siege) {
  output_dir <- paste0("/home/boris/Bureau/scShiny/www/",patient)
  render("/home/boris/Bureau/scShiny/document/script/rapport/rapport.Rmd", output_file = c(paste0("rapport_", patient),paste0("rapport_", patient),paste0("rapport_", patient)), output_format=c("html_document","pdf_document","word_document"), output_dir = output_dir, params = list(output_dir = output_dir, patient = patient))
}
