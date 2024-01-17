library(ggplot2)
library(png)
library(grid)
library(ggimage)

bees <- data.frame(distance = c(0.5, 1.5, 
                                2.5),
                   number = c(15,40,15))

ggplot(data = bees, aes(x = distance, y = number)) +
  geom_point()

ggplot(data = bees, aes(x = distance, y = number)) +
  geom_point() +
  xlab("Distance (km)") +
  ylab("Number of Bees") +
  ylim(0, 45) + xlim(0, 3)

img <- readPNG("/home/boris/Bureau/scShiny/www/bcr.png")


ggplot(data = bees, aes(x = distance, y = number)) +
  annotation_custom(rasterGrob(img, 
                               width = unit(1,"npc"),
                               height = unit(1,"npc")), 
                    -Inf, Inf, -Inf, Inf) +
  geom_point() +
  xlab("Distance (km)") +
  ylab("Number of Bees") +
  ylim(0, 45) + xlim(0, 3)


b <- "Clonotype"

a <- as.character(paste0("Heavy : ", singlet@tools$vloupe$igh_c_genes[1],
            "\n V : ",singlet@tools$vloupe$V_lourde[1], 
            "\n D : ", singlet@tools$vloupe$D_lourde[1], 
            "\n J : ", singlet@tools$vloupe$J_lourde[1]))

c <- as.character(paste0("Light : ",
            "\n V : ", singlet@tools$vloupe$V_legere[1], 
            "\n J : ", singlet@tools$vloupe$J_legere[1]))


bees$image <- c(a,b,c)
bees$image <-  plot_ly(x = singlet@tools$Heavy[[1]], y = singlet@tools$Heavy[[2]], name = "Heavy Chain", type = "bar", hovertemplate = paste0("Locus : %{x}\nProportion : ", round(singlet@tools$Heavy[[2]]/sum(singlet@tools$Heavy[[2]]),3)*100,"%")) %>% layout(title='Frequencies Heavy Chain' , xaxis = list(tickangle = 45), yaxis =list(title="Number of cells"))



ggplot(data = bees, aes(x = distance, y = number)) +
  annotation_custom(rasterGrob(img, 
                               width = unit(1,"npc"),
                               height = unit(1,"npc")), 
                    -Inf, Inf, -Inf, Inf) +
  geom_image(aes(text = image), size = 0.15) +
  xlab("Distance (km)") +
  ylab("Number of Bees") +
  ylim(0, 45) + xlim(0, 3)




ggplot(data = bees, aes(x = distance, y = number)) +
  annotation_custom(rasterGrob(img, 
                               width = unit(1,"npc"),
                               height = unit(1,"npc")), 
                    -Inf, Inf, -Inf, Inf) +
  geom_text(size = 10,aes(label = image)) +
  xlab("Distance (km)") +
  ylab("Number of Bees") +
  ylim(0, 45) + xlim(0, 3) +
  theme(panel.background = element_rect(fill = "white"))+
  theme(axis.line.x = element_line(color = "white"),axis.line.y = element_line(color = "white"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())



load(file = paste0("/home/boris/Bureau/scShiny/datasets/", patient,".RData"))


