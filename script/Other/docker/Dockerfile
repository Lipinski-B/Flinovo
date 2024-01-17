# get shiny server plus tidyverse packages image
FROM rocker/shiny-verse:latest

# system libraries of general use
RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev nano

# install R packages required escape monocle future
# (change it depending on the packages you need)
RUN R -e "install.packages(c('Seurat','plotly','shiny','shinybusy','shinydashboard','dashboardthemes','shinyWidgets'), repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('dittoSeq')"

# Copy configuration files into the Docker image
COPY shiny-server.conf  /etc/shiny-server/shiny-server.conf
COPY /app /srv/shiny-server/
RUN rm /srv/shiny-server/index.html

# Make the ShinyApp available at port 80
EXPOSE 8080
EXPOSE $PORT

# Copy further configuration files into the Docker image
COPY shiny-server.sh /usr/bin/shiny-server.sh
RUN ["chmod", "+x", "/usr/bin/shiny-server.sh"]
CMD ["/usr/bin/shiny-server.sh"]