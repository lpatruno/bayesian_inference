FROM rocker/rstudio:3.4.3

RUN sudo apt-get update
RUN sudo apt-get install -y jags

RUN Rscript -e "install.packages('rjags')"
RUN Rscript -e "install.packages('car')"