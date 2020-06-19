FROM r-base:3.6.2

ADD bin/Plot.R /opt/Plot.R

# Install R packages
RUN Rscript -e "install.packages('rMVP')"