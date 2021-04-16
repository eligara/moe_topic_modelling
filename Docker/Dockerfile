#script to launch simulation on OCCAM facility

# Base image https://hub.docker.com/u/rocker/
FROM rocker/r-base:latest

#RUN apt-get update

RUN useradd -m -g users elisa   

#set workdirectory 
WORKDIR /home/elisa/ 

RUN mkdir /home/elisa/files

ADD myScript.R /home/elisa/

COPY install_packages.R /home/elisa/
COPY countsORs_3D.zip /home/elisa/files/

#run requirements script

RUN chmod a+x /home/elisa/myScript.R
RUN Rscript /home/elisa/install_packages.R



USER elisa
ENTRYPOINT ["Rscript", "/home/elisa/myScript.R"]
