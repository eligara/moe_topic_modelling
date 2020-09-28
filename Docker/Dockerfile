# Base image https://hub.docker.com/u/rocker/
FROM rocker/r-base:latest

#RUN apt-get update

RUN useradd -m -g users elisa   

#set workdirectory 
WORKDIR /home/elisa/ 

#install packages
#RUN apt-get -t testing -y install libcurl4-openssl-dev
#RUN apt-get -t testing -y install libssl-dev   

RUN mkdir /home/elisa/files

ADD myScript.R /home/elisa/

#COPY requirements-bin.txt /home/elisa/
COPY install_packages.R /home/elisa/
COPY countsORs_3D.zip /home/elisa/files/

#run requirements script

#RUN cat /home/elisa/requirements-bin.txt | xargs apt-get install -y -qq


RUN chmod a+x /home/elisa/myScript.R
RUN Rscript /home/elisa/install_packages.R

#COPY . /home/elisa/

USER elisa
ENTRYPOINT ["Rscript", "/home/elisa/myScript.R"]
#CMD ["2", "10"] #argomento letto da sys.argv in run.py