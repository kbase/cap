FROM ubuntu:14.04

RUN apt-get update && apt-get install -y \
	samtools \
	bedtools \
	python \
	r-base && \
	mkdir /cap/

ADD ./install_packages.r /cap/
ADD ./get-rpkm.py /cap/
ADD ./merge.py /cap/
ADD ./core.R /cap/

RUN R --vanilla --slave < /cap/install_packages.r