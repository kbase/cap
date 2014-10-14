FROM ubuntu:14.04

RUN apt-get update && apt-get install -y \
	samtools \
	bedtools \
	python \
	git \
	r-base

RUN git clone https://github.com/kbase/cap.git
ENV PATH $PATH:/cap/scripts

RUN R --vanilla --slave < /cap/dockefiles/cap/install_packages.r

