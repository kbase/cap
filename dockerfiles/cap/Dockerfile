FROM ubuntu:14.04

RUN apt-get update && apt-get install -y python && mkdir /cap/

# instructions for R here

ADD ./get-rpkm.py /cap/
ADD ./merge.py /cap/
ADD ./core.R /cap/