FROM nexus.informatik.haw-hamburg.de/nvidia/cuda:9.0-devel

ADD common.h /app/common.h
ADD simpleDeviceQuery.cu /app/simpleDeviceQuery.cu

ADD geojson/gis_vector_rivers.geojson /app/gis_vector_rivers.geojson
ADD geojson/points.geojson /app/points.geojson

ADD jsmn.h /app/jsmn.h
ADD jsmn.c /app/jsmn.c
ADD kdtree.h /app/kdtree.h
ADD kdtree.c /app/kdtree.c
ADD pq.cuh /app/pq.cuh
ADD pq.cu /app/pq.cu
ADD nn3d.cu /app/nn3d.cu

WORKDIR app

#RUN apt-get update && apt-get install -y \
#    valgrind \
#	vim

RUN nvcc simpleDeviceQuery.cu -o sdq
RUN nvcc --relocatable-device-code true -g -G jsmn.c kdtree.c pq.cu nn3d.cu -o nn3d

ENTRYPOINT ./sdq 
#ENTRYPOINT 
#ENTRYPOINT ./nn3d gis_vector_rivers.geojson 10 