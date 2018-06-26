#!/usr/bin/env bash

# Builds the gis-data JAR, builds a docker image, pushes the image to Nexus and restarts the pod.

DOCKER_REGISTRY="docker-hub.informatik.haw-hamburg.de"
PROJECT="mars/mars-cuda-nn3d"
SERVICE_NAME="cuda-test-svc"

docker build -t ${DOCKER_REGISTRY}/${PROJECT}/${SERVICE_NAME}:dev .
docker push ${DOCKER_REGISTRY}/${PROJECT}/${SERVICE_NAME}:dev

# kubectl delete pod -l service=${SERVICE_NAME} --force
kubectl delete job ${SERVICE_NAME} 
kubectl apply -f cuda-test-svc.yml
#kubectl delete pod -l service=${SERVICE_NAME} --force
