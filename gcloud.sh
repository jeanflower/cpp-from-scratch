#!/bin/bash

set -e # exit if a command returns a non-zero exit status

PROJECT_ID="cpp-project-449415"
PROJECT_NUMBER=236401783700
REGION="us-central1"
SERVICE_NAME="cloud-run-cpp"

# Uncomment these lines to refresh or change GCP credentials
#gcloud auth login
#gcloud auth application-default login

#gcloud projects get-iam-policy $PROJECT_ID --flatten="bindings[].members" --format="table(bindings.role)"

gcloud config set project $PROJECT_ID

#gcloud projects add-iam-policy-binding \
#    $PROJECT_ID \
#    --member=user:$(gcloud config get-value account) \
#    --role=roles/serviceusage.serviceUsageConsumer
gcloud auth application-default set-quota-project $PROJECT_ID

#gcloud projects add-iam-policy-binding $PROJECT_ID \
#    --member=serviceAccount:$PROJECT_NUMBER-compute@developer.gserviceaccount.com \
#    --role=roles/run.builder

#gcloud services enable \
#    run.googleapis.com \
#    cloudbuild.googleapis.com

# faster build, 20x more expensive per minute - takes around 5 minutes
time gcloud builds submit \
  --machine-type=e2_highcpu_32 \
  --tag gcr.io/$PROJECT_ID/cloud-run-cpp

# slower build, 20x cheaper per minute - takes around 25 minutes
#time gcloud builds submit \
#  --machine-type=e2_medium \
#  --tag gcr.io/$PROJECT_ID/cloud-run-cpp

gcloud run deploy $SERVICE_NAME \
    --image=gcr.io/$PROJECT_ID/cloud-run-cpp \
    --region $REGION \
    --allow-unauthenticated

