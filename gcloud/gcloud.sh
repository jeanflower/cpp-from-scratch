#!/bin/bash

PROJECT_ID="cpp-project-449415"
PROJECT_NUMBER=236401783700

gcloud auth login
gcloud auth application-default login

#gcloud projects get-iam-policy $PROJECT_ID --flatten="bindings[].members" --format="table(bindings.role)"

gcloud config set project $PROJECT_ID

gcloud projects add-iam-policy-binding \
    $PROJECT_ID \
    --member=user:$(gcloud config get-value account) \
    --role=roles/serviceusage.serviceUsageConsumer
gcloud auth application-default set-quota-project $PROJECT_ID

gcloud projects add-iam-policy-binding $PROJECT_ID \
    --member=serviceAccount:$PROJECT_NUMBER-compute@developer.gserviceaccount.com \
    --role=roles/run.builder

gcloud services enable \
    run.googleapis.com \
    cloudbuild.googleapis.com

gcloud builds submit --machine-type=e2_highcpu_32 --tag gcr.io/$PROJECT_ID/cloud-run-hello-world

gcloud run deploy --image=gcr.io/$PROJECT_ID/cloud-run-hello-world
