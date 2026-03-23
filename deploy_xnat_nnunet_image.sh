version=0.03
docker tag xnat-env-nnunet:latest mmilch01/xnat-env-nnunet:latest
docker login docker.io
docker push mmilch01/xnat-env-nnunet:latest
