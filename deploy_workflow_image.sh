version=0.02
docker tag xnat-env-bootstrap:latest mmilch01/xnat-ai-workflow:latest
docker login docker.io
docker push mmilch01/xnat-ai-workflow:latest
