version=0.01
docker tag xnat_env_bootstrap:latest mmilch01/xnat-ai-workflow:$version
docker login docker.io
docker push mmilch01/xnat-ai-workflow:$version
