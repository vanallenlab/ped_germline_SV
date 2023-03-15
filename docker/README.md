## Docker images

This subdirectory contains the Dockerfiles to build the Docker images specific to this study.  

Example command for building the `PedSV` image, which should be executed from the directory that contains a local copy of this repo:  

```
export my_tag=my_test_image_hash
docker build -f ped_germline_SV/docker/PedSV/Dockerfile --progress plain --tag vanallenlab/pedsv:$my_tag `pwd`
```

After building, you will (usually) want to push the image to Dockerhub and also update the `latest` tag:  

```
docker push vanallenlab/pedsv:$my_tag
docker tag vanallenlab/pedsv:$my_tag vanallenlab/pedsv:latest
docker push vanallenlab/pedsv:latest
```
