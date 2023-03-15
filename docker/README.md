## Docker images

This subdirectory contains the Dockerfiles to build the Docker images specific to this study.  

Example command for building the `PedSV` image:  

```
$ ./build_docker.sh my_tag
```

Where `my_tag` is the tag to be applied to the rebuilt Docker image, and should ideally be something specific like the Github commit hash.  
