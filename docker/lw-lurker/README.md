# LOLAweb lurker

First build the image based on instructions:

```
docker build --no-cache -t lw-lurker .
```

Next make sure that `LWREF` and `LWLOCAL` environment variables are set:

```
export LWREF="your/path/to/reference/data/"
export LWLOCAL="your/path/to/cache/log/temp/data"
```

Run the `lw-db` image:

```
docker pull mongo
docker run -d --name lw-db -p 27017:27017 -v ~/data:/data/db mongo --bind_ip_all
```

Run the frontend with volumes mounted and linked to `lw-db`:

```
docker run -d \
  --name lw-lurker \
  -e LWREF=$LWREF \
  -e LWLOCAL=$LWLOCAL \
  --volume ${LWLOCAL}:${LWLOCAL} \
  --volume ${LWREF}:${LWREF} \
  --volume ${LWLOCAL}/temp:/tmp \
  --link lw-db:lw-lurker \
  --memory="8000m" \
  lw-lurker
```

