#!/usr/bin/env sh

for release in $(seq 3300 10 3380);
do
  echo $release
  docker build -t gfe-db:$release --build-arg IMGT="${release}" --build-arg AN=True .
done
