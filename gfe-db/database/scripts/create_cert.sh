#!/bin/bash

# Part of user data, to be run on the database instance on initialization

echo "Provisioning SSL certificate..."
export NEO4J_HOME=/opt/bitnami/neo4j
export HOST_DOMAIN=gfe-db.cloudftl.com
export ADMIN_EMAIL=gclindsey@gmail.com

certbot certonly -n \
  -d $HOST_DOMAIN \
  --standalone \
  -m $ADMIN_EMAIL \
  --agree-tos \
  --redirect

chgrp -R neo4j /etc/letsencrypt/*
chmod -R g+rx /etc/letsencrypt/*
mkdir -p $NEO4J_HOME/certificates/{bolt,cluster,https}/trusted

for certsource in bolt cluster https; do
  sudo ln -sf "/etc/letsencrypt/live/$HOST_DOMAIN/fullchain.pem" "$NEO4J_HOME/certificates/$certsource/neo4j.cert"
  sudo ln -sf "/etc/letsencrypt/live/$HOST_DOMAIN/privkey.pem" "$NEO4J_HOME/certificates/$certsource/neo4j.key"
  sudo ln -sf "/etc/letsencrypt/live/$HOST_DOMAIN/fullchain.pem" "$NEO4J_HOME/certificates/$certsource/trusted/neo4j.cert"
done

sudo chgrp -R neo4j $NEO4J_HOME/certificates/*
sudo chmod -R g+rx $NEO4J_HOME/certificates/*

exit 0
