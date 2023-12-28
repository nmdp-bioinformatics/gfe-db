#!/bin/bash -x

# Part of user data, to be run on the database instance on initialization, or later for renewal

echo "Renewing SSL certificate..."
export NEO4J_HOME=/var/lib/neo4j

# Passed from command line
DOMAIN=$1
# ADMIN_EMAIL=$2

# certbot certonly -n \
#   -d $DOMAIN \
#   --standalone \
#   -m $ADMIN_EMAIL \
#   --agree-tos \
#   --redirect

certbot renew

chgrp -R neo4j /etc/letsencrypt/*
chmod -R g+rx /etc/letsencrypt/*
mkdir -p $NEO4J_HOME/certificates/{bolt,cluster,https}/trusted

for certsource in bolt cluster https; do
  ln -sf "/etc/letsencrypt/live/$DOMAIN/fullchain.pem" "$NEO4J_HOME/certificates/$certsource/neo4j.cert"
  ln -sf "/etc/letsencrypt/live/$DOMAIN/privkey.pem" "$NEO4J_HOME/certificates/$certsource/neo4j.key"
  ln -sf "/etc/letsencrypt/live/$DOMAIN/fullchain.pem" "$NEO4J_HOME/certificates/$certsource/trusted/neo4j.cert"
done

chgrp -R neo4j $NEO4J_HOME/certificates/*
chmod -R g+rx $NEO4J_HOME/certificates/*

exit 0
