#!/bin/bash -x


gfe_bucket=$(aws ssm get-parameters \
    --region $(curl -s http://169.254.169.254/latest/meta-data/placement/availability-zone | sed 's/\(.*\)[a-z]/\1/') \
    --names "/gfe-db/dev/us-east-1/DataBucketName" \
    | jq -r '.Parameters | map(select(.Version == 1))[0].Value')

if [[ -z $gfe_bucket ]]; then
    echo "S3 bucket not found."
    exit 1
else
    echo "Found S3 bucket: $gfe_bucket"
fi

# TODO update start/stop commands for Bitnami AMI
# Stop, perform dump and restart database
echo "Stopping Neo4j..."
sudo systemctl stop neo4j
echo "Backing up graph data..."

# mkdir -p /var/lib/neo4j/backups # This directory is created by the user data script on initial boot to avoid permissions issues
# TODO update Neo4j paths for Bitnami AMI
sudo neo4j-admin dump \
    --database=neo4j \
    --to=/var/lib/neo4j/backups/gfedb.dump

# Confirm that Neo4j has restarted
echo "Restarting database..."
sudo systemctl start neo4j
until $(curl --output /dev/null --silent --head --fail http://localhost:7474); do \
    printf '.' ; \
    sleep 1 ; \
    done
printf "%s\n" " "
echo "Neo4j is ready"

# Copy Neo4j startup logs
echo "Writing logs..."
journalctl -e -u neo4j > /tmp/logs/neo4j-s3-backup.log
echo "Copying logs to S3..."
aws s3 cp /tmp/logs/neo4j-s3-backup.log s3://$gfe_bucket/backups/database/$(date +'%Y/%m/%d/%H')/neo4j-s3-backup.log

# Copy to S3
echo "Copying data to S3..."
aws s3 cp /var/lib/neo4j/backups/gfedb.dump s3://$gfe_bucket/backups/database/$(date +'%Y/%m/%d/%H')/gfedb.dump
echo "Done"

exit 0
