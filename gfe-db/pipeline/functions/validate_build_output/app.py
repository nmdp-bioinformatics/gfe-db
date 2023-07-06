import os
if __name__ != "app":
    import sys

    # for dev, local path to gfe-db modules
    # ./gfe-db/pipeline/lambda_layers/gfe_db_models (use absolute path)
    sys.path.append(os.environ["GFEDBMODELS_PATH"])
import logging
import re
import json
import boto3
import polars as pl
from gfedbmodels.constants import (
    session,
    infra,
    pipeline,
    database
)

logger = logging.getLogger()
logger.setLevel(logging.INFO)

s3 = session.client("s3", region_name=os.environ["AWS_REGION"])

# SSM Parameters
data_bucket_name = infra.params.DataBucketName

# TODO use Powertools event parser
# @event_parser(model=Order)
def lambda_handler(event, context):
    logger.info(json.dumps(event))

    # assert that Contents is in event.s3_response
    if "Contents" not in event["s3_response"]:
        raise Exception(f"No files found in S3 {data_bucket_name}")

    # Extract file names from event.s3_response.Contents[].Key using jmespath
    output_file_paths = [file["Key"] for file in event["s3_response"]["Contents"]]

    for file_path in output_file_paths:
        logger.info(f"Validating {file_path}")

        # get schema name from file name
        schema = file_path.split("/")[-1].split(".")[0]

        # validate_file_name: Validate file name using regex for the pattern "^[a-z_]+\.csv$"
        pattern = f"^{schema}.{event['version']}.csv$"
        if not re.match(pattern, file_path.split("/")[-1]):
            raise Exception(f"File name {file_path} does not match pattern {pattern}")

        # validate_schema: validate CSV headers by comparing to schema map
        df = load_csv_from_s3(file_path)
        if not set(df.columns) == set(output_headers[schema]):
            raise Exception(f"Columns in {file_path} do not match schema for {schema}")

        # validate_rows: validate CSV data by confirming that rows exist
        if len(df) == 0:
            raise Exception(f"No rows in {file_path}")

        logger.info(
            json.dumps(
                {
                    "file_path": file_path,
                    "schema": schema,
                    "num_cols": len(df.columns),
                    "num_rows": len(df),
                }
            )
        )

    # return payload to step functions
    del event["s3_response"]
    event["file_paths"] = [
        f"s3://{data_bucket_name}/{path}" for path in output_file_paths
    ]
    return event


# Schema map
output_headers = {
    "all_cds": ["gfe_name", "bp_seq_id", "bp_sequence", "aa_seq_id", "aa_sequence"],
    "all_features": [
        "accession",
        "hash_code",
        "locus",
        "rank",
        "sequence",
        "term",
        "gfe_name",
        "allele_id",
        "hla_name",
        "imgt_release",
    ],
    "all_groups": [
        "gfe_name",
        "allele_id",
        "hla_name",
        "ard_id",
        "ard_name",
        "locus",
        "imgt_release",
    ],
    "gfe_sequences": [
        "gfe_name",
        "allele_id",
        "locus",
        "hla_name",
        "seq_id",
        "sequence",
        "length",
        "imgt_release",
    ],
}


# TODO Use S3 streaming https://docs.powertools.aws.dev/lambda/python/latest/utilities/streaming/#getting-started
def load_csv_from_s3(path: str) -> dict:
    data = (
        s3.get_object(Bucket=data_bucket_name, Key=path)["Body"].read().decode("utf-8")
    )

    # save to /tmp
    output_path = f"/tmp/{path.split('/')[-1]}.csv"
    with open(output_path, "w") as file:
        file.write(data)

    # use string for dtype of bp_seq_id
    df = pl.read_csv(output_path, infer_schema_length=0)
    return df


if __name__ == "__main__":
    from pathlib import Path

    path = Path(__file__).parent / "event.json"

    with open(path, "r") as file:
        event = json.load(file)

    lambda_handler(event, "")
