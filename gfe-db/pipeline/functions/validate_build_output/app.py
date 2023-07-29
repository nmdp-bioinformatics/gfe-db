import os
import logging
from datetime import datetime
from dateutil import tz
import re
import json
import s3fs
import boto3
import polars as pl

logger = logging.getLogger()
logger.setLevel(logging.INFO)

STAGE = os.environ["STAGE"]
APP_NAME = os.environ["APP_NAME"]
AWS_REGION = os.environ["AWS_REGION"]

fs = s3fs.S3FileSystem()

session = boto3.Session(region_name=AWS_REGION)
ssm = session.client("ssm")
s3 = session.client("s3")

# get data bucket name
data_bucket_name = ssm.get_parameter(
    Name=f"/{APP_NAME}/{STAGE}/{AWS_REGION}/DataBucketName"
)["Parameter"]["Value"]

def lambda_handler(event, context):

    logger.info(json.dumps(event))
    execution_start_time = datetime.strptime(event['execution_context']['Execution']['StartTime'], '%Y-%m-%dT%H:%M:%S.%fZ').replace(tzinfo=tz.tzutc())

    # TODO get the expected input from execution context and validate against this
    # TODO Remove this and use the output of validation against the expected input
    # Remove failed executions
    parsed_inputs = get_successful_inputs(event['build_output']['input'])
    if not parsed_inputs:
        raise Exception("No inputs found")

    releases = [ item['RELEASES'] for item in parsed_inputs ]
    errors = []
    reports = []
    for release in releases:
        
        # build path to CSV directory
        csv_dir = f"data/{release}/csv"

        # get list of CSV files in directory
        csv_file_objs = list_s3_objects(data_bucket_name, csv_dir)

        # validate all files were created after the execution start time
        num_errors = 0
        for key, obj in csv_file_objs.items():

            obj['data'] = pl.read_csv(f"s3://{data_bucket_name}/{csv_dir}/{key}", infer_schema_length=0)
            obj["details"] = {}

            # Validate the timestamp is after the execution start time
            obj["details"]["is_valid_timestamp"] = obj['created_utc'] > execution_start_time
            if not obj["details"]["is_valid_timestamp"]:
                error_msg = f"CSV file timestamp ({str(obj['created_utc'])}) preceeds execution start time ({execution_start_time}): {key}"
                logger.error(error_msg)
                errors.append(error_msg)
                num_errors += 1

            # Validate the file name is correct
            obj["details"]["is_valid_csv_filename"] = is_valid_csv_filename(key, release)
            if not obj["details"]["is_valid_csv_filename"]:
                error_msg = f"CSV file name is incorrect: {key}"
                logger.error(error_msg)
                errors.append(error_msg)
                num_errors += 1

            # # Validate the headers are correct
            obj["details"]["is_valid_csv_headers"] = set(obj['data'].columns) == set(csv_headers[key.split('.')[0]])
            if not obj["details"]["is_valid_csv_headers"]:
                error_msg = f"CSV file headers are incorrect: {key}"
                logger.error(error_msg)
                errors.append(error_msg)
                num_errors += 1

            # # Validate rows are present
            obj["details"]["is_valid_csv_rows"] = obj['data'].shape[0] > 0
            if not obj["details"]["is_valid_csv_rows"]:
                error_msg = f"CSV file has no rows: {key}"
                logger.error(error_msg)
                errors.append(error_msg)
                num_errors += 1

            # TODO get the expected total number of rows from hla.dat and validate the row count against that
            # # Validate expected number of rows

            # TODO Validate the data types of each column
            # # Validate data types

            report = {
                "num_errors": num_errors,
                "file_path": f"s3://{data_bucket_name}/{csv_dir}/{key}",
                "schema": key.split('.')[0],
                "cols": obj['data'].columns,
                "num_rows": len(obj['data']),
                **{k: v.strftime('%Y-%m-%dT%H:%M:%S.%f')[:-3]+"Z" if isinstance(v, datetime) else v for k, v in obj.items() if k != 'data'}
            }
            reports.append(report)
    
    logger.info(json.dumps(reports))
    if errors:
        logger.error(json.dumps(errors))
        raise Exception(f'Validation failed:\n{json.dumps(errors, indent=2)}')

    return {
        "input": parsed_inputs,
        "reports": reports,
    }

# TODO implement Pydantic classes for each CSV schema
# Schema map
csv_headers = {
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
        "acc_name",
        "locus",
        "hla_name",
        "seq_id",
        "sequence",
        "length",
        "imgt_release",
    ],
}


def is_valid_csv_filename(key, release):
    return re.match(re.compile(f'^{key.split("/")[-1].split(".")[0]}.{release}.csv$'), key.split("/")[-1]) is not None


def get_successful_inputs(event_input):
    return [ item for item in event_input if "Error" not in item.keys() ]


def list_s3_objects(bucket_name: str, prefix: str) -> list:
    
    # list objects in bucket
    objs = s3.list_objects_v2(
            Bucket=bucket_name,
            Prefix=prefix,
    )['Contents']
    return { obj['Key'].split('/')[-1]: {"created_utc": obj['LastModified']} for obj in objs }

# Needed to serialize datetime objects in JSON responses
class DatetimeEncoder(json.JSONEncoder):
    """
    Helps convert datetime objects to pure strings in AWS service API responses. Does not
    convert timezone information.

    Extend `json.JSONEncoder`. 
    """

    def default(self, obj):
        try:
            return super().default(obj)
        except TypeError:
            return str(obj)


if __name__ == "__main__":
    from pathlib import Path

    event_path = Path(__file__).parent / "event.json"

    with open(event_path, "r") as file:
        event = json.load(file)

    lambda_handler(event,"")
