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
    """Validates the build output artifacts against the original execution input object."""

    # TODO can specify execution input fields instead of passing entire object
    logger.info(json.dumps(event))
    execution_start_time = datetime.strptime(event['execution_context']['Execution']['StartTime'], '%Y-%m-%dT%H:%M:%S.%fZ').replace(tzinfo=tz.tzutc())

    # TODO get the expected input from execution context and validate against this
    # TODO Remove this and use the output of validation against the expected input
    # # In other words, the validation should be "blind" to what the build output is, it should only
    # # tells us how it's different from the expected. Therefore the expected output is constructed from $$.Execution.Input

    # expected input is the execution input
    execution_input = event['execution_context']['Execution']['Input']['input']

    # TODO clean up; left over from previous iteration that processed multipe releases per execution
    # releases = list(set([ item['RELEASES'] for item in execution_input ]))
    releases = [execution_input['version']]

    # reports for all release builds
    reports = []

    # errors for all release builds used only for logging, not used for validation logic
    errors = []
    for release in releases:

        release_report = {
            "release": release,
            "details": [],
            "errors": []
        }

        # TODO hard-coded S3 paths should be replaced with SSM parameters at some point
        # build path to CSV directory
        csv_dir = f"data/{release}/csv"
        
        # expected artifacts from state machine input
        release_report["expected_artifacts"] = [ f"{key}.{release}.csv" for key in csv_headers.keys() ]

        # Validate that the S3 prefix exists and has data
        try:
            csv_file_objs = list_s3_objects(data_bucket_name, csv_dir)
        except KeyError:
            error_msg = f"CSV directory does not exist: {csv_dir}"
            logger.error(error_msg)
            release_report["errors"].append(error_msg)
            release_report["is_valid_build"] = False
            reports.append(release_report)
            errors.append(error_msg)
            continue

        # Validate that all expected files are present
        if set(csv_file_objs.keys()) != set(release_report["expected_artifacts"]):
            error_msg = f"CSV files do not match expected artifacts: {csv_dir}"
            logger.error(error_msg)
            release_report["errors"].append(error_msg)
            release_report["is_valid_build"] = False
            reports.append(release_report)
            errors.append(error_msg)
            continue

        for key, obj in csv_file_objs.items():
            obj_errors = []

            obj['data'] = pl.read_csv(f"s3://{data_bucket_name}/{csv_dir}/{key}", infer_schema_length=0)
            obj["details"] = {}

            # Validate the file's timestamp is after the execution start time
            obj["details"]["is_valid_csv_timestamp"] = obj['created_utc'] > execution_start_time
            if not obj["details"]["is_valid_csv_timestamp"]:
                error_msg = f"CSV file timestamp ({str(obj['created_utc'])}) preceeds execution start time ({execution_start_time}): {key}"
                logger.error(error_msg)
                obj_errors.append(error_msg)

            # Validate the file name is correct
            obj["details"]["is_valid_csv_filename"] = is_valid_csv_filename(key, release)
            if not obj["details"]["is_valid_csv_filename"]:
                error_msg = f"CSV file name is incorrect: {key}"
                logger.error(error_msg)
                obj_errors.append(error_msg)

            # Validate the headers are correct
            obj["details"]["is_valid_csv_headers"] = set(obj['data'].columns) == set(csv_headers[key.split('.')[0]])
            if not obj["details"]["is_valid_csv_headers"]:
                error_msg = f"CSV file headers are incorrect: {key}\n\tExpected: {csv_headers[key.split('.')[0]]}\n\tActual: {obj['data'].columns}"
                logger.error(error_msg)
                obj_errors.append(error_msg)

            # Validate rows are present
            obj["details"]["is_valid_csv_rows"] = obj['data'].shape[0] > 0
            if not obj["details"]["is_valid_csv_rows"]:
                error_msg = f"CSV file has no rows: {key}"
                logger.error(error_msg)
                obj_errors.append(error_msg)

            # TODO get the expected total number of rows from hla.dat and validate the row count against that
            # # # Validate expected number of rows and replace is_valid_csv_rows
            # obj["details"]["is_valid_csv_num_rows"] = None

            # TODO Validate the data types of each column
            # # # Validate data types
            # obj["details"]["is_valid_csv_data_types"] = None

            obj_report = {
                "schema": key.split('.')[0],
                "release": release,
                "file_path": f"s3://{data_bucket_name}/{csv_dir}/{key}",
                "cols": obj['data'].columns,
                "num_rows": len(obj['data']),
                **{k: v.strftime('%Y-%m-%dT%H:%M:%S.%f')[:-3]+"Z" if isinstance(v, datetime) else v for k, v in obj.items() if k != 'data'},
                "num_errors": len(obj_errors),
                "is_valid_csv": all(obj["details"].values()),
            }
            if obj_errors:
                obj_report["errors"] = obj_errors
                obj_report["num_errors"] = len(obj_errors)
                errors.extend(obj_errors)

            release_report["details"].append(obj_report)
            if obj_errors:
                release_report["errors"].append(obj_errors)
    
        release_report["is_valid_build"] = all([obj_report["is_valid_csv"] for obj_report in release_report["details"]])
        reports.append(release_report)

    logger.info(json.dumps(reports))
    if errors:
        error_msg = f'Validation failed:\n{json.dumps(errors)}'
        logger.error(error_msg)

    ### Update output ###
    valid_release_builds = [ release['release'] for release in reports if release['is_valid_build'] ]

    # TODO clean up; left over from previous iteration that processed multipe releases per execution
    # payload = list(filter(lambda x: x["RELEASES"] in valid_release_builds, execution_input))

    payload = list(filter(lambda x: x["version"] in valid_release_builds, [execution_input]))

    result = {
        "validated": payload,
        "build_details": reports,
        "has_valid_payload": len(payload) > 0,
    }

    return result

# TODO implement Pydantic classes for managing and validating CSV schemas
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
