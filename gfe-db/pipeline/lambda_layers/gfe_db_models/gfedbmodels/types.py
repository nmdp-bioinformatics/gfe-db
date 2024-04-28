import re
from datetime import datetime
from typing import Optional
from enum import Enum
from pydantic import BaseModel, validator, root_validator
import jmespath
from .utils import restore_nested_json, filter_nested_nulls


class ExecutionStatus(str, Enum):
    """
    ExecutionStatus is synced using the Step Functions DynamoDB integration:
    NOT_PROCESSED: never processed (set by CheckSourceUpdate) ✅
    SKIPPED: never processed (set by CheckSourceUpdate) ✅
    PENDING: state machine execution started (set by CheckSourceUpdate) ✅
    BUILD_IN_PROGRESS: build started (set by State Machine) ✅
    BUILD_SUCCESS: build succeeded (set by State Machine) ✅
    LOAD_IN_PROGRESS: load started (set by State Machine) ✅
    LOAD_SUCCESS: load succeeded (set by State Machine) ✅
    LOAD_FAILED: load failed (set by State Machine) ✅
    LOAD_INVALID: load invalid from query results (set by State Machine) ✅
    LOAD_SKIPPED: load skipped (set by State Machine) ✅
    BUILD_FAILED: build failed (set by State Machine) ✅
    EXECUTION_FAILED: build or load failed (set by State Machine) ✅
    ABORTED: build or load aborted (set by UpdateExecutionState) ✅
    """

    NOT_PROCESSED = "NOT_PROCESSED"
    SKIPPED = "SKIPPED"
    PENDING = "PENDING"
    BUILD_IN_PROGRESS = "BUILD_IN_PROGRESS"
    BUILD_SUCCESS = "BUILD_SUCCESS"
    BUILD_FAILED = "BUILD_FAILED"
    LOAD_IN_PROGRESS = "LOAD_IN_PROGRESS"
    LOAD_COMPLETE = "LOAD_COMPLETE"
    LOAD_SUCCESS = "LOAD_SUCCESS"
    LOAD_FAILED = "LOAD_FAILED"
    LOAD_INVALID = "LOAD_INVALID"
    LOAD_SKIPPED = "LOAD_SKIPPED"
    EXECUTION_FAILED = "EXECUTION_FAILED"
    ABORTED = "ABORTED"

    @classmethod
    def __contains__(cls, item):
        return item in cls.__members__


def str_to_datetime(v, fmt="%Y-%m-%dT%H:%M:%S.%fZ"):
    return datetime.strptime(v, fmt)


def str_from_datetime(v, fmt="%Y-%m-%dT%H:%M:%SZ"):
    return v.strftime(fmt)


# validate that date field is ISO 8601 format with timezone
def date_is_iso_8601_with_timezone(v):
    # Check if the date is already in the desired ISO 8601 format with 3 milliseconds
    if re.match(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d{3}Z$", v):
        return v

    # Check if the date is in ISO 8601 format with fractional seconds (arbitrary number of digits)
    match = re.match(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d+Z$", v)
    if match:
        fractional_seconds = v.split(".")[1].split("Z")[0]
        # Truncate or pad fractional seconds to 3 digits
        truncated_fractional_seconds = fractional_seconds[:3].ljust(3, "0")
        return v.replace(fractional_seconds, truncated_fractional_seconds)

    # Check if the date is in ISO 8601 format without fractional seconds
    if re.match(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$", v):
        # Add milliseconds and return
        return v[:-1] + ".000Z"

    raise ValueError("Date must be in ISO 8601 format with timezone")


# validate that url is a valid URL
def url_is_valid(v):
    if not re.match(r"^https?://", v):
        raise ValueError("Url must be a valid URL")
    return v


release_version_re = r"^[1-9][0-9]{1,2}0$"

def version_is_valid(v, return_bool=False):

    if not return_bool:
        if not re.match(release_version_re, str(v)):
            raise ValueError(
                f"Release version must match regex pattern `{release_version_re}`"
            )
        return v
    
    elif return_bool:
        return bool(re.match(release_version_re, str(v)))


# validate that commit sha is a 40 character hex string
def commit_sha_is_hex(v):
    if not re.match(r"^[0-9a-f]{40}$", v):
        raise ValueError(f"Commit sha must be a 40 character hex string matching regex pattern `^[0-9a-f]{40}$`: {v}")
    return v


def s3_path_is_valid(v):
    if not re.match(r"^s3://", v):
        raise ValueError("S3 path must be a valid S3 path (s3://<bucket>/<key>)")
    return v


### Source Config Models ###
class Commit(BaseModel):
    sha: str
    date_utc: str
    message: Optional[str] = None
    html_url: str
    # # TODO the url field returns a JSON object with a list of files; see all-branches for more info
    # url: str

    @validator("sha")
    def _commit_sha_is_hex(cls, v):
        return commit_sha_is_hex(v)

    # implement jmespath mapping to create the commit Class
    @classmethod
    def from_response_json(cls, response_json):
        return cls(
            sha=jmespath.search("sha", response_json),
            date_utc=jmespath.search("commit.committer.date", response_json),
            message=jmespath.search("commit.message", response_json),
            html_url=jmespath.search("html_url", response_json),
        )

    # validate that date is ISO 8601 format with timezone
    @validator("date_utc")
    def date_utc_is_iso_8601_with_timezone(cls, v):
        return date_is_iso_8601_with_timezone(v)



class InputParameters(BaseModel):
    align: bool
    kir: bool
    mem_profile: bool
    limit: Optional[int] = -1
    use_existing_build: Optional[bool] = False
    skip_load: Optional[bool] = False

    # validate that limit is an integer equalt to -1 or greater than 0 but not equal to 0
    @validator("limit")
    def limit_is_valid(cls, v):
        if v == -1:
            return v
        elif v > 0:
            return v
        else:
            raise ValueError("Limit must be an integer equal to -1 or greater than 0")


class ExcludedCommitShas(BaseModel):
    description: Optional[str] = None
    values: list[str]

    @validator("values")
    def _commit_shas_are_hex(cls, v):
        for sha in v:
            sha = commit_sha_is_hex(sha)
        return v


class TrackedAssetsConfig(BaseModel):
    description: Optional[str] = None
    values: list[str]


class TargetMetadataConfigItem(BaseModel):
    description: Optional[str] = None
    asset_path: str  # Path (in remote git repository, on GitHub) to the file to check for strings matching a release version
    metadata_regex: str  # Contextual metadata for the commit. This is the regex to match the release version string that might be found in the assets being checked


class TargetMetadataConfig(BaseModel):
    description: Optional[str] = None
    items: list[TargetMetadataConfigItem]


class RepositoryConfig(BaseModel):
    owner: str
    name: str
    description: Optional[str] = None
    url: str
    tracked_assets: Optional[TrackedAssetsConfig] = None
    target_metadata_config: Optional[TargetMetadataConfig] = None
    excluded_commit_shas: Optional[ExcludedCommitShas] = None
    default_input_parameters: Optional[InputParameters] = None

    # validate that the url is a valid URL
    @validator("url")
    def url_is_valid(cls, v):
        return url_is_valid(v)


class ExecutionDetailsConfig(BaseModel):
    id: str = (
        None  # Refers to execution id in Step Functions and is only set if the commit is processed
    )
    invocation_id: str = (
        None  # One invocation can have multiple executions depending on how many release versions are given
    )
    version: int
    status: str
    date_utc: Optional[str] = None
    input_parameters: Optional[InputParameters] = None
    s3_path: Optional[str] = None
    receipt_handle: Optional[str] = None
    # TODO move execution errors here

    @validator("status")
    def status_is_valid(cls, v):
        if v not in ExecutionStatus.__members__:
            raise ValueError(
                f"Status must be one of {[value.value for value in ExecutionStatus.__members__.values()]}"
            )
        return v

    # validate that version is a 4 digit number, position 0 is a number between 1 and 9, and position 1:2 is a number between 0 and 99 and position 3 is 0
    @validator("version")
    def _version_is_valid(cls, v):
        return version_is_valid(v)


class SourceConfig(BaseModel):
    created_utc: Optional[str] = None
    updated_utc: Optional[str] = None
    repositories: dict[str, RepositoryConfig]

    # validate dates are ISO 8601 format with timezone for created_utc, updated_utc
    @validator("created_utc", "updated_utc")
    def date_utc_is_iso_8601_with_timezone(cls, v):
        return date_is_iso_8601_with_timezone(v)


# Uses similar schema from Step Functions Fail state
class ExecutionError(BaseModel):
    message: str
    cause: str


# One item in the ExecutionState table
# The Primary Key is commit.sha, or commit__sha in the table
class ExecutionStateItem(BaseModel):
    created_utc: Optional[str] = (
        None  # Partial updates may not be able to include timestamps
    )
    updated_utc: Optional[str] = (
        None  # Partial updates may not be able to include timestamps
    )
    repository: Optional[RepositoryConfig]
    commit: Commit
    execution: ExecutionDetailsConfig
    error: Optional[ExecutionError] = None
    s3_path: Optional[str] = None

    @classmethod
    def from_execution_state_item_json(cls, execution_state_item: dict):
        # Items from table are separated by "__" because "." is not allowed in DynamoDB
        execution_state_item = restore_nested_json(execution_state_item, split_on="__")
        return cls(**execution_state_item)

    # validate s3 path uses s3://<bucket>/<key> format
    @validator("s3_path")
    def s3_path_is_valid(cls, v):
        return s3_path_is_valid(v)


class ExecutionState(BaseModel):
    created_utc: str
    items: list[ExecutionStateItem]

    @root_validator(pre=True)
    def set_items_created_utc(cls, values):
        timestamp_utc = values.get("created_utc")
        items = values.get("items", [])
        try:
            for item in items:
                item.created_utc = timestamp_utc
        except:
            for item in items:
                item["created_utc"] = timestamp_utc
        return values

    # validate that items is sorted by commit.date_utc descending
    @validator("items")
    def execution_state_is_sorted(cls, v):
        if not all(
            v[i].commit.date_utc >= v[i + 1].commit.date_utc for i in range(len(v) - 1)
        ):
            raise ValueError(
                "Execution history must be sorted by commit.date_utc descending"
            )
        return v

    # Releases are formatted as a 3 or 4 digit integer incrementing by 10 with a lower bound of 300
    # Based on the formatting described, validate that no releases are missing from items
    # Remember that the items is sorted by commit.date_utc descending and that release versions decrement by 10
    @validator("items")
    def execution_state_has_no_missing_releases(cls, items):
        actual_unique_release_versions = sorted(
            list(set([item.execution.version for item in items])), reverse=True
        )

        first_version = 300
        expected_version = items[0].execution.version
        for idx, actual_version in enumerate(actual_unique_release_versions):

            # Make custom adjust for releases under 3100 to account for inconsistent versioning
            if actual_version == 390:
                expected_version = int(
                    str(expected_version)[:1] + str(expected_version)[2:]
                )

            # If the version is not the expected version, raise an error
            if actual_version != expected_version:
                raise ValueError(
                    f"Execution history is missing version {expected_version}"
                )

            # If the version is the first version, then the expected version should be the last version
            if actual_version == first_version:
                if idx != len(actual_unique_release_versions) - 1:
                    raise ValueError(
                        f"Execution history has an unexpected version {actual_version}"
                    )
                break

            # Since release versioning increments by 10 in a complete dataset (all version), get the expected value using math
            expected_version -= 10

        return items


class ExecutionPayloadItem(BaseModel):
    id: str
    version: int
    invocation_id: str
    commit_sha: str
    input_parameters: InputParameters
    s3_path: str

    @validator("version")
    def _version_is_valid(cls, v):
        return version_is_valid(v)

    @validator("commit_sha")
    def _commit_sha_is_hex(cls, v):
        return commit_sha_is_hex(v)

    # create the ExecutionPayloadItem from ExecutionStateItem
    @classmethod
    def from_execution_state_item(cls, execution_state_item):
        return cls(
            id=execution_state_item.execution.id,
            version=execution_state_item.execution.version,
            invocation_id=execution_state_item.execution.invocation_id,
            commit_sha=execution_state_item.commit.sha,
            input_parameters=execution_state_item.execution.input_parameters,
            s3_path=execution_state_item.execution.s3_path,
        )
