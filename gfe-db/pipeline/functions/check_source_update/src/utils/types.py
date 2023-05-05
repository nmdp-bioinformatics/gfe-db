import re
from datetime import datetime
from typing import Optional
from pydantic import BaseModel, validator

# EventBridge Rules trigger UpdateStatus Lambda
# SKIPPED: state machine execution started and skipped 
# PENDING: state machine execution started
# IN_PROGRESS: batch build job triggered 
# SUCCESS: state machine execution succeeded
# FAILED: state machine execution failed
valid_statuses = ["SKIPPED", "PENDING", "IN_PROGRESS", "SUCCESS", "FAILED", None]

def to_datetime(v, fmt="%Y-%m-%dT%H:%M:%SZ"):
    return datetime.strptime(v, fmt)

# validate that date field is ISO 8601 format with timezone
def date_is_iso_8601_with_timezone(v):
    if not re.match(r'^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$', v):
        raise ValueError("Date must be in ISO 8601 format with timezone")
    return v

# validate that url is a valid URL
def url_is_valid(v):
    if not re.match(r'^https?://', v):
        raise ValueError("Url must be a valid URL")
    return v

### Source Config Models ###
class Commit(BaseModel):
    sha: str
    date_utc: str
    message: Optional[str]
    html_url: str

    # validate that commit_sha is a 40 character hex string
    @validator('sha')
    def commit_sha_is_hex(cls, v):
        if not re.match(r'^[0-9a-f]{40}$', v):
            raise ValueError("Commit sha must be a 40 character hex string")
        return v

    # validate that date is ISO 8601 format with timezone
    @validator('date_utc')
    def date_utc_is_iso_8601_with_timezone(cls, v):
        return date_is_iso_8601_with_timezone(v)
    
    # TODO: validate that html_url is a valid URL for a commit

class InputParameters(BaseModel):
    align: bool
    kir: bool
    mem_profile: bool
    limit: Optional[int]

class ExcludedCommitShas(BaseModel):
    description: Optional[str]
    values: list[str]

    # TODO validate that values are hex strings

class TrackedAssetsConfig(BaseModel):
    description: Optional[str]
    values: list[str]

class TargetMetadataConfigItem(BaseModel):
    description: Optional[str]
    asset_path: str
    metadata_regex: str

class TargetMetadataConfig(BaseModel):
    description: Optional[str]
    items: list[TargetMetadataConfigItem]

class RepositoryConfig(BaseModel):
    owner: str
    name: str
    description: Optional[str]
    url: str
    tracked_assets: Optional[TrackedAssetsConfig]  
    target_metadata_config: Optional[TargetMetadataConfig]
    excluded_commit_shas: Optional[ExcludedCommitShas]
    default_input_parameters: InputParameters

    # validate that the url is a valid URL
    @validator('url')
    def url_is_valid(cls, v):
        return url_is_valid(v)

class SourceConfig(BaseModel):
    created_utc: Optional[str]
    updated_utc: Optional[str]
    repositories: dict[str, RepositoryConfig]

    # validate dates are ISO 8601 format with timezone for created_utc, updated_utc
    @validator('created_utc', 'updated_utc')
    def date_utc_is_iso_8601_with_timezone(cls, v):
        return date_is_iso_8601_with_timezone(v)

class ExecutionStateItem(BaseModel):
    repository: RepositoryConfig
    version: int
    execution_date_utc: Optional[str]
    commit: Commit
    input_parameters: Optional[InputParameters]
    status: Optional[str]

    @validator('status')
    def status_is_valid(cls, v):
        if v not in valid_statuses:
            raise ValueError(f"Status must be one of {valid_statuses}")
        return v

    # validate that version is a 4 digit number, position 0 is a number between 1 and 9, and position 1:2 is a number between 0 and 99 and position 3 is 0
    @validator('version')
    def version_is_valid(cls, v):
        if not re.match(r'^[1-9][0-9]{2}0$', str(v)):
            raise ValueError("Version must match '^[1-9][0-9]{2}0$'")
        return v
    
class ExecutionState(BaseModel):
    created_utc: str
    updated_utc: str
    items: list[ExecutionStateItem]

    # validate that items is sorted by commit.date_utc descending
    @validator('items')
    def execution_state_is_sorted(cls, v):
        if not all(v[i].commit.date_utc >= v[i+1].commit.date_utc for i in range(len(v)-1)):
            raise ValueError("Execution history must be sorted by commit.date_utc descending")
        return v
    
    # Releases are formatted as a 4 digit integer incrementing by 10 with a lower bound of 3170
    # Based on the formatting described, validate that no releases are missing from items
    # Remember that the items is sorted by commit.date_utc descending, so release versions will decrement by 10
    @validator('items')
    def execution_state_has_no_missing_releases(cls, v):

        unique_release_versions = sorted(list(set([item.version for item in v])), reverse=True)

        first_version = 3170 
        expected_version = v[0].version
        for version in unique_release_versions:
            if version != expected_version:
                raise ValueError(f"Execution history is missing version {expected_version}")
            expected_version -= 10

            if version == first_version:
                break

        return v 