import re
from pydantic import BaseModel, validator

### Source Config Models ###
class InputParameterConfig(BaseModel):
    align: bool
    kir: bool
    mem_profile: bool
    limit: str

    # validate that limit is an integer
    @validator('limit')
    def limit_is_integer(cls, v):
        if not v.isnumeric():
            raise ValueError('limit must be an integer')
        return v

class ExecutionHistoryConfig(BaseModel):
    version: int
    commit_sha: str
    date: str
    status: str
    input_parameters: InputParameterConfig

    # validate that version is a 4 digit number, position 0 is a number between 1 and 9, and position 1:2 is a number between 0 and 99 and position 3 is 0
    @validator('version')
    def version_is_valid(cls, v):
        if not re.match(r'^[1-9][0-9]{2}0$', str(v)):
            raise ValueError("version must match '^[1-9][0-9]{2}0$'")
        return v

    # validate that commit_sha is a 40 character hex string
    @validator('commit_sha')
    def commit_sha_is_hex(cls, v):
        if not re.match(r'^[0-9a-f]{40}$', v):
            raise ValueError('commit_sha must be a 40 character hex string')
        return v

    # validate that date field is ISO 8601 format with timezone
    @validator('date')
    def date_is_iso_8601_with_timezone(cls, v):
        if not re.match(r'^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$', v):
            raise ValueError('date must be in ISO 8601 format with timezone')
        return v
    
    # validate the status is on of 'SUCCESS', 'FAILURE', or 'IN_PROGRESS'
    @validator('status')
    def status_is_valid(cls, v):
        if v not in ['SUCCESS', 'FAILURE', 'IN_PROGRESS']:
            raise ValueError('status must be one of "SUCCESS", "FAILURE", or "IN_PROGRESS"')
        return v

class RepositoryConfig(BaseModel):
    owner: str
    name: str
    url: str
    tracked_assets: list[str]
    default_input_parameters: InputParameterConfig
    execution_history: list[ExecutionHistoryConfig]

    # validate that the url is a valid URL
    @validator('url')
    def url_is_valid(cls, v):
        if not re.match(r'^https?://', v):
            raise ValueError('url must be a valid URL')
        return v
    
    # validate that execution_history is sorted by date descending
    @validator('execution_history')
    def execution_history_is_sorted_by_date_descending(cls, v):
        for i in range(0, len(v) - 1):
            if v[i].date < v[i+1].date:
                raise ValueError('execution_history must be sorted by date descending')
        return v

class SourceConfig(BaseModel):
    repositories: dict[str, RepositoryConfig]
