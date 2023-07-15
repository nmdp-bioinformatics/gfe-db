import os
from gfedbmodels.constants import (
    session,
    infra,
    pipeline
)

# Environment
PIPELINE_SOURCE_CONFIG_S3_PATH = os.environ["PIPELINE_SOURCE_CONFIG_S3_PATH"]

(
    GITHUB_REPOSITORY_OWNER, 
    GITHUB_REPOSITORY_NAME, 
    execution_state_table_name, 
    gfedb_processing_queue_url, 
    execution_state_table_fields
) = (
    pipeline.params.GitHubSourceRepository["owner"],
    pipeline.params.GitHubSourceRepository["name"],
    pipeline.params.ExecutionStateTableName,
    pipeline.params.GfeDbProcessingQueueUrl,
    pipeline.params.ExecutionStateTableFields
)

data_bucket_name = infra.params.DataBucketName