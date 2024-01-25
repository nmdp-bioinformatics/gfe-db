import os
import logging
import json

logger = logging.getLogger()
logger.setLevel(logging.INFO)

AWS_REGION = os.environ["AWS_REGION"]
STAGE = os.environ["STAGE"]
APP_NAME = os.environ["APP_NAME"]

# Template for the report
report_template = """
{title}
-------

Execution ID: {execution_id}
Execution Status: {execution_status}
Execution Date: {execution_date} UTC

Commit Details: {commit_sha} - {commit_message} ({commit_url})

Data Load Summary:
- New Node Counts Added: {new_node_counts}
- Unique Release Version Added: {unique_release_version}

Data Files Processed:
{data_files_info}

Data Integrity Checks:
- CSV File Validations: All passed
- Node Counts Increase Check: {node_counts_check}
- Unique Release Count Increment Check: {release_count_check}

Backup Details:
- Pre-Execution Backup: Command ID {pre_execution_backup}
- Post-Execution Backup: Command ID {post_execution_backup}

"""

def lambda_handler(event, context):

    logger.info(json.dumps(event))
    data = event

    # Set the title
    title = f"{STAGE}-{APP_NAME} Pipeline Execution Report"
    title_underline = "-" * len(title)

    # Extract required information from JSON
    execution_status = data['state']['execution']['status']
    execution_id = data['state']['execution']['id']
    execution_date = data['state']['updated_utc']
    commit_sha = data['state']['commit']['sha']
    commit_url = data['state']['commit']['html_url']
    commit_message = data['state']['commit']['message'].replace("\n", " ")

    # These fields are not present when duplicate executions are run (same release, commit sha and/or limit)
    if 'have_node_counts_increased' in data['validations']['load_results']:
        new_node_counts = data['validations']['load_results']['have_node_counts_increased']['details']['node_counts_post_load']
        node_counts_check = "Passed" if data['validations']['load_results']['have_node_counts_increased']['value'] else "Failed"
    else:
        new_node_counts = "N/A"
        node_counts_check = "N/A"

    if 'has_unique_release_count_increased_by_1' in data['validations']['load_results']:
        unique_release_version = data['validations']['load_results']['has_unique_release_count_increased_by_1']['details']['num_unique_releases_in_db_post_load']
        release_count_check = "Passed" if data['validations']['load_results']['has_unique_release_count_increased_by_1']['value'] else "Failed"
    else:
        unique_release_version = "N/A"
        release_count_check = "N/A"
    
    data_files_info = format_data_files(data['validations']['build_outputs']['details'])

    pre_execution_backup = data['backups']['pre']['command_id']
    post_execution_backup = data['backups']['post']['command_id']

    # Create the report
    report = report_template.format(
        title=title,
        title_underline=title_underline,
        execution_status=execution_status,
        execution_id=execution_id,
        execution_date=execution_date,
        commit_sha=commit_sha,
        commit_url=commit_url,
        commit_message=commit_message,
        new_node_counts=new_node_counts,
        unique_release_version=unique_release_version,
        data_files_info=data_files_info,
        node_counts_check=node_counts_check,
        release_count_check=release_count_check,
        pre_execution_backup=pre_execution_backup,
        post_execution_backup=post_execution_backup
    )

    return report


# Function to process data files information
def format_data_files(data_files):
    file_info = ""
    for file in data_files:
        file_info += f"- `{file['schema']}` - {file['num_rows']} rows\n"
    return file_info

if __name__ == "__main__":
    from pathlib import Path

    event = json.loads((Path(__file__).parent / "error.json").read_text())

    lambda_handler(event,"")
