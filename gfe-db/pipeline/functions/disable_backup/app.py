import os
import logging
import json
import boto3

logger = logging.getLogger()
logger.setLevel(logging.INFO)

# variables
AWS_REGION = os.environ['AWS_REGION']
STAGE = os.environ['STAGE']
APP_NAME = os.environ['APP_NAME']

session = boto3.Session(region_name=AWS_REGION)
ssm = session.client('ssm')

# use path '/${AppName}/${Stage}/${AWS::Region}/Neo4jBackupMaintenanceWindowId'
maintenance_window_id = ssm.get_parameter(
    Name=f'/{APP_NAME}/{STAGE}/{AWS_REGION}/Neo4jBackupMaintenanceWindowId'
)['Parameter']['Value']

def lambda_handler(event, context):
    logger.info(json.dumps(event, indent=4))

    # parse event to determine whether to disable or enable the backup window
    alarm_state = json.loads(event['Records'][0]['Sns']['Message'])['NewStateValue']

    if alarm_state == 'ALARM':
        logger.info(f'Disabling maintenance window for backup process')
        response = ssm.update_maintenance_window(
            WindowId=maintenance_window_id,
            Enabled=False
        )
        if response['Enabled']:
            raise ValueError(f'Failed to disable maintenance window')
    elif alarm_state == 'OK':
        logger.info(f'Enabling maintenance window for backup process')
        response = ssm.update_maintenance_window(
            WindowId=maintenance_window_id,
            Enabled=True
        )
        if not response['Enabled']:
            raise ValueError(f'Failed to enable maintenance window')
    else:
        raise ValueError(f'Unknown state: {alarm_state}')
    
    return


if __name__ == "__main__":
    from pathlib import Path

    event_path = Path(__file__).resolve().parent / "enable-window-event.json"

    with open(event_path, "r") as file:
        event = json.load(file)

    lambda_handler(event,"")
