import os

__all__ = ['operators','evodevo','rencode','utils','arn']

#check for SGE environment vars
if not os.getenv('JOB_NAME'):
    os.environ['JOB_NAME'] = __name__
if not os.getenv('SGE_TASK_ID'):
    os.environ['SGE_TASK_ID'] = '0'
