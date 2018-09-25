import os
import logging

logging.getLogger('code').addHandler(logging.NullHandler())
FORMAT = '%(name)s:%(levelname)s:  %(message)s'
logging.basicConfig(format=FORMAT)
#log.info('importing WriteRotateFileHandler...')
#from utils.filestuff import WriteRotateFileHandler

__all__ = ['operators','evodevo','rencode','utils','arn']

#check for SGE environment vars
if not os.getenv('JOB_NAME'):
    os.environ['JOB_NAME'] = __name__
if not os.getenv('SGE_TASK_ID'):
    os.environ['SGE_TASK_ID'] = '0'

notice = "CoDe  Copyright (C) 2014 Rui L. Lopes\n" + \
         "This program comes with ABSOLUTELY NO WARRANTY.\n"

print notice
