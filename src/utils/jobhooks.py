'''Commonly-used project-level operation hooks for Signac Projects (https://docs.signac.io/en/latest/hooks.html)'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Callable, Optional
import time

import json
from pathlib import Path

from signac.job import Job
from signac.project import Project


# JSON FORMATTING
def reindent_job_json_factory(filename : str) -> Callable[[Job, Optional[int]], None]:
    '''Build a function which accepts a Job and indents the file'''
    if not filename.endswith('.json'): # TODO - add support/handling for Path objects
        raise ValueError('Filename must point to file with .json extension')
    
    def reindent_job_json(job : Job, indent : Optional[int]=4) -> None:
        # prechecks
        if (indent is not None) and indent < 0:
            raise ValueError(f'Must provide positive indent amount (not {indent})')
        
        if not job.isfile(filename):
            raise FileNotFoundError(job.fn(filename)) # TODO - find way to log this failure
            
        # read/re-write
        with open(job.fn(filename), 'r') as file:
            data = json.load(file)

        with open(job.fn(filename), 'w') as file:
            json.dump(data, file, indent=indent)
    
    return reindent_job_json

reindent_json_statepoint : Callable[[Job, Optional[int]], None] = reindent_job_json_factory(Job.FN_STATE_POINT)
reindent_json_document   : Callable[[Job, Optional[int]], None] = reindent_job_json_factory(Job.FN_DOCUMENT)

# JOB HOOKS
class ProjectHooks:
    def __init__(
        self,
        operation_times_attr : str='operation_times_sec',
        indent_amount : Optional[int]=4,
    ) -> None:
        self.operation_times_attr = operation_times_attr
        self.indent_amount = indent_amount
        
    # OPERATION TIMING
    def operation_start_time_hook(self) -> Callable[[str, Job], None]:
        '''Define a hook which records operation start times'''
        def record_operation_start_time(operation_name : str, job : Job) -> None:
            '''Record into the job document when a particular operation began'''
            job.doc.setdefault(self.operation_times_attr, {})
            job.doc[self.operation_times_attr].update({f'{operation_name}_start' : time.time()})
        return record_operation_start_time
        
    def operation_duration_hook(self) -> Callable[[str, Job], None]:
        '''Define a hook which records operation completion times, assuming that a start time has been recorded'''
        def record_operation_duration(operation_name : str, job : Job) -> None:
            '''Record into the job document how long a particular operation took'''
            op_times = job.doc.get(self.operation_times_attr, {})
            start_time = op_times.pop(f'{operation_name}_start') # look up and withdraw start time
            if start_time is None: # NOTE: dicts in Signac documents are NOT pure Python dicts; their "pop()" method returns None by default and never raises KeyError
                raise ValueError(f'No start time recorded for operation "{operation_name}"; cannot calculate operation duration')
            
            job.doc[self.operation_times_attr].update({operation_name : (time.time() - start_time)})
        return record_operation_duration
    
    # JSON FORMATTING
    def operation_indent_statepoint_hook(self)-> Callable[[str, Job], None]:
        def indent_statepoint(operation_name : str, job : Job) -> None:
            '''Apply indentation to Job statepoint output once an operation completes'''
            reindent_json_statepoint(job, indent=self.indent_amount)
        return indent_statepoint

    def operation_indent_document_hook(self)-> Callable[[str, Job], None]:
        def indent_document(operation_name : str, job : Job) -> None:
            '''Apply indentation to Job document output once an operation completes'''
            reindent_json_document(job, indent=self.indent_amount)
        return indent_document
    
    # HOOK INSTALLATION
    def install_hooks(self, project : Project) -> Project:
        '''Inject hooks defined at the project level'''
        project.project_hooks.on_start.extend([
            self.operation_start_time_hook(),
        ])
        project.project_hooks.on_exit.extend([
            self.operation_duration_hook(), # for now, record times for failing and succeeding operations the same
            self.operation_indent_document_hook(),
            self.operation_indent_statepoint_hook(),
        ]) 
        # project.project_hooks.on_success = ...
        # project.project_hooks.on_exception = ... # TODO: add distinct timing mode for cases where operations fail before exit

        return project
