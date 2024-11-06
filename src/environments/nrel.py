'''Custom ComputeEnvironment definitions for the National Renewable Energy Laboratory (NREL) supercomputing'''

from flow.environment import DefaultSlurmEnvironment, get_environment, _PartitionConfig


# TODO: get SLURM scheduler directives and node specs from PolyID team
class NRELEnvironment(DefaultSlurmEnvironment):
    '''Environment boilerplate common to NREL HPC clusters'''
    hostname_pattern = r'.*hpc\.nrel\.gov'
    template = 'nrel.sh'

class KestrelEnvironment(NRELEnvironment):
    '''
    Environment for the NREL Kestrel cluster
    https://nrel.github.io/HPC/Documentation/Systems/Kestrel/
    '''
    hostname_pattern = r'.*\.kestrel.*?hpc\.nrel\.gov'
    template = 'kestrel.sh'

class VermillionEnvironment(NRELEnvironment):
    '''
    Environment for the NREL Vermillion cluster
    https://nrel.github.io/HPC/Documentation/Systems/Vermilion/
    '''
    hostname_pattern = r'.*\.vermillion.*?hpc\.nrel\.gov'
    template = 'vermillion.sh'

class SwiftEnvironment(NRELEnvironment):
    '''
    Environment for the NREL Swift cluster
    https://nrel.github.io/HPC/Documentation/Systems/Swift/
    '''
    hostname_pattern = r'.*\.swift.*?hpc\.nrel\.gov'
    template = 'swift.sh'