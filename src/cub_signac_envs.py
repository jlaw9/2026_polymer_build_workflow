'''Custom ComputeEnvironment definitions for CU Boulder supercomputing'''

from flow.environment import DefaultSlurmEnvironment, get_environment, _PartitionConfig


class CURCEnvironment(DefaultSlurmEnvironment):
    '''Environment boilerplate common to University of Colorado Boulder Research Computing clusters'''
    hostname_pattern = r'.*\.rc(\.int)?\.colorado\.edu' # also works on interactive nodes
    template = 'curc.sh'

    @classmethod
    def add_args(cls, parser):
        '''Inject Blanca-specific directive into submit command'''
        super().add_parser(parser)
        parser.add_argument(
            '--nodes',
            type=int,
            default=1,
            help='Number of compute nodes (bundles of cores) requested',
        )
        parser.add_argument(
            '--ntasks_per_node',
            type=int,
            default=1,
            help='Number of CPU/GPU cores to request on each node',
        )
        parser.add_argument(
            '-t',
            '--time',
            type=float,
            help='Job walltime, in number of hours (as float)'
        )
        parser.add_argument(
            '--gres',
            choices=['gpu'],
            default=None,
            help='Optional generic request for service; used to request GPUs',
        )
        parser.add_argument(
            '--mail-user',
            default=None,
            help='The email address to which job event emails should be sent'
        )
        parser.add_argument(
            '--mail-type',
            choices=[
                'NONE',
                'BEGIN',
                'END',
                'FAIL',
                'REQUEUE',
                'ALL',
            ],
            default='NONE',
            help='Optional, event on which to send an email to the user',
        )

class CUAlpineEnvironment(CURCEnvironment):
    '''
    Environment for the University of Colorado Boulder Alpine cluster
    
    https://curc.readthedocs.io/en/latest/clusters/alpine/alpine-hardware.html
    '''
    hostname_pattern = r'[^b].*\.rc(\.int)?\.colorado\.edu' # overtly DOESN'T start with "b", also works on interactive nodes
    template = 'alpine.sh'

    # partition
    _partition_config = _PartitionConfig(
        cpus_per_node={
            'amilan' : 32,
            'amem'   : 48, 
            'ami100' : 64,
            'aa100'  : 64,
            'csu'    : 32,  
            'amc'    : 64,
            'acompile' : 64,
            'atesting' : 64,
            'atesting_a100' : 64,
            'atesting_a100mig' : 64,
            'atesting_mi100' : 64,
        },
        gpus_per_node={
            'ami100' : 24,
            'aa100'  : 33,
            'amc'    : 12,
        },
    )

    @classmethod
    def add_args(cls, parser):
        '''Inject Blanca-specific directive into submit command'''
        super().add_parser(parser)
        parser.add_argument(
            '--account',
            default='ucb-general',
            help='The name of the allocation being submitted under',
        )
        parser.add_argument(
            '--qos',
            choices=[
                'normal',
                'long',
                'mem',
            ],
            default='normal',
            help='Quantlity of Service modifier for requesting special config e.g. longer runtime',
        )

class CUBlancaShirtsEnvironment(CURCEnvironment):
    '''Environment for the Shirts Research Group's allocation on the University of Colorado Boulder Blanca condo cluster'''
    hostname_pattern = r'b.*\.rc(\.int)?\.colorado\.edu' # starts with "b" also works on interactive nodes
    template = 'blanca_shirts.sh'

    @classmethod
    def add_args(cls, parser):
        '''Inject Blanca-specific directive into submit command'''
        super().add_parser(parser)
        parser.add_argument(
            '--account',
            default='blanca-shirts',
            help='The name of the allocation being submitted under',
        )
        parser.add_argument(
            '--qos',
            default='blanca-shirts',
            help='Quantlity of Service modifier for requesting special config e.g. longer runtime',
        )
        parser.add_argument(
            '--gres',
            choices=['gpu'],
            default='gpu',
            help='Optional generic request for service; used to request GPUs',
        )



test_hostnames = [
    'login13.rc.colorado.edu',              # Alpine login node
    'c3cpu-c15-u32-3.rc.int.colorado.edu',  # Ainteractive
    'bgpu-shirts3.rc.int.colorado.edu',     # Blanca OnDemand
    'bgpu-shirts1.rc.int.colorado.edu',     # Blanca sbatch job
]

if __name__ == '__main__':
    curr_env = get_environment()
    print(curr_env.__name__, curr_env.template)