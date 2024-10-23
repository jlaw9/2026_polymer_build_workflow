'''Custom ComputeEnvironment definitions for CU Boulder supercomputing'''

from flow.environment import DefaultSlurmEnvironment, get_environment


class CUAlpineEnvironment(DefaultSlurmEnvironment):
    '''Environment for the University of Colorado Boulder Alpine cluster'''
    hostname_pattern = r'[^b].*\.rc(\.int)?\.colorado\.edu' # overtly DOESN'T start with "b", also works on interactive nodes
    template = 'alpine.sh'

class CUBlancaEnvironment(DefaultSlurmEnvironment):
    '''Environment for the University of Colorado Boulder Alpine cluster'''
    hostname_pattern = r'b.*\.rc(\.int)?\.colorado\.edu' # starts with "b" also works on interactive nodes
    template = 'blanca.sh'

test_hostnames = [
    'login13.rc.colorado.edu',              # Alpine login node
    'c3cpu-c15-u32-3.rc.int.colorado.edu',  # Ainteractive
    'bgpu-shirts3.rc.int.colorado.edu',     # Blanca OnDemand
    'bgpu-shirts1.rc.int.colorado.edu',     # Blanca sbatch job
]

if __name__ == '__main__':
    curr_env = get_environment()
    print(curr_env.__name__, curr_env.template)