'''For creating a clone of a PolymerBuildProject containin gonly those jobs which successfully produced LAMMPS files'''

import logging
logging.basicConfig(level=logging.INFO)

from argparse import ArgumentParser, Namespace
from shutil import copytree
from pathlib import Path

# TODO: expand this to support passing in an arbitrary condition
try: # call as python module
    from .project import PolymerBuildProject, exported_to_lammps
except ImportError: # call as script file
    from project import PolymerBuildProject, exported_to_lammps


def prepare_LAMMPS_project(args : Namespace) -> None:
    '''Prechecks and project setup for successful-LAMMPS-only co-projects'''
    if not args.postfix:
        raise ValueError('Must provide non-empty postfix for cloned project')
    
    # setup clone project directory and initialize project
    lmp_proj_path : Path = Path(f'{args.project_path.stem}_{args.postfix}')
    lmp_proj_path.mkdir(exist_ok=args.allow_overwrites)
    lmp_project = PolymerBuildProject.init_project(lmp_proj_path)

    # copy valid jobs to new project
    orig_project = PolymerBuildProject.get_project(args.project_path) # this will raise Exception if an invalid path is provided
    n_jobs_found : int = 0
    for job in orig_project:
        if exported_to_lammps(job):
            logging.info(f'Found LAMMPS files for job "{job.id}"')
            copytree(job.path, f'{lmp_project.workspace}/{job.id}')
            n_jobs_found += 1
    logging.info(f'Found and copied {n_jobs_found} job(s) satifsying LAMMPS file requirement')
    logging.info(f'Results can be found in {lmp_proj_path.resolve()}')

def main() -> None:
    parser = ArgumentParser()
    parser.add_argument(
        '-proj',
        '--project_path',
        type=Path,
        required=True,
        help='The path to the pre-existing project from which jobs should be drawn',
    )
    parser.add_argument(
        '-pf',
        '--postfix',
        default='LAMMPS',
        help='The postfix to attach to the original project\'s name to enforce uniqueness',
    )
    parser.add_argument(
        '-aow',
        '--allow_overwrites',
        action='store_true',
        help='Whether to permit creation of a new project on top of one which already exists'
    )
    
    args = parser.parse_args()
    prepare_LAMMPS_project(args)

if __name__ == '__main__':
    main()