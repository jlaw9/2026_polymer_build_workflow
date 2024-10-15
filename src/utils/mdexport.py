'''Wrappers for exporting parameterized systems to MD file outputs'''

from typing import Optional

from pathlib import Path
from openmm import Integrator, Context

from openff.interchange import Interchange
from openff.interchange.interop.openmm._positions import to_openmm_positions
from openff.interchange.components.mdconfig import MDConfig

from polymerist.mdtools.openmmtools import serialization
from polymerist.mdtools.openmmtools.forcegroups import impose_unique_force_groups


def interchange_to_lammps(interchange : Interchange, lmp_data_path : Path, lmp_input_path : Path) -> None:
    '''Produce LAMMPS input and data files from an OpenFF Interchange'''
    interchange.to_lammps(lmp_data_path) # MD data file
    mdc = MDConfig.from_interchange(interchange)
    # mdc.write_lammps_input(lmp_input_path) # input directive file
    mdc.write_lammps_input(input_file=lmp_input_path, interchange=interchange) # input directive file

    # replacing generic lmp file with data file from above
    with lmp_input_path.open('r') as in_file:
        in_file_block = in_file.read()

    with lmp_input_path.open('w') as in_file:
        in_file.write(
            in_file_block.replace('out.lmp', f'"{lmp_data_path}"') # need surrounding double quotes to allow LAMMPS to read special symbols in filename (if present)
        )

def interchange_to_openmm(interchange : Interchange, integrator : Integrator, omm_top_path : Path, omm_sys_path : Path, omm_state_path : Path, state_params : Optional[dict[str, bool]]=None) -> Context:
    '''Produce OpenMM System and State .xml files from an OpenFF Interchange'''
    if state_params is None:
        state_params = {
        'getPositions'  : True,
        'getVelocities' : True,
        'getForces'     : True,
        'getEnergy'     : True,
        'getParameters' : True,
        'getParameterDerivatives' : True,
        'getIntegratorParameters' : True
    }
    
    system  = interchange.to_openmm(combine_nonbonded_forces=False)
    topology = interchange.to_openmm_topology()
    positions = to_openmm_positions(interchange, include_virtual_sites=True)

    impose_unique_force_groups(system)
    context = Context(system, integrator)
    context.setPositions(positions)

    ## writing OpenMM files
    serialization.serialize_system(omm_sys_path, system)
    serialization.serialize_state_from_context(omm_state_path, context, state_params=state_params)
    serialization.serialize_openmm_pdb(omm_top_path, topology, positions)

    return context