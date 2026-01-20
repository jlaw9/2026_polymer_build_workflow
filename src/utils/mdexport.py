'''Wrappers for exporting parameterized systems to MD file outputs'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from typing import Optional
from pathlib import Path

from openmm import Integrator, Context, XmlSerializer

from openff.interchange import Interchange
from openff.interchange.interop.openmm._positions import to_openmm_positions

from polymerist.mdtools.openmmtools import serialization
from polymerist.mdtools.openmmtools.forces import impose_unique_force_groups


def interchange_to_openmm(
    interchange : Interchange,
    integrator : Integrator,
    omm_top_path : Path,
    omm_sys_path : Path,
    omm_state_path : Path,
    omm_integ_path : Path,
    state_params : Optional[dict[str, bool]]=None
) -> Context:
    '''Produce OpenMM System and State .xml files from an OpenFF Interchange'''
    # sanitizing inputs
    if isinstance(omm_top_path, str):
        omm_top_path = Path(omm_top_path) # TODO: add xml extension assertion
    if isinstance(omm_sys_path, str):
        omm_sys_path = Path(omm_sys_path) # TODO: add xml extension assertion
    if isinstance(omm_state_path, str):
        omm_state_path = Path(omm_state_path) # TODO: add xml extension assertion
    if isinstance(omm_integ_path, str):
        omm_integ_path = Path(omm_integ_path) # TODO: add xml extension assertion

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
    
    # creating OpenMM simulation components
    system  = interchange.to_openmm(combine_nonbonded_forces=False)
    topology = interchange.to_openmm_topology()
    positions = to_openmm_positions(interchange, include_virtual_sites=True)

    impose_unique_force_groups(system)
    context = Context(system, integrator)
    context.setPositions(positions)

    # writing OpenMM files
    with open(omm_integ_path, 'w') as integ_file:
        integ_file.write(XmlSerializer.serialize(integrator))

    with open(omm_sys_path, 'w') as sys_file:
        sys_file.write(XmlSerializer.serialize(system))

    serialization.serialize_state_from_context(omm_state_path, context, state_params=state_params)
    serialization.serialize_openmm_pdb(omm_top_path, topology, positions)

    return context