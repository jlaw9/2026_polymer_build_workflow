A software workflow for automated, high-throughput generation of molecular dynamics (MD) inputs for polymer systems. Works with on databases of monomers (given via SMILES) and system size specifications without needing to provide MD engine-specific inputs

# Installation
It is assumed you have access to a package and environment manager like `mamba` for this installation
To acquire, clone this repo and create a compatible virtual environment viz:
```sh
git clone https://github.com/timbernat/NREL_polymers/
cd NREL_polymers
mamba env create -f nrel_reqs.yml
mamba activate nrel-polymers
```
If installing on HPC, your may require use of a higher-memory node to complete the environment solve in this installation.

# Configuring chemical and structural inputs
## Reaction templates
This toolkit ships with 8 classes of polymerization mechanism pre-defined by default, namely:
* Polyesters
* Polyamides
* Polyimides
* Polycarbonates
* Polycarbonates (non-phosgene route)
* Polyvinyls
* Polyurethanes
* Polyurethanes (non-isocyanate)
  
These are defined in [src.reactions](./src/reactions.py), and can be appended to insert other mechanisms not included here, if such chemistries are of interest. For details on how to define these reaction inputs, see the [`polymerist` reaction assembly tutorials](https://github.com/timbernat/polymerist_examples/tree/main/1-polymerization)

Once you're satisfied with the mechanisms defined, initialized the SMARTS templates for these reaction definitions by running:
```sh
python -m src.reactions
```

## Parameters
System size and force field parameters for each system build job are configured in [src.parameters](src/parameters.py), and are broken down into two types:
* ParametersSwept: each of these fields is a range of values which will be iterated over in all combinations (i.e. in Cartesian product). These include:
  *  Degree of polymerization (peroligomer)
  *  Max number of atoms (per box)
  *  Partial charge method
* ParametersConfig: these are single parameters common to all jobs and are related to configuring MD parameters, including choice of Sage forcefield and nonbonded cutoffs.

Once you're satisfied with the parameters, initialize and cache them by running:
```sh
python -m src.parameters
```

## Supplying monomer data
Formatting for monomer data inputs is (by design) very tolerant, and requires only a handful of criteria to be met to use as the basis for a polymer project. Namely, a monomer data input file must consist of:
* A tabular file in either .csv or .xlsx format
* Containing one column (field) titled any of the following:
  * `smiles_original`
  * `smiles_monomer`
  * `monomer_smiles`
  * `monomer`
  * `monomers`
  * `Monomer`
  * `Monomers`
* With records whose value for that field consist of either
  * A tuple of SMILES strings for each distinct monomer (e.g. for PET, have `('COC(=O)c1ccc(cc1)C(=O)OC', 'OCCO')`)
  * A single SMILES string with a [disconnection](https://www.daylight.com/meetings/summerschool98/course/dave/smiles-disco.html) (single period character) separating the distinct monomers (e.g. for PET, have `'COC(=O)c1ccc(cc1)C(=O)OC.OCCO'`)

Any additional fields in each record of the data file can contain arbitrary data related to the monomer preparation e.g. name for resulting polymer, expected polymer density, labelled mechanism of polymerization, etc. These additional fields are transferred to the `document` portion of any job acting on the specified monomer chemistry in that record.

Once you have supplied you monomer data file(s), you can preprocess them to ensure formatting compliance with the `src.format_data` util. This supports two formatting modes:
* `Merge`: combines one or more data files into a single, formatted "master" file
* `Sequential`: takes one of more data files and formats each separately into the same number of formatted datafiles

Monomer data files can be identified by either of the following mutually-exclusive options:
* Names of files: pass as list after `-mdat`/`--monomer-paths` flag
* Search pattern: pass as file regex pattern after `-g`/`--glob` flag

Many formatting utilities are supplied, including uniqufication of chemistry, subselection of data, etc. For more details on formatting options, run
```sh
python -m src.format_data merge --help
python -m src.format_data sequence --help
```

As a quickstart example, if you have only one monomer data file, you should run in `sequential` mode as:
```sh
python -m src.format_data sequence -mdat <path to datafile> -od <output-directory> --postfix "fmt"
```

# Project
## Initialization
Once you've specified your project config and chemical data, create a new project using:
```sh
src/init_signac.py -mdat <path_to_monomer_data> --project-name <pick_a_name> -od <output_directory> # TODO: specify statepoint parameters config file explicitly?
```

Optionally, you may want to create a project from only a subset of the provided chemistries. You may do this either by selecting the first N monomer chemistries via:
```sh
src/init_signac.py -mdat <path_to_monomer_data> --project-name <pick_a_name> -od <output_directory> -num $n_sampled
```

or selecting N chemistries randomly via:
```sh
src/init_signac.py -mdat <path_to_monomer_data> --project-name <pick_a_name> -od <output_directory> -num $n_sampled --random
```

## Managing a project
The operations available to you on a polymer build project are largely the same as those on any other Signac FlowProject. Complete documentation on available commands is available in the [Signac-flow docs on command line interface](https://signac.readthedocs.io/projects/flow/en/latest/project-cli.html), but we highlight a few of the most common here for reference.

To check the status of a project, run
```sh
python -m src.project -path <path_to_project> -rxns <path_to_rxns_json> status
```
where `<path_to_project>` is the path to the directory created during project initialization and `<path_to_rxns_json>` is the path to the JSON file created during reaction initilaization.

To run structure build jobs locally, run:
```sh
python -m src.project -path <path_to_project> -rxns <path_to_rxns_json> run -o <opgrp>
```
Valid operation group names (to be substituted for `<opgrp>` above) for the PolymerBuildProject defined here are, in order):
* everything (runs all of below)
* validate_chemistry
* perceive_mechanism
* fragment
* oligomerize
* pack_lattice
* to_interchange
* md_export
  * openmm_export
  * lammps_export

Operations later in this list cannot be requested to run before operations earlier in the list

To submit build jobs to cluster, run:
```sh
python -m src.project \
  -path <path_to_project> \
  submit \
  --parallel \
  --bundle 6 \
  -o <opgrp> \
  <any other args you define in your ComputeEnvironment template> \
  --job-output <logfilename.out> \
```

*NOTE*: the above submission requires you define a ComputeEnvironment compatible with your local HPC cluster in [src/templates](src/templates); see [Signac-flow docs on ComputeEnvironments](https://signac.readthedocs.io/en/latest/environments.html) for more info

# Other utils
Beyond managing parameter config, project setup, status tracking, and job submission, the source code provided with this repository also provides you with a few other useful utilities including:

## Automated ring piercing detection
See `python -m src.analyze_piercing --help`

## Bond length distribution compilation and plotting
See `python -m src.bond_lengths --help`

## Extracting subset of jobs matching criterion to new project
See `python -m src.sieve --help`

## Polymer SMILES compilation for machine learning training
See `python -m src.compile_ML_smiles --help`