# polymer-build-workflow

A software workflow for automated, high-throughput generation of molecular
dynamics (MD) inputs for polymer systems. Works from databases of monomers
(given via SMILES) and system size specifications, without needing to provide MD
engine-specific inputs.

This is the light, data-free distribution of the build workflow: the source
shell plus a small worked example, with the study-specific datasets and
production projects removed. The full research repository, including the
monomer databases and the projects used for published work, is
[NREL_polymers](https://github.com/timbernat/NREL_polymers).

# Quickstart

After installing the environment (below), build three polymer chemistries
end-to-end from the bundled example monomer set:

```sh
bash examples/quickstart.sh
```

That runs the four steps the rest of this README documents in detail:
initialize reaction templates, write build parameters, create a Signac project
from `examples/monomers_example.csv`, and run the build. MD inputs land in
`quickstart_project/workspace/<job-id>/`. Pass a chemistry count and degree of
polymerization to vary it, e.g. `bash examples/quickstart.sh 5 10`.

`examples/monomers_example.csv` holds ten common polymers spanning six of the
eight supported mechanisms (PET, PBT, nylon-6,6, nylon-6,10, Kapton, BPA
polycarbonate, an HDI/BDO polyurethane, polystyrene, PMMA, PVAc), and doubles as
a template for the input format described under "Monomer data" below.


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

## Monomer data
Chemically, each distinct polymer chemistry in a project is encoded by its monomer feedstocks, provided as SMILES strings. Formatting for monomer dataset is (by design) very tolerant, and requires only a handful of criteria to be met to use as the basis for a polymer project. Namely, a monomer data input file must consist of:
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

Once you have supplied your monomer data file(s), you can preprocess them to ensure formatting compliance with the `src.format_data` util. This supports two formatting modes:
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

## Parameters
Non-chemical parameters for specifying system size and force field behavior are required to specify a polymer build project. For each parameter type, one can specify a set of possible values which will be swept over; these define the state space on which a polymer build project operates. To each monomer chemistry provided, a build job for every point from the Cartesian product of the sets of build parameters will be initialized for that chemistry. 
See `python -m src.parameters --help` for defaults and parameter options.

To write a new system build configuration to disc, one can run
```sh
python -m src.parameters write ...
```
with arguments specifying the values to sweep over for each parameter (E.g. can specify `--DOP 3 5 10` to indicate trimer, pentamer, and decamer versions of each chemistry should be built).

### Forcefields
This workflow current only supports [SMIRNOFF-style forcefields](https://docs.openforcefield.org/projects/toolkit/en/stable/users/smirnoff.html), as supported by the OpenFF toolkit. In principle, arbitrary common force fields (e.g. GAFF, CHARMM) can be employed, as long as one can provide a local SMIRNOFF port of the desired force field. Any locally-defined force fields should be placed in `src/forcefields/<your-ff-name>.offxml` to be recognized by the workflow; these can then be referenced in parameter files as `$FORCEFIELDS/<your-ff-name>.offxml`. Force fields without this prefix will be assumed to be shipped as flagship OpenFF forcefields, and will be searched for in any installed [`openforcefields`](https://github.com/openforcefield/openff-forcefields/tree/main/openforcefields) libraries.

Choices of combined force fields can be injected into a build workflow via the `-ffs/--forcefields` keyword of the parameters write **as a JSON-parsable string of a list of lists**, like:
To be correctly parsed, the argument following the `--forcefields` flag must be enclosed in *SINGLE QUOTES*, with each individual force field name enclosed in *DOUBLE QUOTES*. As an example, the following would be a valid argument:
```sh
python -m src.parameters write ... --forcefields '[["openff_unconstrained-2.0.0.offxml"], ["$FORCEFIELDS/diatomic_gases_IFF.offxml.offxml", "$FORCEFIELDS/CO2_TRAPPE-flex.offxml"]]'
```
Force fields grouped together by inner list will be combined into a single forcefield upon workflow run; these groups can be as small as a single force field, and will be treated as another state variable. E.g. in the example above, we define two force field groups,
* one for the unconstrained Sage 2.0.0 force field (installed and searched for in `openforcefields`) 
* one for a combination of diatomic gas and CO2 force fields, derived from IFF and TraPPE, respectively (searched for locally in `src/forcefields `)

With this specification, all polymer build jobs will be incarnated in two copies; one with the Sage FF, and the other with the composite gas force field.

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