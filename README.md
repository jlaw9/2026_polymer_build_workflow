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

# Config
## Reaction templates
This toolkit ships with 8 classes of polymerization mechanism pre-defined by default, namely:
* Polyesters
* Polyamides
* Polyimides
* Polycarbonates
* Polycarbonates (non-phosgene route)
* Polyvinyls
* Polyuerthanes
* Polyurethanes (non-isocyanate)
  
These are defined in [src.reactions](./src/reactions.py), and can be appended to insert other mechanisms not included here, if such chemistries are of interest. For details on how to define these reaction inputs, see the [`polymerist` reaction examples](https://github.com/timbernat/polymerist_examples/tree/main/1-polymerization)

Once you're satisfied with the mechanisms defined, initialized the SMARTS templates for these reaction definitions by running:
```sh
bash src/reactions.py
```

## Parameters
System size and force field parameters for each system build job are configured in [src.parameters](src/parameters.py), and are broken down into two types:
* ParametersSwept: each of these fields is a range of values which will be iterated over in all combinations (i.e. in Cartesian product). These include:
  *  Degree of polymerization (peroligomer)
  *  Max number of atoms (per box)
  *  Partial charge method
* ParametersConfig: these are single parameters common to all jobs and are related to configuring MD parameters, including choice of Sage forcefield and nonbonded cutoffs.

Once you;re satisfied with the parameters, initialize and cache them by running:
```sh
bash src/reactions.py
```

# Supplying chemical data
Finally, provide your monomer data as a csv with the column containing you monomer smiles strings as in a column labelled. <add details on formatting>

On first-time use, the parameter sets and reaction mechanism templates defined for your project need to be initialized. This is as easy as running:
```sh
bash format_mdat.sh # this is actually pretty specific to the PolyID study; want to extend before final publication
```

# Project
Once you've specified your project config and chemical data, create a new project using:
```sh
```