# Codes that reproduce the computational model of Bilash et al manuscript.

Enter `_codes` directory.

## mechanisms

This directory contains all .mod files that are used in the hoc neuronal models.

## poisson inputs

Python script that generates theta-filtered poisson inputs.


## Python scripts

### `{cellname}_cell.py`: validation of the {cell_name} cell

### `{cellname}_synapses.py`: validation of the {cell_name} synapses

`{cellname}` can take the values:
- pyramidal
- cck
- olm
- vipcr
- vipcck


### The rest files are described below:

#### analysis_canonical_circuit.py
Analysis script that analyzes the voltage traces of the single presynaptic spike simulation.

#### analysis_interneurons_synapses.py
Script that analyses the outputs from `{interneuron_name}_synapses.py` scritps.

#### analysis_poisson_inputs.py
Analysis script that analyzes the voltage traces of the poisson inputs simulation.

#### analysis_single_cell.py
Analysis script that analyzes the voltage traces of the single cell experiment.

#### canonical_circuit.py
Script that simulates the canonical circuit. Stimulation of all EC at t=700ms.

#### canonical_circuit_poisson_inputs.py
Script that simulates the canonical circuit. Stimulation of all EC with theta-filtered poisson inputs.

#### canonical_ccsv_files_.py
Transform the output in csv files to import them on GraphPad.

#### cell_models.py
Classes of all neuronal models used here.

#### opt.py
Plot functions and NEURON initialization function.

#### plots_poisson_inputs.py
Plotting of poisson input experiments.

#### soma_spiking_analysis.py
Somatic firing and ISI for the poisson input experiment.

#### synaptic_metrics.py
Function that calculate all synaptic metrics, such as rise time and decay time.



Author: Spiros Chavlis, PhD (schavlis [AT] imbb [DOT] forth [DOT] gr)