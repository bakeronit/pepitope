# pEpitope Calculator

[![Python](https://img.shields.io/badge/python-3.11%2B-blue)](https://www.python.org/downloads/) [![GitHub Pages](https://img.shields.io/badge/demo-GitHub%20Pages-blue)](https://bakeronit.github.io/pepitope)  [![License](https://img.shields.io/github/license/bakeronit/pepitope)](https://github.com/bakeronit/pepitope/blob/main/LICENSE)

A Python implementation of the [pEpitope Calculator](https://www.mathworks.com/matlabcentral/fileexchange/83698-pepitope-calculator), originally developed in MATLAB by [Melia Bonomo](https://meliabonomo.com/portfolio/influenza/)

This tool estimates influenza vaccine efficacy/effectiveness by calculating the pEpitope distance between vaccine strains and dominant circulating strains based on hemagglutinin (HA) protein sequences.


## Features

- Calculate pEpitope scores between vaccine and circulating strains
- Predict vaccine efficacy (VE) with standard error estimates
- Identify the dominant epitope and amino acid substitutions
- Support for CDC and Northern Hemisphere (NH) epidemiological datasets
- Both CLI and web interface (via PyScript)


## Installation

Requires Python 3.11+

```bash
# Clone the repository
git clone https://github.com/bakeronit/pepitope.git
cd pepitope

# Install dependencies
pip install -e .
# or using uv
uv sync
```

## Usage

### Command Line Interface

```bash
python cli.py <dominant_fasta> <vaccine_fasta> [OPTIONS]
```

**Arguments:**
- `dominant_fasta`: Path to FASTA file containing the dominant circulating strain HA sequence
- `vaccine_fasta`: Path to FASTA file containing the vaccine strain HA sequence

**Options:**
- `--vaccine-type`: Influenza subtype (default: `H3N2`)
- `--dataset`: Epidemiological dataset, `CDC` or `NH` (default: `CDC`)
- `--mode`: Prediction mode, `Efficacy` or `Effectiveness` (default: `Efficacy`)
- `--outfmt`: Output format, `txt` or `json` (default: `txt`)

**Example:**

```bash
python cli.py test/test1.fasta test/vaccine.fasta --dataset CDC --mode Efficacy
```

### Web Interface

Open [pepitope](https://bakeronit.github.io/pepitope) in a browser to use the PyScript-powered web interface. Paste HA sequences directly into the text areas and configure options via the UI.

## Input Format

Input files must be in FASTA format containing the full hemagglutinin (HA) protein sequence:

```
>A/Switzerland/9715293/2013
MKTIIALSYILCLVFAQKIPGNDNSTATLCLGHHAVPNGTIVKTITNDRIEVTNATELVQ...
```

## References

> M.E. Bonomo, R.Y. Kim, M.W. Deem, “Modular epitope binding predicts influenza quasispecies dominance and vaccine effectiveness: Application to 2018/19 season,” Vaccine, 37(24): 3154-3158, 2019.