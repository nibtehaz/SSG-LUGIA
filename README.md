# SSG-LUGIA
## Single Sequence based Genome Level Unsupervised Genomic Island Prediction Algorithm

This repository contains the original implementation of SSG-LUGIA, an unsupervised learning based tool to predict genomic islands.

## Project Website

[SSG-LUGIA Project Website](https://nibtehaz.github.io/SSG-LUGIA)

## Publication

## Codes

The codes for SSG-LUGIA are written in python and can be found [here](https://github.com/nibtehaz/SSG-LUGIA/tree/master/codes)


## Requirements

**Note** : In the latest version of Scikit-Learn the implementation of *EllipticEnvelope* has been changed, so please use the specified version to obtain reproducible results.

> * numpy==1.17.0
> * biopython==1.70
> * tqdm==4.19.5
> * scikit-learn==0.19.1


## Usage

1. Clone this repository

```$ git clone https://github.com/nibtehaz/SSG-LUGIA.git```

2. Install the requirements

```$ pip3 install -r requirements.txt```

3. Navigate to the ```/codes``` directory

4. Launch Python CLI

``` $ python3 ```

5. Import the SSG-LUGIA pipeline

```from main import SSG_LUGIA```

6. Execute it with a genome sequence fasta file and a standard model name from ```SSG-LUGIA-F```, ```SSG-LUGIA-R```, ```SSG-LUGIA-P```


```SSG_LUGIA(sequence_fasta_file_path='sample_data/NC_003198.1.fasta',model_name='SSG-LUGIA-F')```
