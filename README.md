# SSG-LUGIA
## Single Sequence based Genome Level Unsupervised Genomic Island Prediction Algorithm

This repository contains the original implementation of SSG-LUGIA, an unsupervised learning based tool to predict genomic islands.

<a href="https://colab.research.google.com/drive/1sT7jEKJMgDVCHjkxX1vt-utx9GyoU7cY?usp=sharing" target="_blank"><img src="https://colab.research.google.com/assets/colab-badge.svg"></a>



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

7. Alternatively, the model name can be omitted and the user can set the parameters interactively

```SSG_LUGIA(sequence_fasta_file_path='sample_data/NC_003198.1.fasta')```

8. Alternatively, the user can input a custom model as dictionary

```SSG_LUGIA(sequence_fasta_file_path='sample_data/NC_003198.1.fasta',model_parameters=custom_model)```

9. Alternatively, the user can create a model based on their requirement, save it as a json file and input the path to the json file

```SSG_LUGIA(sequence_fasta_file_path='sample_data/NC_003198.1.fasta',model_name='path-to-json')```


## Model Parameters

SSG-LUGIA combines several sequence based features to infer GIs using an unsupervised anomaly detection pipeline. The various model parameters can be found in [SSG-LUGIA Model Parameters](https://nibtehaz.github.io/SSG-LUGIA/params). Users can develop custom model variants by changing these parameters and also save the model as json for future use.


## Citation Request

If you use ***SSG-LUGIA*** in your project, please cite the following paper

```
@article{ibtehaz2021ssg,
  title={SSG-LUGIA: Single Sequence based Genome Level Unsupervised Genomic Island Prediction Algorithm},
  author={Ibtehaz, Nabil and Ahmed, Ishtiaque and Ahmed, Md Sabbir and Rahman, M Sohel and Azad, Rajeev K and Bayzid, Md Shamsuzzoha},
  journal={Briefings in Bioinformatics},
  year={2021}
}
```
