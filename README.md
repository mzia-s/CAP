# CAP

**Title** - CAP: Commutative Algebra Prediction of Protein-Nucleic Acid Binding Affinities.

**Authors** - Mushal Zia, Faisal Suwayyid, Yuta Hozumi, JunJie Wee, Hongsong Feng, and Guo-Wei Wei.

---

## Table of Contents

- [Table of Contents](#table-of-contents)
- [Introduction](#introduction)
- [Model Interpretability](#model-Interpretability)
- [Prerequisites](#prerequisites)
- [Datasets](#datasets)
- [Modeling with PSRT-based features](#Modeling-with-PSRT-based-features)
    - [Generation of PSRT-based features for protein-nucleic acid complex](#II-Generation-of-PSRT-based-features-for-protein-ligand-complex)

- [Results](#results)
- [License](#license)
- [Citation](#citation)

---

## Introduction

An accurate prediction of protein-nucleic acid binding affinity is vital for deciphering genomic processes, yet existing approaches often struggle in reconciling high accuracy with interpretability and computational efficiency. In this study, we introduce CAP, a sequence-centric framework which couples persistent Stanley-Reisner theory with advanced representation learning for predicting protein-nucleic acid binding affinities. CAP encodes proteins through transformer-learned embeddings that retain long-range evolutionary context and represents DNA with $\textit{k}$-mer-based topological embeddings derived from persistent facet ideals, which capture fine-scale nucleotide geometry. We demonstrate that CAP surpasses the SVSBI protein-nucleic acid benchmark \cite{shen2023svsbi} and, in a further test, maintains reasonable performance on a newly curated protein-RNA dataset. Leveraging only primary sequences, CAP generalizes to any protein-nucleic acid pair with minimal preprocessing, enabling genome-scale analyses without 3D structural data and promising faster pathways for drug discovery and protein engineering.

> **Keywords**: Persistent commutative algebra, facet persistence barcodes, persistent ideals, commutative algebra learning, protein-nucleic acid binding.

---

## Model Architecture

An illustration of the filtration process of the persistent commutative algebra is shown below.

![Model Architecture](scheme.png)

Further explain the details in the [paper](https://github.com/WeilabMSU/CAML), providing context and additional information about the architecture and its components.

---

## Prerequisites

- numpy                     1.21.0
- scipy                     1.7.3
- pytorch                   1.10.0 
- pytorch-cuda              11.7
- torchvision               0.11.1
- scikit-learn              1.0.2
- python                    3.10.12
- biopandas                 0.4.1
--- 

## Datasets

A brief introduction about the benchmark datasets.

| Datasets                |Total    | Training Set                 | Test Set                                             |
|-|-----------------------------|------------------------------|------------------------------                        |
| PDBbind-v2016       |4057 [data](./dataset)|3767  [data](./dataset)                        | 290 [data](./dataset)                         |
| Metalloprotein-ligand       |2463 [data](./dataset)|1845  [data](https://weilab.math.msu.edu/Downloads/PSRT/PDBbind.zip)                        | 618 [data](./dataset)                         |


- PDBbind-v2016: the protein-ligand complex structures. Download from [PDBbind database](http://www.pdbbind.org.cn/)
- Metalloprotein-ligand: the metalloprotein-ligand complex structures were complied from PDBbind-v2020 ([PDBbind database](http://www.pdbbind.org.cn/)) by [Jiang2023]
- data: the .csv file, which contains the protein ID and corresponding binding affinity for PDBbind data.
---

## Modeling with PSRT-based features

### I. Build machine learning models using PSRT-based features.
```shell
python src/build_model.py --dataset_name v2016
```
### II. Generation of sequence-based features for protein and small molecules
Protein sequence embeddings were generated with [Transformer Protein language model ESM2](https://github.com/facebookresearch/esm) [Rives2021].

Small molecular sequence embeddings were generated with [Transformer small molecule language model](https://github.com/WeilabMSU/PretrainModels) [Chen2021]. The input small molecular sequence is SMILES string. Instructions on molecular descriptors using this language model is provided in the github.

## Results

### I. Modeling the (Metallo)protein-ligand datasets

|Datasets                                        | Training Set                  | Test Set| PCC | RMSE (kcal/mol) |
|-------------------------------------------------|-------------                  |---------|-    |-                |
| PDBbind-v2016 [result](./Results)      |  3767 | 290 | 0.858 |  1.673|
| Metalloprotein-ligand [result](./Results) |1845| 618 | 0.745/0.755 |  1.947|


Note, twenty gradient boosting regressor tree (GBRT) models were built for each dataset with 20 indenpedent runs with different random numbers. The PSRT-based features and transformer-based features were paired with GBRT, respectively. The predictions can be found in the [results](./Results) folder. 

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Citation

- [Jiang2023] Dejun Jiang, Zhaofeng Ye, Chang-Yu Hsieh, Ziyi Yang, Xujun Zhang, Yu Kang, Hongyan Du, Zhenxing Wu, Jike Wang, Yundian Zeng, et al. Metalprognet: a structure-based deep graph model for metalloprotein–ligand interaction predictions. Chemical Science, 14(8):2054–2069, 2023.
- [Chen2021] Dong Chen, Jiaxin Zheng, Guo-Wei Wei, and Feng Pan. Extracting predictive representations from hundreds of millions of molecules. The Journal of Physical Chemistry Letters, 12(44):10793–10801, 2021.
- [Rives2021] Rives, Alexander, Joshua Meier, Tom Sercu, Siddharth Goyal, Zeming Lin, Jason Liu, Demi Guo et al. "Biological structure and function emerge from scaling unsupervised learning to 250 million protein sequences." Proceedings of the National Academy of Sciences 118, no. 15 (2021): e2016239118.
---
