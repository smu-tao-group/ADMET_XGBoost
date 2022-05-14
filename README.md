# Accurate ADMET Prediction with XGBoost
[![Python](https://img.shields.io/badge/Python-3.7+-blue.svg)](https://www.python.org)
[![LGTM](https://img.shields.io/lgtm/grade/python/github/smu-tao-group/ADMET_XGBoost.svg?style=square)](https://lgtm.com/projects/g/HTian1997/getarticle)
[![DOI](http://img.shields.io/badge/DOI-arXiv:2204.07532-B31B1B.svg)](https://arxiv.org/abs/2204.07532)

## Installation

```bash
git clone https://github.com/smu-tao-group/ADMET_XGBoost
cd ADMET_XGBoost
conda env create -f environment.yml
conda activate tdc
```

## Usage

1. Featurization: run `python featurize.py TASK_NAME` to convert SMILES to features. This step is time consuming and we provide the processed data that can be downloaded [here](https://drive.google.com/file/d/1un1kO5ZoFQ6G7WCbL0SffTiYiBon06bT/view?usp=sharing). 
2. Modeling: run `python model.py TASK_NAME` for model training and prediction. 

## Results

<table>
    <thead>
        <tr>
            <th>Tasks</th>
            <th>Evaluation</th>
            <th>Performance</th>
            <th>Rank</th>
        </tr>
    </thead>
    <tbody>
        <tr>
            <td colspan=3 style="text-align: center;">Absorption</td>
        </tr>
        <tr>
            <td>Caco2</td>
            <td>MAE</td>
            <td>0.288 &#177; 0.011</td>
            <td>1st</td>
        </tr>
        <tr>
            <td>HIA</td>
            <td>AUROC</td>
            <td>0.987 &#177; 0.002</td>
            <td>1st</td>
        </tr>
        <tr>
            <td>Pgp</td>
            <td>AUROC</td>
            <td>0.911 &#177; 0.002</td>
            <td>5th</td>
        </tr>
        <tr>
            <td>Bioav</td>
            <td>AUROC</td>
            <td>0.700 &#177; 0.010</td>
            <td>2nd</td>
        </tr>
        <tr>
            <td>Lipo</td>
            <td>MAE</td>
            <td>0.533 &#177; 0.005</td>
            <td>1st</td>
        </tr>
        <tr>
            <td>AqSol</td>
            <td>MAE</td>
            <td>0.727 &#177; 0.004</td>
            <td>1st</td>
        </tr>
        <tr>
            <td colspan=3 style="text-align: center;">Distribution</td>
        </tr>
        <tr>
            <td>BBB</td>
            <td>AUROC</td>
            <td>0.905 &#177; 0.001</td>
            <td>1st</td>
        </tr>
        <tr>
            <td>PPBR</td>
            <td>MAE</td>
            <td>8.251 &#177; 0.115</td>
            <td>1st</td>
        </tr>
        <tr>
            <td>VDss</td>
            <td>Spearman</td>
            <td>0.612 &#177; 0.018</td>
            <td>1st</td>
        </tr>
        <tr>
            <td colspan=3 style="text-align: center;">Metabolism</td>
        </tr>
        <tr>
            <td>CYP2C9 Inhibition</td>
            <td>AUPRC</td>
            <td>0.794 &#177; 0.004</td>
            <td>3rd</td>
        </tr>
        <tr>
            <td>CYP2D6 Inhibition</td>
            <td>AUPRC</td>
            <td>0.721 &#177; 0.003</td>
            <td>3rd</td>
        </tr>
        <tr>
            <td>CYP3A4 Inhibition</td>
            <td>AUPRC</td>
            <td>0.877 &#177; 0.002</td>
            <td>3rd</td>
        </tr>
        <tr>
            <td>CYP2C9 Substrate</td>
            <td>AUPRC</td>
            <td>0.387 &#177; 0.018</td>
            <td>3rd</td>
        </tr>
        <tr>
            <td>CYP2D6 Substrate</td>
            <td>AUPRC</td>
            <td>0.648 &#177; 0.023</td>
            <td>5th</td>
        </tr>
        <tr>
            <td>CYP3A4 Substrate</td>
            <td>AUPRC</td>
            <td>0.680 &#177; 0.005</td>
            <td>1st</td>
        </tr>
        <tr>
            <td colspan=3 style="text-align: center;">Excretion</td>
        </tr>
        <tr>
            <td>Half Life</td>
            <td>Spearman</td>
            <td>0.396 &#177; 0.027</td>
            <td>1st</td>
        </tr>
        <tr>
            <td>CL-Hepa</td>
            <td>Spearman</td>
            <td>0.420 &#177; 0.011</td>
            <td>2nd</td>
        </tr>
        <tr>
            <td>CL-Micro</td>
            <td>Spearman</td>
            <td>0.587 &#177; 0.006</td>
            <td>2nd</td>
        </tr>
        <tr>
            <td colspan=3 style="text-align: center;">Toxicity</td>
        </tr>
        <tr>
            <td>LD50</td>
            <td>MAE</td>
            <td>0.602 &#177; 0.006</td>
            <td>2nd</td>
        </tr>
        <tr>
            <td>hERG</td>
            <td>AUROC</td>
            <td>0.806 &#177; 0.005</td>
            <td>4th</td>
        </tr>
        <tr>
            <td>Ames</td>
            <td>AUROC</td>
            <td>0.859 &#177; 0.002</td>
            <td>1st</td>
        </tr>
        <tr>
            <td>DILI</td>
            <td>AUROC</td>
            <td>0.933 &#177; 0.011</td>
            <td>1st</td>
        </tr>
    </tbody>
</table>
