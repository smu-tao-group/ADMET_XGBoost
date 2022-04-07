# Accurate Prediction of ADMET properties with XGBoost
[![Python](https://img.shields.io/badge/Python-3.7+-blue.svg)](https://www.python.org)
[![LGTM](https://img.shields.io/lgtm/grade/python/github/smu-tao-group/ADMET_XGBoost.svg?style=square)](https://lgtm.com/projects/g/HTian1997/getarticle)

## Installation

```bash
git clone https://github.com/smu-tao-group/ADMET_XGBoost
cd ADMET_XGBoost
conda env create -f environment.yml
conda activate tdc
```

## Usage

1. Featurization: run `python featurize.py TASK_NAME` to convert SMILES to features. This step is time consuming and we provide the processed data that can be downloaded here. 
2. Modeling: run `python model.py TASK_NAME` for model training and prediction. 

## Results

<table>
    <thead>
        <tr>
            <th>Tasks</th>
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
            <td>0.291 &#177; 0.015</td>
            <td>1st</td>
        </tr>
        <tr>
            <td>HIA</td>
            <td>0.988 &#177; 0.002</td>
            <td>1st</td>
        </tr>
        <tr>
            <td>Pgp</td>
            <td>0.911 &#177; 0.001</td>
            <td>5th</td>
        </tr>
        <tr>
            <td>Bioav</td>
            <td>0.7 &#177; 0.01</td>
            <td>2nd</td>
        </tr>
        <tr>
            <td>Lipo</td>
            <td>0.536 &#177; 0.002</td>
            <td>2nd</td>
        </tr>
        <tr>
            <td>AqSol</td>
            <td>0.734 &#177; 0.006</td>
            <td>1st</td>
        </tr>
        <tr>
            <td colspan=3 style="text-align: center;">Distribution</td>
        </tr>
        <tr>
            <td>BBB</td>
            <td>0.907 &#177; 0.002</td>
            <td>1st</td>
        </tr>
        <tr>
            <td>PPBR</td>
            <td>8.252 &#177; 0.19</td>
            <td>1st</td>
        </tr>
        <tr>
            <td>VDss</td>
            <td>0.627 &#177; 0.009</td>
            <td>1st</td>
        </tr>
        <tr>
            <td colspan=3 style="text-align: center;">Metabolism</td>
        </tr>
        <tr>
            <td>CYP2C9 Inhibition</td>
            <td>0.769 &#177; 0.0</td>
            <td>3rd</td>
        </tr>
        <tr>
            <td>CYP2D6 Inhibition</td>
            <td>0.717 &#177; 0.001</td>
            <td>3rd</td>
        </tr>
        <tr>
            <td>CYP3A4 Inhibition</td>
            <td>0.872 &#177; 0.005</td>
            <td>3rd</td>
        </tr>
        <tr>
            <td>CYP2C9 Substrate</td>
            <td>0.383 &#177; 0.012</td>
            <td>3rd</td>
        </tr>
        <tr>
            <td>CYP2D6 Substrate</td>
            <td>0.636 &#177; 0.012</td>
            <td>5th</td>
        </tr>
        <tr>
            <td>CYP3A4 Substrate</td>
            <td>0.677 &#177; 0.007</td>
            <td>1st</td>
        </tr>
        <tr>
            <td colspan=3 style="text-align: center;">Excretion</td>
        </tr>
        <tr>
            <td>Half Life</td>
            <td>0.383 &#177; 0.058</td>
            <td>2nd</td>
        </tr>
        <tr>
            <td>CL-Hepa</td>
            <td>0.395 &#177; 0.016</td>
            <td>4th</td>
        </tr>
        <tr>
            <td>CL-Micro</td>
            <td>0.588 &#177; 0.003</td>
            <td>2nd</td>
        </tr>
        <tr>
            <td colspan=3 style="text-align: center;">Toxicity</td>
        </tr>
        <tr>
            <td>LD50</td>
            <td>0.601 &#177; 0.004</td>
            <td>1st</td>
        </tr>
        <tr>
            <td>hERG</td>
            <td>0.779 &#177; 0.01</td>
            <td>4th</td>
        </tr>
        <tr>
            <td>Ames</td>
            <td>0.859 &#177; 0.0</td>
            <td>1st</td>
        </tr>
        <tr>
            <td>DILI</td>
            <td>0.925 &#177; 0.012</td>
            <td>1st</td>
        </tr>
    </tbody>
</table>
