## Setup virtual enviroment
```bash
python -m venv venv
source venv/bin/activate
```

## Run after updating the code if you want to run it locally
```bash
python -m pip install --upgrade pip
pip install -r requirements.txt
```

## Code to create a Circos plot of brain functional connectivity from an ADHD fMRI resting-state dataset.

![circos_plot](https://github.com/DoraBeslic/Brain-Connectivity-Viz/assets/122387780/c57ef3d0-87b2-4926-9f03-46e9d203cca5)

### Data: 
- One subject from the **Nitrc adhd resting-state dataset** [[1]](#1), accessed via `nilearn.datasets`.
- **Atlases**: Atlases used for parcellation of fMRI images.
- **roi_list_198.csv**: Dataframe of ROI information from parcellation scheme including atlas, network, label, hemisphere, and coordinates. 

### Code:
- **circos_functions.py**: All functions for the creation of a circos plot.
- **create_fc_matrices.ipynb**: Creates FC matrices (precision matrix and correlation matrix) from fMRI data.
- **plot_corrmat.ipynb**: Creates Circos plot from FC correlation matrix.

### Circos plot:
- Track is a heatmap of weighted connectivity degrees.
- Links represent functional connectivity as Z-scores.

### References
<a id="1">[1]</a>
Nitrc adhd resting-state dataset. ftp://www.nitrc.org/fcon_1000/htdocs/indi/adhd200/sites/ADHD200_40sub_preprocessed.tgz. Accessed: 2021-05-19.

