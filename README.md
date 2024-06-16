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

## Data: 
### One subject from the Nitrc adhd resting-state dataset, accessed via `nilearn.datasets`.
### Atlases: atlases used for parcellation of fMRI images.
### roi_list_19.csv: dataframe of ROI information from parcellation scheme including atlas, network, label, hemisphere, and coordinates. 

## Code:
### circos_functions.py: all functions for the creation of a circos plot.
### create_fc_matrices.ipynb: creates FC matrices (precision matrix and correlation matrix) from fMRI data.
### plot_corrmat.ipynb: creates Circos plot from FC correlation matrix.

## Circos plot:
### Track is a heatmap of weighted connectivity degrees.
### Links represent functional connectivty as Z-scores.


## References
#### Nitrc adhd resting-state dataset. ftp://www.nitrc.org/fcon_1000/htdocs/indi/adhd200/sites/ADHD200_40sub_preprocessed.tgz. Accessed: 2021-05-19.