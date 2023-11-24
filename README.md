# The Hydropower Potential Explorer (HyPE) model
This repository provides the source codes and documentation for the HyPE model developed by Dhaubanjar et al. (2021, 2023a, 2023b) under the SustainIndus project. The HyPE model performs a spatial analysis to explore hydropower potential at user-defined river segmentation through out the basin to evaluate theoretical, technical, financial and sustainable potential. The theoretical potential module is physically-based with the hydropower potential equation applied to evaluate hydropower energy that can theoretically be generated from water flow between two grid cells. The other potential modules use a cost-minimization framework to find optimal combinations of two run-of-river hydropower plant configurations - River power plant (RP) and Diversion canal plant (DP).

The HyPE model was developed for the upper Indus basin (UIB) based on a [global hydropower assessment model](https://github.com/davidgernaat/Hydro_CostCurves) by Gernaat et al. (2019). Concepts on HyPE model development and data pre-processing for sustainable potential and regionalization to the Indus context are in [Dhaubanjar et al. (2021)](https://edu.nl/7cxcx). Detailed documentation on the cost minimization core for technical and financial potential is best found in [Gernaat et al. (2017)](https://www.nature.com/articles/s41560-017-0006-y) and [Gernaat (2019)](https://dspace.library.uu.nl/handle/1874/381146).  Results of HyPE model run for historical and future hydro-climatology of the UIB are presented in [Dhaubanjar et al. (2023a)](???) and [Dhaubanjar et al. (2023b)](10.3389/frwa.2023.1256249) respectively.

## File structure
The **Hydrus** folder contains all the model codes. All functions and scripts for data preprocessing are in **dataprep** folder. Thereafter the **theoreticalPot** folder has scripts for running the theoretical potential evaluation. Technical, financial and sustainable potential classes are run by main model files inside **Hydrus** that calls codes inside **functions** and **scripts**. This model can be run for single scenario using the ***RunOneScenario.m*** script or multiple runs for historical hydrology on HPC using scripts ***RunMainScenarios.m*** and  ***RunSAScenarios.m***. Hydropower potential assessment under future climate scenarios can be run using ***RunFutScenarios.m***. Finally, **postprocess** has the scripts for analysis of outputs and creating plots. It also contains scripts for performance evaluation of future hydropower portfolios. Refer to ***documentation/Code flow for future hydropower potential analysis.docx*** for code flow used for analysis shows in published papers. 

The **data** folder contains the input data in raw and processed forms.

The **documentation** folder contains supporting files documenting analysis workflow, model parameters, input data and model setup.

The **papers** folder contains the published papers supporting this work.

The  **output** folder contains the model outputs as well as processed figures and tables.
```
├── Hydrus
│   ├── dataprep
│   │   ├── data_prep_functions
│   │   ├── data_prep_scripts
│   │   └── downscale_Q
│   ├── functions
│   ├── postprocess
│   ├── scripts
│   └── theoreticalPot
├── data
│   ├── ASIA
│   │   └── Basin_UIB
│   ├── UI
│   │   └── data
│   └── data_prep
│       ├── UIB_outputs_LTMavgs_mmday
│       └── UIB_outputs_monTS_mmmonth
├── documentation
├── papers
└── output
```

## Installation
<!--How does the user access your project? (E.g. download, or clone with git clone…)-->
Create a copy of the project repository on your computer. If you use git, you can navigate to the location where you want to create your project on your computer. Open terminal or Git Bash from the location and clone this repository using:
```
git clone https://github.com/cicatrixx/HydropowerPotentialExplorer.git
```

If you do not use git, simply download a zip file for the repository from above:
- Click on the green Code button 
- Select "Download ZIP"
- Unzip the files in your computer into the loaction of your choice

After you have the file structure, obtain the input data from:https://surfdrive.surf.nl/files/index.php/s/Css4Dq8BYd7Ehjr. Save the data in the data/ASIA/Basin_UIB folder.

## Getting Started
Open MATLAB and navigate to the project repository location. Open the trial example case file `RunOneScenario.m` inside the Hydrus folder. Review the code and run it from the main folder where your repository sits. 

## Scenario Setup
Different scenario setup was used for the historical and future potential analysis. Under historical hydro-climatology, the HyPE model was run to explore two types of hydropower development policy scenarios. The energy focus scenarios (Large/Medium/Mixed) explored the impact of different scales of hydropower development while the geo-hazard risk representation scenarios (Risk-averse/Cost-based/Multi-hazard) evaluated three ways to represent geo-hazard risk in hydropower development policies. Additionally historical analysis quantified technical, financial and sustainable potential under full and remaining cases.

## Model Outputs
Direct model outputs can be downloaded from the ICIMOD RDS for [historical](10.26066/rds.1973690) and [future](10.26066/rds.1973697) potential datasets. Additionally, visualized potential inventory for the upper Indus compiled along side the HyPE model can be found [here]( 10.26066/rds.1973705).

### Interactive Viewers
You can interact with subbasin total results from the HyPE model in the [Indus Knowledge Partnership Platform](https://edu.nl/wqxrh and [our storymap](https://edu.nl/w3bfm). 

## Publication and Media
1. The HyPE model development and its application to the upper Indus are documented in three peer-reviewed papers.
2. Read [our storymap](https://edu.nl/w3bfm) for major highlights of the research.
3. A video is available from the presentation of the conceptual model at [AGU 2020](https://edu.nl/kph88). 
4. Check out the project website for [SustainIndus](https://edu.nl/vngvf).
 
## Citation
Please cite the HyPE model as follows:
Dhaubanjar, S., Lutz, A. F., Gernaat, D. E. H. J., Nepal, S., Smolenaars, W., Pradhananga, S., Biemans, H., Ludwig, F., Shrestha, A. B., & Immerzeel, W. W. (2021). A systematic framework for the assessment of sustainable hydropower potential in a river basin – The case of the upper Indus. Science of the Total Environment, 786. https://doi.org/10.1016/j.scitotenv.2021.147142

## Acknowledgement
This repository was developed for research conducted under the [SustainIndus](https://www.sustaindus.org/) project. The project received funding from the Netherlands Organization for Scientific Research under WOTRO Joint Sustainable Development Goals (SDG) research program (Grant W 07.30318.002). This work was partially supported by Sustainable Development Investment Portfolio (SDIP), the Department of Foreign Affairs and Trade (DFAT), Government of Australia, the Swiss Agency for Development and Cooperation (SDC) and by core funds from ICIMOD contributed by the governments of Afghanistan, Australia, Austria, Bangladesh, Bhutan, China, India, Myanmar, Nepal, Norway, Pakistan, Switzerland and the United Kingdom. 

## References
1. Dhaubanjar S, Lutz AF, Gernaat DEHJ, et al (2021) A systematic framework for the assessment of sustainable hydropower potential in a river basin – The case of the upper Indus. Sci Total Environ 786:. https://doi.org/10.1016/j.scitotenv.2021.147142

1. Gernaat DEHJ (2019) The role of renewable energy in long-term energy and climate scenarios

1. Gernaat DEHJ, Bogaart PW, Vuuren DPV, et al (2017) High-resolution assessment of global technical and economic hydropower potential. Nat Energy 2:821–828. https://doi.org/10.1038/s41560-017-0006-y
