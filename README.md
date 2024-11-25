# PET Pre-processing Tool

A Shiny R application for automated preprocessing of neuroimaging data. This tool supports PET and MRI scans, streamlining brain extraction, MNI template registration, MRI-PET co-registration, COMBAT harmonization, and RAVEL normalization.

## Features

- **A Shiny R Application**: A tool for automated neuroimaging preprocessing with interactive visualization, artifact inspection and error handling capabilities for smooth execution, with clear and informative messages displayed in the output console for user guidance.
  ![Tool Overview](assets/images/Basic_UI.png)
  ![Tool Overview with Error Handling Messages](assets/images/Error_Handling.png)
- **Artifact Viewer**: Launch an external NIfTI artifact viewer for detailed examination of imaging data to include/exclude pre-processing of unwanted neuroimaging data.
  ![Artifact Viewer Example](assets/images/NIfTI_Viewer_Dashboard.png)
- **Interactive Visualization**: Inspect preprocessing results with slice-by-slice navigation and real-time updates.
- **Brain Extraction**: Automated skull-stripping of MRI scans to isolate brain tissue.
  ![Brain Extraction Example](assets/images/Brain_Extraction.png)
- **MNI Template Registration**: Align MRI scans to the MNI standard space for consistent spatial analysis.
  ![MNI Template Registration Example](assets/images/MRI_Template_Registration.png)
- **MRI-PET Co-registration**: Aligns PET scans with corresponding MRI scans for integrated analysis.
  ![MRI-PET Co-Registration Example](assets/images/MRI-PET_CoRegistration.png)
- **COMBAT Harmonization**: Reduces scanner and site effects in neuroimaging data using the Combating Batch Effects (COMBAT) algorithm.
  ![COMBAT Harmonization Example](assets/images/COMBAT_Harmonization.png)
- **RAVEL Normalization**: Normalizes PET scan intensities to reduce unwanted variations using the Removing inter-subject technical variability in magnetic resonance imaging studies (RAVEL) algorithm.
  ![RAVEL Normalization Example](assets/images/RAVEL_Normalization.png)


- ## Folder Structure

- `templates/`: Contains MNI template files for registration.
- `Dashboard/`: Contains a python code to display and export NIfTI image properties.
- `Output_Brain_Extraction/`: Outputs from the brain extraction step.
- `Output_MNI_Template_Registration/`: Outputs from the MNI registration step.
- `Output_MRI_PET_Co-Registration/`: Outputs from the MRI-PET co-registration step.
- `Output_COMBAT_Harmonization/`: Outputs from COMBAT harmonization.
- `Output_RAVEL_Normalization/`: Outputs from RAVEL normalization.

## Dependencies

### R Packages
- **Data Processing**: `neurobase`, `oro.nifti`, `extrantsr`, `neuroCombat`, `RAVEL`, `ANTsR`
- **Visualization**: `shiny`, `grid`, `rgl`, `png`
- **Utilities**: `fslr`, `shinyFiles`, `reticulate`, `readxl`

### External Tools
- **FSL**: Required for image registration using tools like FLIRT.
- **ANTs**: Used for advanced normalization techniques.
- **Python 3**: For the NIfTI artifact viewer and Python-based integrations.

## Installation

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/yourusername/yourrepository.git
2. 

