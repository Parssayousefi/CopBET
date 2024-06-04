from bids import BIDSLayout
from nilearn import datasets
from nilearn import input_data
from nilearn.interfaces.fmriprep import load_confounds
import os
import pandas as pd
import numpy as np
#import scipy.io
#import matlab.engine


# def process_data_bids(bids_root, strategy, atlas_name):
#     print("BIDS root:", bids_root)
#     print("Directories and files at BIDS root:", os.listdir(bids_root))

#     atlases = {
#         'schaefer400': datasets.fetch_atlas_schaefer_2018( n_rois=400),
#         'schaefer1000': datasets.fetch_atlas_schaefer_2018(n_rois=1000),
#         "yeo": datasets.fetch_atlas_yeo_2011()
        
#     }
#     if atlas_name not in atlases:
#         raise ValueError(f"Atlas not supported. Available options: {', '.join(atlases.keys())}")

#     atlas_data = atlases[atlas_name]
#     if atlas_name == 'yeo':
#         # Debugging line to check available keys in the returned atlas data
#         print(atlas_data.keys())
#         # Ensure this key exists in the atlas_data dictionary
#         atlas_filename = atlas_data['thin_17']
#     else:
#         atlas_filename = atlas_data.maps  # Fetching the atlas filename
    
#     # Initialize BIDS Layout with derivatives enabled
#     layout = BIDSLayout(bids_root, validate=False, derivatives=True, absolute_paths=True)
#     #print("Subjects found:", layout.get_subjects())

#     masker = input_data.NiftiLabelsMasker(labels_img=atlas_filename, standardize=True, verbose=2)

#     all_data = [] 
#     for subject_id in layout.get_subjects(): #iterating through all the subjects
#         print("Processing subject:", subject_id)
#         #fetching the functional files
#         func_files = layout.get(subject=subject_id, suffix='bold', extension='nii.gz', return_type='filename')
        
#         func_file = func_files[0] #fetching the first functional file

#         if strategy == 'gsr': #if the strategy is global signal regression
            
#             confounds = load_confounds(func_file, strategy=['global_signal',"high_pass", "wm_csf"])
#             print(type(confounds[0]))
#             print(confounds[0].shape)

#         elif strategy == 'compcor': #if the strategy is compcor
#             confounds = load_confounds(func_file, strategy=['motion',"high_pass", "wm_csf"], motion= "basic", compcor='anat_combined', n_compcor='all')
#             print(type(confounds[0]))
#             print(confounds[0].shape)

#         else: #if the strategy is not supported
#             raise ValueError(f"Unknown strategy: {strategy}", "Available options: gsr, compcor")

#         time_series = masker.fit_transform(func_file, confounds=confounds[0]) #extracting the time series data

#         df = pd.DataFrame(time_series, columns=[f'Region_{i}' for i in range(time_series.shape[1])]) #creating a dataframe with the time series data
        
#         df['Subject'] = subject_id #adding the subject id to the dataframe

#         all_data.append(df) #appending the dataframe to the list of all data
#     return pd.DataFrame() if not all_data else pd.concat(all_data, ignore_index=True)

# bids_root = r"C:\Users\prsyu\OneDrive\Bidlung\University\M.S. Leiden University\M.S. Neuroscience (Research)\Thesis\CopBET\CopBET_Python\nilearn_data\development_fmri"

# strategy = 'compcor'
# atlas_name = "yeo"
# #results_yeo = process_data_bids(bids_root, strategy , atlas_name = 'yeo')
# results = process_data_bids(bids_root, strategy, atlas_name)

# print(results.head)

results = pd.DataFrame()
def df_to_matlab_table(df, eng):
    """Converts a pandas DataFrame to a MATLAB table."""
    # Convert DataFrame to dictionary of lists for MATLAB compatibility
    mat_data = {col: matlab.double(df[col].tolist()) for col in df.columns if df[col].dtype in [float, int]}
    mat_data.update({col: eng.cell(df[col].tolist()) for col in df.columns if df[col].dtype == object})
    
    # Transfer dictionary to MATLAB and convert to table
    mat_dict = eng.struct(**mat_data)
    mat_table = eng.struct2table(mat_dict)
    return mat_table

def save_data_for_matlab(df, save_dir):
    subjects = df['Subject'].unique()
    
    # Start MATLAB engine
    eng = matlab.engine.start_matlab()
    
    for subject in subjects:
        subject_data = df[df['Subject'] == subject].drop('Subject', axis=1)
        mat_table = df_to_matlab_table(subject_data, eng)
        mat_filename = os.path.join(save_dir, f'subject_{subject}_data.mat')
        
        # Save the MATLAB table to a .mat file
        eng.save(mat_filename, 'mat_table')
        print(f"Data for {subject} saved to {mat_filename}")
    
    # Stop MATLAB engine
    eng.quit()

# After processing all data
save_dir = r"C:\Users\prsyu\OneDrive\Bidlung\University\M.S. Leiden University\M.S. Neuroscience (Research)\Thesis\CopBET\CopBET_Python\Leiden_datasets_plug"
save_data_for_matlab(results, save_dir)