%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% User inputs %%
QA_pathname='/home/dawnsong/test/qa.dti/vanderbilt/Demo_DTI_QA_120807/QA_DTI';  % path to QA_DTI, end with QA_DTI, not QA_DTI/ 
input_file = '/home/dawnsong/test/qa.dti/vanderbilt/Demo_DTI_QA_120807/demo_data/inputs/KKI2009-23-DTI.nii'; % location of input DTI data, either filename.nii OR filename.par (filename.PAR)
output_folder = '/home/dawnsong/test/qa.dti/vanderbilt/Demo_DTI_QA_120807/output.test'; % location for desired outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% ADD path to java and matlab 
java_path=sprintf('%s%smulti-atlas%smasi-fusion%sbin%s',QA_pathname,filesep,filesep,filesep,filesep);

addpath(genpath(QA_pathname))
javaaddpath(java_path)


%% call DTI QA function
DTI_QA_Pipeline(input_file, output_folder,QA_pathname)
