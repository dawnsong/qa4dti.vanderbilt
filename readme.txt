Author:   Carolyn B. Lauzon ( clauzon1@gmail.com )
                Medical Image Analysis and Statistical Interpretation Lab
                Vanderbilt University
                Nashville TN, 37235
                https://masi.vuse.vanderbilt.edu/index.php/Main_Page
  Creation Date:  August 07, 2012
  Last Modified:  August 07, 2012
%  For nitric release, version 1.0
% please see https://masi.vuse.vanderbilt.edu/index.php/QA_DTI

 The open source programs provided herein are freely provided by the Dr. Bennett Landman MASI lab at Vanderbilt University: https://masi.vuse.vanderbilt.edu/index.php/Main_Page
 Several re-distributed open source programs are part of this release. The authors thank these contributions and suggest visititing  https://masi.vuse.vanderbilt.edu/index.php/QA_DTI for a larger reference list.
Empirical data for this simulation comes from the Multi-Modal MRI Reproducibility Study [1] available online at: http://www.nitrc.org/projects/multimodal
 [1] Bennett. A. Landman, Alan J. Huang, Aliya Gifford, Deepti S. Vikram, Issel Anne L. Lim, Jonathan A.D. Farrell, John A. Bogovic, Jun Hua, Min Chen, Samson Jarso, Seth A. Smith, Suresh Joel, Susumu Mori, James J. Pekar, Peter B. Barker, Jerry L. Prince, and Peter C.M. van Zijl. \'93Multi-Parametric Neuroimaging Reproducibility: A 3T Resource Study\'94, NeuroImage. (2010) NIHMS/PMC:252138 doi:10.1016/j.neuroimage.2010.11.047\

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %%% Inventory of Demo_DTI_QA_120807 tar file %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The zip file contains five files/folders. (1) this readme.txt, (2) release.txt, (3) runme.m, (4) demo_data, and (5) QA_DTI
(1) this readme.txt contains the instructions for the program
(2) release.txt with release information
(3) runme.m is the user called function to run the pipeline
(4) demo_data contains the demonstration data inputs and expected-outputs in two folders. 
     (4.1) inputs - this folder contains input data from scan number 23 of the multi-modal study cited above.
     (4.2) expected_outputs - this folder contains two items:
             (4.2.1) KKI*.pdf-pipeline report for demo data
             (4.2.2) 'Matlab_terminal_output.txt'-terminal outputs from running the pipeline on the demo data. This examples typical 'errors/warnings' that may occur even for a successful dataset.         
(5) QA_DTI contains the programs necessary for the pipeline. The folder contains two .m functions and three subfolders.
    (5.1) QA_DTI_Pipeline.m - the main program for the pipeline
    (5.2) QAmfiles - this folder contains additional matlab programs and is broken into two sub-folders
              (5.3.1) homemade_mfiles - these are files written by this author or others editing the pipeline program
              (5.3.2) matlab_central_files - these are community written open-source files downloaded from various sources, including mathworks 
    (5.3) multi-atlas - this folder contains the files for multi-atlas segmentation, as well as necessary functions from the nifti toolbox from mathworks


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           %%% Instructions  %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Requirements %%%
(1) This program makes system calls to a UNIX command line and therefore requires UNIX/LINUX/MAC O.S.
(2) This program requires installation of CAMINO: http://cmic.cs.ucl.ac.uk/camino/index.php?n=Main.HomePage
(3) This program requires installation of FSL: http://www.fmrib.ox.ac.uk/fsl/
(5) This program has been tested with CAMINO version date 2009-09-03 (latest version as of August 07,2012) 
    FSL version 4.1.4 and 
    Matlab 7.11 (R2010). Other versions may effect program performance.
(6) This program processes large DTI datasets. It has been tested with 12 GB of RAM. Minimum RAM requirement is dataset dependent and has not been determined for demo-data.

%%% Installation %%%
(1) Unzip file in desired location. 
(2) FSL and CAMINO both require editing of the .bashrc file to define global variables. Follow software package installation instructions.
(3) CAMINO heap size may require adjusting. The minimum possible heap size for this program is data dependent and has not been determined for demo-data. 
    If possible, set to 8000. Either edit your .bashrc file OR in the shell that you will run QA_DTI from type:
                 export CAMINO_HEAP_SIZE=8000 
   **This command is what our lab uses but may not be compatible with other versions of CAMINO. plese see CAMINO website.
(4) Matlab JAVA heap size may require adjusting.The minimum possible JAVA heap space is data dependent and has not been determined for the demo-data. 
     4a. You can follow instructions found at mathworks: http://www.mathworks.com/support/solutions/en/data/1-18I2C/
     4b. We run our program with a java.opts file alotting 8000mb to Java heap space. If this size is possible for your system you can create a java.opts textfile in the Matlab home directory with the following one line:
                -Xmx8000m
     ** This optino may not be compatible with different matlab versions. Please see maltab website. 

%%% Running QA pipeline on Included Demo Data %%%
Note: Demo data is subject 23 from the multi-modal study referenced above. The gradient table has been changed for compatability with CAMINO. Gradient table and b-value vector have been adjusted and renamed for use with the QA_DTI pipeline.
(1) open runme.m
(2) edit line-4 with the full path name to the demo data: ('/full/path/to/Demo_DTI_QA_120807/demo_data/inputs/KKI2009-23-DTI.nii')
(3) edit line-5 with the full path name to the user specified location for the pipeline output folder. ('/full/path/to/my_outputs')
(4) save changes
(5) In the matlab command line, initiate the pipeline by typing 'runme'
(6) output data is stored in user specified output folder. Ouput PDF can be compared to the expected PDF which is stored in Demo_DTI_QA_120807/demo_data/expected_outputs.

%%% Running QA pipeline on user specified DTI Dataset %%%
Note: We have tried to make the program flexible to different orientation structures, however validation has only been for a limited subset of image orientations. See 'Troubleshooting' for suggestions.
(1) program accepts par/rec (PAR/REC) and nifti files.
(2) Nifti files.  
     2a. Image orientation and gradient orientation of nifti must be compatible with CAMINO. 
         This can be verified by running data directly with CAMINO following the CAMINO DTI tutorial: 
            http://cmic.cs.ucl.ac.uk/camino/index.php?n=Tutorials.DTI
     2b. You must follow the QA pipeline naming convention for nifti files. If your nifti file is named filename.nii, 
         then your gradient table must be stored as filename_bvecs with each gradient direction corresponding to 1-row, 
         and your b-value vector stored as a single column in the file named filename_bvals
(3) open runme.m
(4) edit line-4 with the full path name to your data: (e.g. '/full/path/to/mydata.par')
(3) edit line-5 with the full path name to the location of the pipeline output folder. ('/full/path/to/my_outputs')
(6) save changes
(7) open QA_DTI/QA_DTI_Pipeline.m
(8) edit line 17 'n_bo=5';  Set variable no_bo to the number of averaged Bo volumes used in your data aquisition. If unknown, set n_bo=1.
    The n_bo parameter only effects bias estimation. n_bo=1 will likely be the 'worst case' as more modern DTI protocols include an averaged Bo.
(9) save changes. In the matlab command line, initiate the pipeline by typing 'runme'
(10) output data is stored in output_folder.

%%% Troubleshooting %%%
(1)	To troubleshoot, first double check the correctiness of your input file in a robust viewing environment (such as MIPAV). 
      Especially check image voxel resolution.
                       image orientation (make sure L/R A/P S/I are labelled correctly in viewer). 
(2) If it exists, check correctness of pipeline output registered image, RegIm.nii, in robust viewing environment (such as MIPAV).
      RegIm.nii is located in either ('/full/path/to/my_outputs/temp_folder') OR ('/full/path/to/my_outputs/extra')
(2) If you use par/rec data, you can try editing the function QA_DTI_nitrc01/QA_DTI/QAmfiles/homemade/getParinfoV4.m
(3) If you use nii data, you can try editing the function QA_DTI_nitrc01/QA_DTI/QAmfiles/homemade/getNIIinfo.m
(4) For faster runtimes during trouble-shooting, open QA_DTI/DTI_QA_Pipeline_nitrc01.m Comment 'NORMAL RUN VALUES' (l. 11-14). Uncomment 'TROUBLESHOOTING VALUES' (l. 6-8).
(5) visit for more information:  https://masi.vuse.vanderbilt.edu/QAserver/
    
