#! /bin/bash
#
# runme4quest.sh
#% Init: 2013-11-12 19:53
#% Copyright (C) 2013~2020 Xiaowei.Song <dawnwei.song@gmail.com>
#% http://restfmri.net
#% Distributed under terms of the AFL (Academy Free license).
#
CALLDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
elog(){ echo "$@" 1>&2 ; }
 echo "$(date) | $0 '"$@"'" >> $(basename $0).LOG
trim(){ trimmed=$@
    trimmed="${trimmed#"${trimmed%%[![:space:]]*}"}"   # remove leading whitespace characters
    trimmed="${trimmed%"${trimmed##*[![:space:]]}"}"   # remove trailing whitespace characters
    echo $trimmed ; }
usage(){ printf "
Usage: (%s)
    ${0} [options] x.dcm
    ###${0} [options] x.nii x.qa4dti.output
    x.nii/x_bvecs/x_bvals must be paired exist at the same path
    options:
" "$(grep -m1 '#% Version:' `readlink -f ${BASH_SOURCE[0]}`)" 1>&2 ; }
if [ $# -eq 0 ]; then usage; exit 0; fi 
verbose=0;  
while [ $# -gt 0 ]; do 
    case "$1" in
        '-v') verbose=1 ; shift ;;    
        '-wd') WORKINGDIR=$2 ; shift 2 ;;    
        '-*') elog 'Unknown parameters'; usage; exit 0 ;;
        *) break ;;
    esac
done

runxdummy -x && source ~/.xdummy
scrot $(getmp ).scr.png

dcm2nii -d n -g n -p n $1  
nii=$(ls *.nii)
prfx=${nii%.nii}
mv ${prfx}.bvec ${prfx}_bvecs
mv ${prfx}.bval ${prfx}_bvals
mv ${prfx}.nii o${prfx}.nii

outdir=$(readlink -f ${WORKINGDIR:-qa4dti.vanderbilt.output})
mkdir -p $outdir

export FSLOUTPUTTYPE=NIFTI_GZ
matlab -nodesktop -nosplash <<eom
QA_pathname='/home/xst833/fmri/vanderbilt/qa4dti/QA_DTI'
java_path=sprintf('%s%smulti-atlas%smasi-fusion%sbin%s',QA_pathname,filesep,filesep,filesep,filesep)

addpath(genpath(QA_pathname));
javaaddpath(java_path);

reslice_nii('o${prfx}.nii', '${prfx}.nii');

infile='$(readlink -f ${prfx}.nii)'; %.nii/_bvals/_bvecs must be paried
outdir='$outdir';


%Run
DTI_QA_Pipeline_page12(infile, outdir, QA_pathname);
eom

