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
    ${0} [options] x.nii x.qa4dti.output
    x.nii/x_bvecs/x_bvals must be paired exist at the same path
    options:
" "$(grep -m1 '#% Version:' `readlink -f ${BASH_SOURCE[0]}`)" 1>&2 ; }
if [ $# -eq 0 ]; then usage; exit 0; fi 
verbose=0;  
while [ $# -gt 0 ]; do 
    case "$1" in
        '-v') verbose=1 ; shift ;;    
        '--template') T1TPL=$2 ; shift 2 ;;    
        '-*') elog 'Unknown parameters'; usage; exit 0 ;;
        *) break ;;
    esac
done

runxdummy -x && source ~/.xdummy
scrot $(getmp ).scr.png



infile=$(readlink -f $1);
outdir=$(readlink -f ${2:-qa4dti.vanderbilt.output})
mkdir -p $outdir

export FSLOUTPUTTYPE=NIFTI_GZ
matlab -nodesktop -nosplash <<eom
QA_pathname='/home/xst833/fmri/vanderbilt/qa4dti/QA_DTI'
java_path=sprintf('%s%smulti-atlas%smasi-fusion%sbin%s',QA_pathname,filesep,filesep,filesep,filesep)

addpath(genpath(QA_pathname));
javaaddpath(java_path);

infile='$infile'; %.nii/_bvals/_bvecs must be paried
outdir='$outdir';


%Run
DTI_QA_Pipeline(infile, outdir, QA_pathname);
eom

