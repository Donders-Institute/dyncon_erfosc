#!/bin/sh

export FREESURFER_HOME=/opt/freesurfer/5.3
export SUBJECTS_DIR=$1

source $FREESURFER_HOME/SetUpFreeSurfer.sh

# rest of autorecon2
recon-all -autorecon2 -subjid $2
recon-all -fill -subjid $2
recon-all -tessellate -subjid $2
recon-all -smooth1 -subjid $2
recon-all -inflate1 -subjid $2
recon-all -qsphere -subjid $2
recon-all -fix -subjid $2
recon-all -white -subjid $2
recon-all -smooth2 -subjid $2
recon-all -inflate2 -subjid $2

# running a full autorecon3 requires the file rawavg.mgz to exist
cp $SUBJECTS_DIR/mni_resliced.mgz $SUBJECTS_DIR/rawavg.mgz

recon-all -autorecon3 -subjid $2

#recon-all -gcareg -subjid $2
#recon-all -canorm -subjid $2
#recon-all -careg -subjid $2
#recon-all -careginv -subjid $2
#recon-all -calabel -subjid $2
#recon-all -normalization2 -subjid $2
#recon-all -maskbfs -subjid $2
#recon-all -segmentation -subjid $2
#recon-all -fill -subjid $2

#recon-all -tessellate -subjid $2
#recon-all -smooth1 -subjid $2
#recon-all -inflate1 -subjid $2
#recon-all -qsphere -subjid $2
#recon-all -fix -subjid $2
##cp brain.mgz brain.finalsurfs.mgz
##recon-all -finalsurfs -subjid $2
#recon-all -white -subjid $2
#recon-all -smooth2 -subjid $2
#recon-all -inflate2 -subjid $2
#recon-all -sphere -subjid $2
#recon-all -surfreg -subjid $2
