sleepSeconds=300
nrWells=384
indir="/collab-ag-fischer/PROMISE/incoming/PROMISE-10x-4t-c-16z"
hdf5dir="/collab-ag-fischer/PROMISE/data-10x-4t-c-16z/hdf5dir"
hdf5validationdir="/collab-ag-fischer/PROMISE/data-10x-4t-c-16z/hdf5validation"
hdf5projection="/collab-ag-fischer/PROMISE/data-10x-4t-c-16z/hdf5projection"
htmldir="/collab-ag-fischer/PROMISE/data-10x-4t-c-16z/htmldir"
layoutdir="/collab-ag-fischer/PROMISE/layouts"
featuresdir="/collab-ag-fischer/PROMISE/data-10x-4t-c-16z/features"
qsub="qsub -q ag_fischer -l walltime=18:00:00,select=1:ncpus=1:mem=10GB -W umask=0011 -M jan.sauer@dkfz-heidelberg.de -m bea"
startScript="library(PROMISE);startWorkflow();warnings()"
fctcall="PROMISEworkflow"
channels=c('dsRed','FITC','DAPI')
exciteWavelengths=c('Green', 'Blue', 'UV')
fields=c('1','2','3','4')
fields_layout=c('1','2','4','3')
stacks = c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16')
hdf5_compression=3
hdf5_chunksize=256
