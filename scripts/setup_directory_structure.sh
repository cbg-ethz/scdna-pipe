#!/bin/bash

# Setup the directory structure for running the TuPro DNA pipeline. 
#
# Example: ./setup_directory_structure.sh /cluster/work/bewi/ngs/projects/tumorProfiler/analysis/trial_melanoma/ MADEGOD
#

CheckSuccess(){
  if [ $? -eq 0 ]; then
    echo "OK"
  else
    echo "FAIL"
  fi
}

echo "Usage: ./setup_directory_structure.sh base_path sample_id"

if [ $# -ne 2 ]
then
  echo "Please provide 2 argumens."
  exit 1
fi

base_path=$1
sample_id=$2
sample_id_dir=${sample_id}"-T"

echo
echo -n "Create base path... "
cd ${base_path}
CheckSuccess
echo 

echo -n "Create directory for the sample id... "
mkdir ${sample_id_dir}
CheckSuccess
echo

echo -n "Copy folder structure... "
cp -R /cluster/work/bewi/ngs/projects/tumorProfiler/analysis/exampleFolder/* ${sample_id_dir}/.
CheckSuccess
echo
