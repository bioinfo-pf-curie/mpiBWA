#! /bin/bash

type mpiBWAIdx

if [[ ! $? -eq 0 ]]
then
    echo "ERROR: mpiBWAIdx not found in your PATH"
    exit 1
fi

output_dir=${HOME}/mpiBWAExample

if [[ ! -d ${output_dir} ]]
then
    mkdir -p ${output_dir}
    echo "INFO: the directory ${output_dir} has been created to store the results generated by mpiBWAExample"
fi


echo "INFO: untar the (small) reference human genome with bwa index"
tar zxvf data/hg19.small.tar.gz --directory ${output_dir}

echo "INFO: creating the binary image of the (small) reference genome with mpiBWAIdx"
mpiBWAIdx  ${output_dir}/hg19.small.fa

