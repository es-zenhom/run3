#!/bin/bash

# Check if analyzer argument is provided
if [ -z "$1" ]; then
    echo "Error: Analyzer not specified. Please provide the analyzer as the first argument."
    exit 1
fi

# Assign the provided arguments
analyzer=$1         # e.g. "Run2AllHad3FJ"
output_dir_tag=$2   # optional tag for output directory name (e.g. the date)

function runPreselection() {

    # Set the output directory path
    output_dir="output_"${analyzer}
    if [[ -n "$output_dir_tag" ]]; then
        output_dir="${output_dir}_${output_dir_tag}"
    fi
    output_dir="/ceph/cms/store/user/$USER/analysisNtuples/"${output_dir}"/"

    # Configuration files
    config_file_sig="input/input_sig_jg.json" 
    config_file_bkg="input/input_bkg_jg.json"
    config_file_data="input/input_data_jg.json"

    # Run the analysis commands in the background
    ./bin/runAnalysis -n 2 -i ${config_file_sig} -o sig --outdir ${output_dir} --ana ${analyzer} --cutflow &
    ./bin/runAnalysis -n 2 -i ${config_file_data} -o data --outdir ${output_dir} --ana ${analyzer} --cutflow &
    ./bin/runAnalysis -n 2 -i ${config_file_bkg} -o bkg --outdir ${output_dir} --ana ${analyzer} --cutflow &

    ./bin/runAnalysis -n 2 -i ${config_file_sig} --met --var up -o sig_met_unclustered_up --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --met --var down -o sig_met_unclustered_down --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_Absolute --var up -o sig_jec_absolute_up --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_Absolute --var down -o sig_jec_absolute_down --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_Absolute_2016preVFP --var up -o sig_jec_absolute_2016preVFP_up --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_Absolute_2016preVFP --var down -o sig_jec_absolute_2016preVFP_down --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_Absolute_2016postVFP --var up -o sig_jec_absolute_2016postVFP_up --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_Absolute_2016postVFP --var down -o sig_jec_absolute_2016postVFP_down --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_Absolute_2017 --var up -o sig_jec_absolute_2017_up --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_Absolute_2017 --var down -o sig_jec_absolute_2017_down --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_Absolute_2018 --var up -o sig_jec_absolute_2018_up --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_Absolute_2018 --var down -o sig_jec_absolute_2018_down --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_BBEC1 --var up -o sig_jec_bbec1_up --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_BBEC1 --var down -o sig_jec_bbec1_down --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_BBEC1_2016preVFP --var up -o sig_jec_bbec1_2016preVFP_up --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_BBEC1_2016preVFP --var down -o sig_jec_bbec1_2016preVFP_down --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_BBEC1_2016postVFP --var up -o sig_jec_bbec1_2016postVFP_up --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_BBEC1_2016postVFP --var down -o sig_jec_bbec1_2016postVFP_down --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_BBEC1_2017 --var up -o sig_jec_bbec1_2017_up --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_BBEC1_2017 --var down -o sig_jec_bbec1_2017_down --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_BBEC1_2018 --var up -o sig_jec_bbec1_2018_up --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_BBEC1_2018 --var down -o sig_jec_bbec1_2018_down --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_EC2 --var up -o sig_jec_ec2_up --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_EC2 --var down -o sig_jec_ec2_down --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_EC2_2016preVFP --var up -o sig_jec_ec2_2016preVFP_up --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_EC2_2016preVFP --var down -o sig_jec_ec2_2016preVFP_down --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_EC2_2016postVFP --var up -o sig_jec_ec2_2016postVFP_up --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_EC2_2016postVFP --var down -o sig_jec_ec2_2016postVFP_down --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_EC2_2017 --var up -o sig_jec_ec2_2017_up --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_EC2_2017 --var down -o sig_jec_ec2_2017_down --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_EC2_2018 --var up -o sig_jec_ec2_2018_up --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_EC2_2018 --var down -o sig_jec_ec2_2018_down --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_FlavorQCD --var up -o sig_jec_flavorqcd_up --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_FlavorQCD --var down -o sig_jec_flavorqcd_down --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_HF --var up -o sig_jec_hf_up --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_HF --var down -o sig_jec_hf_down --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_HF_2016preVFP --var up -o sig_jec_hf_2016preVFP_up --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_HF_2016preVFP --var down -o sig_jec_hf_2016preVFP_down --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_HF_2016postVFP --var up -o sig_jec_hf_2016postVFP_up --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_HF_2016postVFP --var down -o sig_jec_hf_2016postVFP_down --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_HF_2017 --var up -o sig_jec_hf_2017_up --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_HF_2017 --var down -o sig_jec_hf_2017_down --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_HF_2018 --var up -o sig_jec_hf_2018_up --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_HF_2018 --var down -o sig_jec_hf_2018_down --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_RelativeBal --var up -o sig_jec_relativebal_up --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_RelativeBal --var down -o sig_jec_relativebal_down --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_RelativeSample_2016preVFP --var up -o sig_jec_relativesample_2016preVFP_up --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_RelativeSample_2016preVFP --var down -o sig_jec_relativesample_2016preVFP_down --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_RelativeSample_2016postVFP --var up -o sig_jec_relativesample_2016postVFP_up --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_RelativeSample_2016postVFP --var down -o sig_jec_relativesample_2016postVFP_down --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_RelativeSample_2017 --var up -o sig_jec_relativesample_2017_up --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_RelativeSample_2017 --var down -o sig_jec_relativesample_2017_down --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_RelativeSample_2018 --var up -o sig_jec_relativesample_2018_up --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type Regrouped_RelativeSample_2018 --var down -o sig_jec_relativesample_2018_down --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type RelativeSample --var up -o sig_jec_relative_up --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jec --jec_type RelativeSample --var down -o sig_jec_relative_down --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jer --var up -o sig_jer_up --outdir ${output_dir} --ana ${analyzer} &
    ./bin/runAnalysis -n 2 -i ${config_file_sig} --jer --var down -o sig_jer_down --outdir ${output_dir} --ana ${analyzer} &
    # ./bin/runAnalysis -n 2 -i ${config_file_sig} --jms --var up -o sig_jms_up --outdir ${output_dir} --ana ${analyzer} &
    # ./bin/runAnalysis -n 2 -i ${config_file_sig} --jms --var down -o sig_jms_down --outdir ${output_dir} --ana ${analyzer} &
    # ./bin/runAnalysis -n 2 -i ${config_file_sig} --jmr --var up -o sig_jmr_up --outdir ${output_dir} --ana ${analyzer} &
    # ./bin/runAnalysis -n 2 -i ${config_file_sig} --jmr --var down -o sig_jmr_down --outdir ${output_dir} --ana ${analyzer} & 

    # Wait for the background process to finish
    wait

}

# Call the function to run the preselection
runPreselection
