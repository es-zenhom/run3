import os
import sys
import uproot
import numpy as np
from glob import glob
import json
import argparse


skims_base_dir = f"/ceph/users/eslam.zenhom/skimdata_test_run3_framework/skimdata"


def extract_sample_year(sample):
    if "UL16" in sample and "APV" in sample:
        return "2016preVFP"
    elif "UL16" in sample and not "APV" in sample:
        return "2016postVFP"
    elif "UL17" in sample or "UL2017" in sample or "Run2017" in sample:
        return "2017"
    elif "UL18" in sample or "UL2018" in sample or "Run2018" in sample:
        return "2018"
    elif ("Run2016B" in sample or "Run2016C" in sample or "Run2016D" in sample or "Run2016E" in sample or "Run2016F" in sample) and "HIPM" in sample:
        return "2016preVFP"
    elif ("Run2016F" in sample or "Run2016G" in sample or "Run2016H" in sample) and not "HIPM" in sample:
        return "2016postVFP"
    else:
        print("Error: year not found for " + sample)
        return None

def yearInName(year, name):
    inName = False
    if year=="2016preVFP":
        inName = ("UL16" in name and "APV" in name) | (("Run2016B" in name or "Run2016C" in name or "Run2016D" in name or "Run2016E" in name or "Run2016F" in name) and "HIPM" in name)
    elif year=="2016postVFP":
        inName = ("UL16" in name and not "APV" in name) | (("Run2016F" in name or "Run2016G" in name or "Run2016H" in name) and not "HIPM" in name)
    elif year=="2017":
        inName =  "UL17" in name or "UL2017" in name
    elif year=="2018":
        inName = "UL18" in name or "UL2018" in name
    else:
        sys.exit(f"Error: year not found {year}")
    return inName

def get_genEventSumw(sample_dir_path, file_regex="output"):
    file_list = glob(sample_dir_path+"/*"+file_regex+"*.root")
    if (len(file_list) == 0):
        sys.exit("Error: len(file_list)==0. Check input files regex?")
    genSumw = 0
    for file_name in file_list:
        with uproot.open(file_name) as f:
            print(f"file_name: {file_name}")
            arr_genEventSumw = f["Runs"].arrays()["genEventSumw"]
            genSumw += np.sum(arr_genEventSumw)
    if genSumw==0:
        print("In get_genEventSumw()")
        print(f"sample_dir_path: {sample_dir_path}")
        print(f"file_regex: {file_regex}")
        print(f"file_list: {file_list}")
        sys.exit(f"Error: no genSumw?")
    else: print(f" --> Retrieved genSumw: {genSumw}")
    return np.sum(genSumw)

def get_xsec(sample_name):
    xsecs_file = os.environ.get('ANALYSISPATH', './')+"/data/xsecs.json"
    with open(xsecs_file, 'r') as j:
        xsecs = json.loads(j.read())
    xsec = [val for key, val in xsecs.items() if key in sample_name]
    if len(xsec) != 1: sys.exit(f"Error: no (or too many) xsec for {sample_name}.")
    else: print(f" --> Retrieved xsec={xsec[0]} for {sample_name}")
    return xsec[0]

def extract_mc_sample_type(sample_name):
    if "DY" in sample_name:
        return "DY"
    elif sample_name.startswith("TTTo"):
        return "ttbar"
    elif "TT" in sample_name or "tt" in sample_name:
        return "ttx"
    elif "ST" in sample_name:
        return "ST"
    elif "WJets" in sample_name:
        return "WJets"
    elif "ZJets" in sample_name:
        return "ZJets"
    elif "EWK" in sample_name:
        return "EWK"
    elif "QCD" in sample_name:
        return "QCD"
    else:
        return "Other"

def get_lumi(year):
    lumi = {
            "2016preVFP": 19.52,
            "2016postVFP": 16.81,
            "2017": 41.529,
            "2018": 59.7
        }
    return lumi[year]

def make_config(args):
    config_name = "input"
    if args.output is None:
        if args.sample_year is not None:
            config_name += "_"+args.sample_year
        if len(args.categories.split(","))<3:
            config_name += "_"+(args.categories).replace(",","_")
    else:
        config_name = args.output

    with open(config_name+".json", "w+") as f:
        try:
            config = json.load(f)
        except json.JSONDecodeError:
            config = None
        if config is None:
            config = {"samples": {}}

        # Background
        # -------------------------------------------------
        if "bkg" in args.categories.split(","):
            skims_job_name = "bkg_0lep_2ak4_2ak8_ttH"
            file_regex = "merged" #"output"
            skims_dir = skims_base_dir+"/"+skims_job_name+"/"
            for sample_dir in os.listdir(skims_dir):
                print('sample_dir : ', sample_dir)
                process = os.path.basename(sample_dir).split("_Tune")[0]
                print('process: ', process)
                sample_year = extract_sample_year(sample_dir)
                if args.sample_year is not None and sample_year != args.sample_year:
                    continue
                xsec = get_xsec(process)
                xsec_sumws =  get_genEventSumw(skims_dir+sample_dir, file_regex)
                config["samples"].update(
                    {
                        process + "_" + sample_year: {
                            "trees": ["Events"],
                            "files": [skims_dir+sample_dir+"/"+file_regex+"*.root"],
                            "metadata": {
                                "sample_category": "bkg",
                                "sample_year" : sample_year,
                                "sample_type": extract_mc_sample_type(process),
                                "xsec": float(xsec),
                                "lumi": get_lumi(sample_year),
                                "sumws": float(xsec_sumws)
                            }
                        }
                    }
                )
        # Data
        # -------------------------------------------------
        if "data" in args.categories.split(","):
            skims_job_name = "data_0lep_2ak4_2ak8_ttH"
            file_regex = "merged"
            skims_dir = skims_base_dir+"/"+skims_job_name+"/"
            for sample_dir in os.listdir(skims_dir):
                print('sample_dir : ', sample_dir)
                sample_name = "Data_"+sample_dir.split("_MiniAOD")[0]
                sample_year = extract_sample_year(sample_dir)
                if args.sample_year is not None and sample_year != args.sample_year:
                    print("Year mismatch")
                    continue

                data_sample_type = "JetHT" if "JetHT" in sample_dir else ""
                config["samples"].update(
                    {
                        sample_name + "_" + sample_year: {
                            "trees": ["Events"],
                            "files": [skims_dir+sample_dir + "/merged*.root"],
                            "metadata": {
                                "sample_category": "data",
                                "sample_year" : sample_year,
                                "sample_type": data_sample_type,
                                "xsec": 1.0,
                                "lumi": 1.0,
                                "sumws": 1.0
                            }
                        }
                    }
                )
        # Signal
        # -------------------------------------------------
        if "sig" in args.categories.split(","):
            skims_job_name = "sig_0lep_2ak4_2ak8_ttH"
            file_regex = "merged"
            skims_dir = skims_base_dir+"/"+skims_job_name+"/"
            for sample_dir in os.listdir(skims_dir):
                print(f"sample_dir : {sample_dir}")
                if not os.path.isdir(skims_dir+sample_dir):
                    continue
                # FIXME: central and c2wc2z skims produced by jguiang are stored in same directory
                # if "VBSVVHSkim" in skims_base_dir:
                if "Private_" in sample_dir:
                    print(f"Skipping private/C2W_C2Z sample: {sample_dir}")
                    continue
                process = (os.path.basename(sample_dir).split("_Tune")[0]).replace("OSWW","OS").replace("SSWW","SS")
                print('process: ', process)
                sample_year = extract_sample_year(sample_dir)
                if args.sample_year is not None and sample_year != args.sample_year:
                    continue
                xsec = get_xsec(process.replace("_MJJ-100_4f","_VBSCuts"))
                xsec_sumws =  get_genEventSumw(skims_dir+sample_dir, file_regex)
                config["samples"].update(
                    {
                        process + "_" + sample_year: {
                            "trees": ["Events"],
                            "files": [skims_dir+sample_dir+"/"+file_regex+"*.root"],
                            "metadata": {
                                "sample_category": "sig",
                                "sample_year" : sample_year,
                                "sample_type": process,
                                "xsec": xsec,
                                "lumi": get_lumi(sample_year),
                                "sumws": xsec_sumws
                            }
                        }
                    }
                )
        json.dump(config, f, indent=4)

if __name__ == "__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--categories", help="categories", default="bkg,data,sig")
    argparser.add_argument("--sample_year", help="sample year", default=None)
    argparser.add_argument("--output", help="output config file without extension", default=None)

    args = argparser.parse_args()

    config = make_config(args)
