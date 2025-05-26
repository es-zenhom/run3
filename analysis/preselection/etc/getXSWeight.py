import os 
import uproot
import numpy as np
import json 
import sys
# Debugging script

skims_job_name = "skims_v2_haveLHEbranch"
skims_base_dir = "/data/userdata/mmazza/vvhjj/skims/"+skims_job_name+"/"
sample_name = "VBSWZH"
year = "20UL16"


    
def get_xsec(sample_name):
    xsecs_file = os.environ.get('ANALYSISPATH', './')+"/data/xsecs.json"
    with open(xsecs_file, 'r') as j:
        xsecs = json.loads(j.read())
        print(xsecs)
    xsec = None
    for sample, xsec_tmp in xsecs.items():
        if sample in sample_name:
            xsec = xsec_tmp
            break
    if not xsec: sys.error(f"Error: xsec not found for {sample_name}.")
    else: print(f" --> Retrieved xsec={xsec} for {sample_name}.")
    return xsec

get_xsec("ZJetsToQQ_HT-400to600")

# for subdir in os.listdir(skims_base_dir):
#     if sample_name in subdir and year in subdir:
#         skim_file_name = skims_base_dir+subdir+"/merged.root"
#         with uproot.open(skim_file_name) as f:
#           xsec_weight_sum = f["Runs"].arrays()["genEventSumw"]
#           print(xsec_weight_sum)
#           print(f"{skim_file_name}: {np.sum(xsec_weight_sum)}")
      