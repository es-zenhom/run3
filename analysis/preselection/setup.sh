cd /cvmfs/cms.cern.ch/el8_amd64_gcc12/cms/cmssw/CMSSW_14_0_5; cmsenv; cd -
export ANALYSISPATH=$PWD
RHELVER=8

# Export some useful CMSSW paths
CMSSW_EXT=$CMSSW_BASE/../../../external
if [[ $RHELVER == 7 ]]; then
    export CORRECTIONLIBDIR=$CMSSW_EXT/py3-correctionlib/2.0.0-0c4f44c8dd5561d8c0660135feeb81f4/lib/python3.9/site-packages/correctionlib
    export BOOSTDIR=$CMSSW_EXT/boost/1.67.0
elif [[ $RHELVER == 8 ]]; then
    export CORRECTIONLIBDIR=$CMSSW_EXT/py3-correctionlib/2.1.0-6711d4d6283f175d9b6fd011bbaad506/lib/python3.9/site-packages/correctionlib
    export BOOSTDIR=$CMSSW_EXT/boost/1.80.0-a1544032d9d65904ac2112b6d35bba55
fi
