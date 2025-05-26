#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

#include <vector>
#include <iostream>

using namespace std;

int main(int argc, char** argv) {
    TFile f(argv[1], "update");
    
    bool new_v;
    auto events = f.Get<TTree>("Events");
    vector<const char *> branches = {"HLT_PFHT1050", "HLT_PFHT800", "HLT_PFHT900"};
    vector<TBranch *> newBranches;
    for (auto branch : branches) {
        if (!events->GetBranch(branch)) {
            cout << "Adding branch " << branch << endl;
            newBranches.push_back(events->Branch(branch, &new_v));
        }
    }
    Long64_t nentries = events->GetEntries();
 
    for (Long64_t i = 0; i < nentries; i++) {
        new_v = false;
        for (auto branch : newBranches) {
            branch->Fill();
        }
    }
 
    events->Write("", TObject::kOverwrite);
}