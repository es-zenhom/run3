#ifndef RUN2ALLHAD2FJ_H
#define RUN2ALLHAD2FJ_H

#include "ROOT/RDataFrame.hxx"
#include "utils.h"
#include "corrections.h"

using RNode = ROOT::RDF::RNode;

namespace Run2AllHad2FJ {

// -----------------------------------------------------------------------------
// Helper function declarations
// -----------------------------------------------------------------------------

// For picking the V fat jet (among 2 AK8), skipping the Higgs index
int pickVFatJetIdx(const ROOT::VecOps::RVec<bool> &goodAK8Jets, int hidx);

// For selecting V→qq from leftover AK4 (Mjj>50, minimal ΔR)
ROOT::VecOps::RVec<int> findVqqPair(
    const ROOT::VecOps::RVec<bool> &mask,
    const ROOT::VecOps::RVec<float> &pt,
    const ROOT::VecOps::RVec<float> &eta,
    const ROOT::VecOps::RVec<float> &phi,
    const ROOT::VecOps::RVec<float> &mass);

// For picking the 2 VBS jets from leftover AK4 (pT>30, max Mjj)
ROOT::VecOps::RVec<int> findVBSjetsMaxM(
    const ROOT::VecOps::RVec<bool> &mask,
    const ROOT::VecOps::RVec<int>  &vqqPairIdx,
    const ROOT::VecOps::RVec<float> &pt,
    const ROOT::VecOps::RVec<float> &eta,
    const ROOT::VecOps::RVec<float> &phi,
    const ROOT::VecOps::RVec<float> &mass);

// -----------------------------------------------------------------------------
// Individual selection steps
// -----------------------------------------------------------------------------

// Trigger selection
RNode triggerSelections(RNode df);

// Require exactly 2 fat jets
RNode requireExactly2FatJets(RNode df);

// Select Higgs candidate by highest Xbb score
RNode selectHiggsCandidate(RNode df);

// Select V candidate (the other AK8)
RNode selectVCandidate(RNode df);

// Remove overlap from AK4 jets
RNode removeOverlapAK4(RNode df);

// Require ≥ 4 leftover AK4
RNode requireAtLeast4AK4(RNode df);

// Select V→qq from leftover AK4
RNode selectVqqFromAK4(RNode df);

// Select the 2 VBS jets from leftover (pT>30, max Mjj)
RNode selectVBSJets(RNode df);

// -----------------------------------------------------------------------------
// Main analysis
// -----------------------------------------------------------------------------
RNode runAnalysis(RNode df);

} // namespace Run2AllHad2FJ

#endif // RUN2ALLHAD2FJ_H
