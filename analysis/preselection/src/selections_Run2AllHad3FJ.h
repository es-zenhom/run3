#ifndef SELECTIONS_H
#define SELECTIONS_H

#include "ROOT/RDataFrame.hxx"

#include "utils.h"
#include "corrections.h"

using RNode = ROOT::RDF::RNode;

namespace Run2AllHad3FJ
{
  RNode triggerSelections(RNode df);
  RNode bosonsReconstruction(RNode df);
  RNode selectHiggsCandidate(RNode df_);
  RNode selectVCandidates(RNode df_);
  auto fFindV1CandidateIdx(const ROOT::RVecF &WqqScore, size_t &HiggsIdx);
  auto fFindV2CandidateIdx(const ROOT::RVecF &WqqScore, size_t &HiggsIdx, size_t &V1Idx);
  RNode AK4JetsSelections(RNode df);
  RNode VBSJetsReconstruction(RNode df);
  RNode eventObjectsReconstruction(RNode df);
  RNode runABCDNet(RNode df);
  RNode runAnalysis(RNode df);

} // Run2AllHad3FJ

#endif // SELECTIONS_H