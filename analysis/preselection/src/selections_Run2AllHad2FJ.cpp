#include "selections_Run2AllHad2FJ.h"
#include "commonSelections.h"
#include "TLorentzVector.h"
#include <iostream>

#include <ROOT/RVec.hxx>
#include <Math/VectorUtil.h> 

namespace Run2AllHad2FJ {

// -----------------------------------------------------------------------------
// 1) InvariantMassF for two single 4-vectors, used in vqq pair nad vbs qq 
// -----------------------------------------------------------------------------
static inline float InvariantMassF(
    float pt1,  float eta1, float phi1, float m1,
    float pt2,  float eta2, float phi2, float m2)
{
    // Use ROOT::VecOps::InvariantMass for better performance
    return ROOT::VecOps::InvariantMass({pt1, pt2}, {eta1, eta2}, {phi1, phi2}, {m1, m2});
}

// -----------------------------------------------------------------------------
// 2) pickVFatJetIdx(...)  (Returns int, uses RVec<bool> for goodAK8Jets)
// in the anaylsis we get one Higgs and one V from 2ak8, the way it's selcted is first we select
// the higgs candidate by highest Xbb score, then we select the V from the ak8 as the remaiang one  and we require events to ahve 2k8 jets (semi merged)
// -----------------------------------------------------------------------------
int pickVFatJetIdx(const ROOT::VecOps::RVec<bool> &goodAK8Jets, int hidx)
{
    // Suppose we want exactly 2 "true" jets, the other is V
    ROOT::VecOps::RVec<int> idxs;
    for (size_t i = 0; i < goodAK8Jets.size(); i++) {
        if (goodAK8Jets[i]) idxs.push_back(i);
    }
    if (idxs.size() != 2) return -1;
    if (hidx < 0 || (size_t)hidx >= goodAK8Jets.size()) return -1;

    // the "other" index is whichever is not hidx
    for (auto idx : idxs) {
        if (idx != hidx) return idx;
    }
    return -1;
}

// -----------------------------------------------------------------------------
// 3) findVqqPair(...)  (Returns RVec<int>, uses RVec<bool> for leftover mask)
// we have 4 or more ak4 jets, we first choose the pair for the vqq boson by getting the pair with the smallest delta r 
// thre is a flag for the method proosed in run2 revew if we need to remove pairs with inv mass of less than 50
// -----------------------------------------------------------------------------
ROOT::VecOps::RVec<int> findVqqPair(
    const ROOT::VecOps::RVec<bool>  &mask,
    const ROOT::VecOps::RVec<float> &pt,
    const ROOT::VecOps::RVec<float> &eta,
    const ROOT::VecOps::RVec<float> &phi,
    const ROOT::VecOps::RVec<float> &mass)
{
    // Flag to control whether to apply Mjj > 50 GeV cut on V->qq jets (suggestion in l2 or l3 review but didn't go with it since improvment was small)
    const bool APPLY_VQQ_MASS_CUT = false; // Set to true to require Mjj > 50 GeV

    ROOT::VecOps::RVec<int> idxs;
    for (size_t i = 0; i < mask.size(); i++) {
        if (mask[i]) idxs.push_back(i);
    }
    if (idxs.size() < 2) return {};

    float bestDR = 999.f;
    int   bestA  = -1, bestB = -1;

    // require minimal ΔR and Mjj > 50
    for (size_t a = 0; a < idxs.size(); a++) {
        for (size_t b = a+1; b < idxs.size(); b++) {
            int iA = idxs[a];
            int iB = idxs[b];

            // Calculate invariant mass using ROOT::VecOps
            float m_inv = ROOT::VecOps::InvariantMass({pt[iA], pt[iB]}, 
                                                      {eta[iA], eta[iB]}, 
                                                      {phi[iA], phi[iB]}, 
                                                      {mass[iA], mass[iB]});
            
            if (APPLY_VQQ_MASS_CUT && m_inv < 50) continue; // Skip pairs with Mjj < 50 GeV only if flag is true

            float dR = ROOT::VecOps::DeltaR(eta[iA], eta[iB], phi[iA], phi[iB]);
            if (dR < bestDR) {
                bestDR = dR;
                bestA  = iA;
                bestB  = iB;
            }
        }
    }
    if (bestA < 0) return {};

    // sort by pT
    if (pt[bestA] >= pt[bestB]) return {bestA, bestB};
    else                        return {bestB, bestA};
}

// -----------------------------------------------------------------------------
// 4) findVBSjetsMaxDeltaEta(...)  (Returns RVec<int>, skip V->qq jets, pick 2 leftover w/ max delta eta)
// we have 4 or more ak4 jets, we first remove the vqq pair, then we select the 2 VBS jets from the remaining ak4 jets by getting the pair with the maximum delta eta
// -----------------------------------------------------------------------------
ROOT::VecOps::RVec<int> findVBSjetsMaxDeltaEta(
    const ROOT::VecOps::RVec<bool> &mask,
    const ROOT::VecOps::RVec<int>  &vqqPairIdx,
    const ROOT::VecOps::RVec<float>&pt,
    const ROOT::VecOps::RVec<float>&eta,
    const ROOT::VecOps::RVec<float>&phi,
    const ROOT::VecOps::RVec<float>&mass)
{
    ROOT::VecOps::RVec<int> idxs;
    for (size_t i=0; i<mask.size(); i++){
        if (!mask[i]) continue;
        
        if (pt[i] < 30) continue;
        
        if (fabs(eta[i]) >= 4.7) continue;

        bool usedInVqq = false;
        for (auto j : vqqPairIdx) {
            if ((int)i == j) { usedInVqq = true; break; }
        }
        if (usedInVqq) {
            continue;
        }
        idxs.push_back(i);
    }
    
    if (idxs.size() < 2) return {};

    float bestDeltaEta = -999.f;
    int   bestA = -1, bestB = -1;

    // Select pair with maximum delta eta
    for (size_t a=0; a<idxs.size(); a++) {
        for (size_t b=a+1; b<idxs.size(); b++) {
            int iA = idxs[a];
            int iB = idxs[b];
            float deltaEta = fabs(eta[iA] - eta[iB]);
            if (deltaEta > bestDeltaEta) {
                bestDeltaEta = deltaEta;
                bestA = iA;
                bestB = iB;
            }
        }
    }
    if (bestA < 0) return {};

    // sort by pT
    if (pt[bestA] >= pt[bestB]) return {bestA, bestB};
    else                        return {bestB, bestA};
}

// ============================================================================
// Implementation of each selection step
// ============================================================================

// -----------------------------------------------------------------------------
// triggerSelections
// -----------------------------------------------------------------------------
ROOT::RDF::RNode triggerSelections(ROOT::RDF::RNode df)
{
    auto df_out = df.Define("Pass_Triggers", 
        [](bool is2016, bool is2017, bool is2018, bool isData,
           bool HLT_PFHT800, bool HLT_PFHT900, bool HLT_PFHT1050) {
            bool passed = false;
            if (is2016) {
                passed = HLT_PFHT800 || HLT_PFHT900;
            } else if (is2017 || is2018) {
                passed = HLT_PFHT1050;
            }
            return passed;
        },
        {"is2016", "is2017", "is2018", "isData", 
         "HLT_PFHT800", "HLT_PFHT900", "HLT_PFHT1050"});
    return df_out;
}

// -----------------------------------------------------------------------------
// requireExactly2FatJets
// goodak8jets from the RNode AK8JetsPreselection(RNode df) in common sleections
// it has fatjet_pt>300,  abs(FatJet_eta) < 2.5, CorrFatJet_mass > 50, FatJet_msoftdrop > 40, FatJet_jetId > 0
// look https://github.com/jkguiang/vbs/blob/3269f28e41e7b3e2113ac9756db5aa06db48d4fe/analysis/include/core/cuts.h#L487-L491
// -----------------------------------------------------------------------------
ROOT::RDF::RNode requireExactly2FatJets(ROOT::RDF::RNode df)
{
    // goodAK8Jets is RVec<bool>, "Sum(goodAK8Jets)" is the # of true jets
    auto df_out = df.Define("Pass_Exactly2FatJets", "Sum(goodAK8Jets) == 2");
    return df_out;
}

// -----------------------------------------------------------------------------
// selectHiggsCandidate
// -----------------------------------------------------------------------------
ROOT::RDF::RNode selectHiggsCandidate(ROOT::RDF::RNode df)
{
    auto df_out = df
      .Define("HiggsScoreIdx_ulong",
              "(goodAK8Jets_HbbScore.size() > 0) ? ArgMax(goodAK8Jets_HbbScore) : 999UL")
      .Define("HiggsScoreIdx",
              "static_cast<int>(HiggsScoreIdx_ulong)")  
      .Define("HiggsScore",
              "(HiggsScoreIdx != 999) ? goodAK8Jets_HbbScore[HiggsScoreIdx] : -999")
      .Define("Pass_HiggsCandidateExists",
              "(HiggsScoreIdx != 999) && (HiggsScore > 0)")
      .Define("Higgs_pt",
              "Pass_HiggsCandidateExists ? goodAK8Jets_pt[HiggsScoreIdx] : -999")
      .Define("Higgs_eta",
              "Pass_HiggsCandidateExists ? goodAK8Jets_eta[HiggsScoreIdx]     : -999")
      .Define("Higgs_phi",
              "Pass_HiggsCandidateExists ? goodAK8Jets_phi[HiggsScoreIdx]     : -999")
      .Define("Higgs_mass",
              "Pass_HiggsCandidateExists ? goodAK8Jets_mass[HiggsScoreIdx] : -999")
      .Define("Higgs_msoftdrop",
              "Pass_HiggsCandidateExists ? goodAK8Jets_msoftdrop[HiggsScoreIdx] : -999")
      .Define("Higgs_xbb",
              "Pass_HiggsCandidateExists ? goodAK8Jets_HbbScore[HiggsScoreIdx] : -999");
    return df_out;
}
// -----------------------------------------------------------------------------
// selectVCandidate
//   Instead of .Define<int>("VfatJetIdx", pickVFatJetIdx, {...}),
//  pickVFatJetIdx(...) is a lambda so that we get it exactly "int". note: the pickVFatJetIdx function
// gives the index of V since we choose the higgs first in th eorder df5 and df6, so we have the h index and we choose the other index of 1k8 to be the higgs 
// -----------------------------------------------------------------------------
ROOT::RDF::RNode selectVCandidate(ROOT::RDF::RNode df)
{
    auto df_out = df
      .Define("VfatJetIdx",
        [](const ROOT::VecOps::RVec<bool> &goodAK8Jets, int hidx){
            return pickVFatJetIdx(goodAK8Jets, hidx);
        },
        {"goodAK8Jets", "HiggsScoreIdx"}
      )
      .Define("Pass_VfatCandidateExists", "VfatJetIdx >= 0")
      .Define("Vfat_pt",
              "Pass_VfatCandidateExists ? goodAK8Jets_pt[VfatJetIdx] : -999")
      .Define("Vfat_eta",
              "Pass_VfatCandidateExists ? goodAK8Jets_eta[VfatJetIdx] : -999")
      .Define("Vfat_phi",
              "Pass_VfatCandidateExists ? goodAK8Jets_phi[VfatJetIdx] : -999")
      .Define("Vfat_mass",
              "Pass_VfatCandidateExists ? goodAK8Jets_mass[VfatJetIdx] : -999")
      .Define("Vfat_msoftdrop",
              "Pass_VfatCandidateExists ? goodAK8Jets_msoftdrop[VfatJetIdx] : -999")
      .Define("Vfat_xqq",
              "Pass_VfatCandidateExists ? goodAK8Jets_WqqScore[VfatJetIdx] : -999");

    return df_out;
}
// -----------------------------------------------------------------------------
// removeOverlapAK4
// -----------------------------------------------------------------------------
ROOT::RDF::RNode removeOverlapAK4(ROOT::RDF::RNode df)
{
    auto df_out = df
      .Define("AK4HDeltaR", VfDeltaR, {"Jet_eta", "Jet_phi", "Higgs_eta", "Higgs_phi"})
      .Define("AK4VDeltaR", VfDeltaR, {"Jet_eta", "Jet_phi", "Vfat_eta", "Vfat_phi"})
      // Convert goodAK4Jets from int to bool for consistency
      .Define("goodAK4Jets_bool", "ROOT::VecOps::Map(goodAK4Jets, [](int x){ return x != 0; })")
      .Redefine("goodAK4Jets_bool", "goodAK4Jets_bool && "
                                     "AK4HDeltaR >= 0.8 && "
                                     "AK4VDeltaR >= 0.8")
      // Convert back to int for compatibility
      .Define("goodAK4Jets_noFatOverlap", "ROOT::VecOps::Map(goodAK4Jets_bool, [](bool x){ return x ? 1 : 0; })");
    // HEM veto check 
    df_out = df_out.Define("passHEMVeto",
      [](bool is2018, bool isData, unsigned int run, unsigned long long event,
         const ROOT::VecOps::RVec<float> &pt,
         const ROOT::VecOps::RVec<float> &eta,
         const ROOT::VecOps::RVec<float> &phi,
         const ROOT::VecOps::RVec<int> &jetId,
         const ROOT::VecOps::RVec<int> &puId)
      {
        // Only apply HEM veto for 2018
        if (!is2018) return true;
        
        // Check if this event should have HEM veto applied
        bool applyHEM = (isData && run >= 319077) || (!isData && event % 1961 < 1286);
        if (!applyHEM) return true;
        
        // Check ALL jets, not just good jets
        for (size_t i = 0; i < pt.size(); i++) {
          // Check if jet is in HEM-affected region with pt > 15
          if (pt[i] > 15 && 
              phi[i] > -1.57 && phi[i] < -0.87 &&
              eta[i] > -3.20 && eta[i] < -1.30) {
            
            // Check jet ID requirements (same as in Core::SelectJets for run2)
            bool passesIDs = true;
            if (jetId[i] < 2) passesIDs = false;  // 2017/2018 require jetId >= 2
            if (pt[i] < 50 && puId[i] == 0) passesIDs = false;
            
            // Veto event if this jet passes the IDs
            if (passesIDs) return false;
          }
        }
        return true;
      },
      {"is2018", "isData", "run", "event",
       "CorrJet_pt", "Jet_eta", "Jet_phi", "Jet_jetId", "Jet_puId"}
    );

    return df_out;
}

// -----------------------------------------------------------------------------
// requireAtLeast4AK4, in run2 that step had the jec and jer applied
// -----------------------------------------------------------------------------
ROOT::RDF::RNode requireAtLeast4AK4(ROOT::RDF::RNode df)
{
    // n_goodAK4Jets_noFatOverlap is already defined in hemVeto function
    auto df_out = df.Define("Pass_AtLeast4AK4Jets", "n_goodAK4Jets_noFatOverlap >= 4");
    
    return df_out;
}

// -----------------------------------------------------------------------------
// selectVqqFromAK4
//   Instead of .Define< RVec<int> >(... findVqqPair, {...})
//   we wrap findVqqPair in a lambda
// -----------------------------------------------------------------------------
ROOT::RDF::RNode selectVqqFromAK4(ROOT::RDF::RNode df)
{
    auto df_out = df
      .Define("VqqPairIdx",
        [](const ROOT::VecOps::RVec<bool>  &mask,
           const ROOT::VecOps::RVec<float> &pt,
           const ROOT::VecOps::RVec<float> &eta,
           const ROOT::VecOps::RVec<float> &phi,
           const ROOT::VecOps::RVec<float> &mass){
            return findVqqPair(mask, pt, eta, phi, mass);
        },
        {"goodAK4Jets_noFatOverlap", "CorrJet_pt", "Jet_eta", "Jet_phi", "CorrJet_mass"}
      )
      .Define("n_VqqPair", "static_cast<int>(VqqPairIdx.size())")
      .Define("Pass_VqqPair", "n_VqqPair == 2")
      .Define("Vqq_m",
        [](bool pass,
           const ROOT::VecOps::RVec<int> &idxs,
           const ROOT::VecOps::RVec<float> &pt,
           const ROOT::VecOps::RVec<float> &eta,
           const ROOT::VecOps::RVec<float> &phi,
           const ROOT::VecOps::RVec<float> &mass)
        {
            if (!pass) return -1.f;
            return InvariantMassF(
                pt[idxs[0]], eta[idxs[0]], phi[idxs[0]], mass[idxs[0]],
                pt[idxs[1]], eta[idxs[1]], phi[idxs[1]], mass[idxs[1]]
            );
        },
        {"Pass_VqqPair","VqqPairIdx","CorrJet_pt","Jet_eta","Jet_phi","CorrJet_mass"}
      )
      .Define("Vqq_dR",
        [](bool pass,
           const ROOT::VecOps::RVec<int> &idxs,
           const ROOT::VecOps::RVec<float> &eta,
           const ROOT::VecOps::RVec<float> &phi)
        {
            if (!pass) return -1.f;
            return ROOT::VecOps::DeltaR(
                eta[idxs[0]], eta[idxs[1]],
                phi[idxs[0]], phi[idxs[1]]
            );
        },
        {"Pass_VqqPair","VqqPairIdx","Jet_eta","Jet_phi"}
      );

    return df_out;
}

// -----------------------------------------------------------------------------
// selectVBSJets
//   Instead of .Define< RVec<int> >(... findVBSjetsMaxDeltaEta, {...})
//   we wrap findVBSjetsMaxDeltaEta in a lambda
// -----------------------------------------------------------------------------
ROOT::RDF::RNode selectVBSJets(ROOT::RDF::RNode df)
{
    auto df_out = df
      .Define("VBSjetsIdx",
        [](const ROOT::VecOps::RVec<bool> &mask,
           const ROOT::VecOps::RVec<int>  &vqqPairIdx,
           const ROOT::VecOps::RVec<float>&pt,
           const ROOT::VecOps::RVec<float>&eta,
           const ROOT::VecOps::RVec<float>&phi,
           const ROOT::VecOps::RVec<float>&mass)
        {
            return findVBSjetsMaxDeltaEta(mask, vqqPairIdx, pt, eta, phi, mass);
        },
        {"goodAK4Jets_noFatOverlap", "VqqPairIdx",
         "CorrJet_pt", "Jet_eta", "Jet_phi", "CorrJet_mass"}
      )
      .Define("n_VBSjets", "static_cast<int>(VBSjetsIdx.size())")
      .Define("Pass_VBSjets", "n_VBSjets == 2")
      // define vbs1 / vbs2 kinematics
      .Define("vbs1_pt",  "Pass_VBSjets ? CorrJet_pt[VBSjetsIdx[0]] : -999")
      .Define("vbs1_eta", "Pass_VBSjets ? Jet_eta[VBSjetsIdx[0]]     : -999")
      .Define("vbs1_phi", "Pass_VBSjets ? Jet_phi[VBSjetsIdx[0]]     : -999")
      .Define("vbs1_m",   "Pass_VBSjets ? CorrJet_mass[VBSjetsIdx[0]]: -999")
      .Define("vbs2_pt",  "Pass_VBSjets ? CorrJet_pt[VBSjetsIdx[1]]  : -999")
      .Define("vbs2_eta", "Pass_VBSjets ? Jet_eta[VBSjetsIdx[1]]     : -999")
      .Define("vbs2_phi", "Pass_VBSjets ? Jet_phi[VBSjetsIdx[1]]     : -999")
      .Define("vbs2_m",   "Pass_VBSjets ? CorrJet_mass[VBSjetsIdx[1]]: -999")
      // define Mjj, detajj
      .Define("VBS_Mjj",
        [](bool pass,
           float pt1, float eta1, float phi1, float m1,
           float pt2, float eta2, float phi2, float m2)
        {
            if (!pass) return -1.f;
            return InvariantMassF(pt1, eta1, phi1, m1,
                                  pt2, eta2, phi2, m2);
        },
        {"Pass_VBSjets",
         "vbs1_pt","vbs1_eta","vbs1_phi","vbs1_m",
         "vbs2_pt","vbs2_eta","vbs2_phi","vbs2_m"}
      )
      .Define("VBS_detajj",
              "Pass_VBSjets ? fabs(vbs1_eta - vbs2_eta) : -1");

    return df_out;
}

// -----------------------------------------------------------------------------
// OPTIONAL: eventObjectsReconstruction
// -----------------------------------------------------------------------------
ROOT::RDF::RNode eventObjectsReconstruction(ROOT::RDF::RNode df)
{
    std::cout << " -> Run2AllHad2FJ::eventObjectsReconstruction()" << std::endl;

    auto df_out = df
      .Define("puIDJets_pt",
              "CorrJet_pt[ goodAK4Jets_noFatOverlap ]")
      .Define("puIDJets_eta",
              "Jet_eta[ goodAK4Jets_noFatOverlap ]")
      .Define("MET",     "CorrMET_pt")
      .Define("ht_ak4",  "Sum(CorrJet_pt[goodAK4Jets_noFatOverlap])")
      .Define("ht_ak8",  "Sum(goodAK8Jets_pt)")
      .Define("ST",      "Higgs_pt + Vfat_pt + CorrMET_pt");

    return df_out;
}

// -----------------------------------------------------------------------------
// initialJetCut: require at least 2 jets with pt > 30 and max fat jet pt > 550
// -----------------------------------------------------------------------------
ROOT::RDF::RNode initialJetCut(ROOT::RDF::RNode df)
{
    std::cout << " -> Run2AllHad2FJ::initialJetCut()" << std::endl;
    
    auto df_out = df
        .Define("nJets_pt30", "(int) CorrJet_pt[CorrJet_pt > 30].size()")
        .Define("maxFatJetPt", "CorrFatJet_pt.size() > 0 ? Max(CorrFatJet_pt) : 0.0")
        .Define("Pass_InitialJetCut", "nJets_pt30 >= 2 && maxFatJetPt > 550");
    
    return df_out;
}

// -----------------------------------------------------------------------------
// triggerPlateauCut: require leading good fat jet pt > 550
// -----------------------------------------------------------------------------
ROOT::RDF::RNode triggerPlateauCut(ROOT::RDF::RNode df)
{
    std::cout << " -> Run2AllHad2FJ::triggerPlateauCut()" << std::endl;
    
    auto df_out = df
        .Define("leadingGoodFatJetPt", 
                "goodAK8Jets_pt.size() > 0 ? goodAK8Jets_pt[0] : 0.0")
        .Define("Pass_TriggerPlateau", "leadingGoodFatJetPt > 550");
    
    return df_out;
}

// -----------------------------------------------------------------------------
// runAnalysis
// -----------------------------------------------------------------------------
ROOT::RDF::RNode runAnalysis(ROOT::RDF::RNode df)
{
    std::cout << " -> Running analysis (SemiMerged)" << std::endl;
    // Step 0: common pre-selections
    auto df0 = runCommonSelections(df);

    // Step 1: initial jet cut
    auto df1 = initialJetCut(df0);

    // Step 2: triggers
    auto df2 = triggerSelections(df1);

    // Step 3: trigger plateau cut
    auto df3 = triggerPlateauCut(df2);

    // Step 4: exactly 2 AK8
    auto df4 = requireExactly2FatJets(df3);

    // Step 5: pick Higgs
    auto df5 = selectHiggsCandidate(df4);

    // Step 6: pick V
    auto df6 = selectVCandidate(df5);

    // Step 7: remove overlap
    auto df7 = removeOverlapAK4(df6);

    // Step 8: require >= 4 leftover
    auto df8 = requireAtLeast4AK4(df7);

    // Step 9: V->qq from leftover
    auto df9 = selectVqqFromAK4(df8);

    // Step 10: pick VBS
    auto df10 = selectVBSJets(df9);

    // Optional: define additional columns
    auto df11 = eventObjectsReconstruction(df10);

    // Define final pass cuts, for example
    auto df_out = df11
      .Define("Pass_xwqq",                     "Vfat_xqq > 0.6")  // V fat jet tagging
      .Define("Pass_xbb",                      "Higgs_xbb > 0.8")  // H fat jet tagging
      .Define("passCut1_Semi_EventFilters",    "Pass_EventFilters")
      .Define("passCut2_Semi_InitialJetCut",   "passCut1_Semi_EventFilters && Pass_InitialJetCut")
      .Define("passCut3_Semi_Triggers",        "passCut2_Semi_InitialJetCut && Pass_Triggers")
      .Define("passCut4_Semi_TriggerPlateau",  "passCut3_Semi_Triggers && Pass_TriggerPlateau")
      .Define("passCut5_Semi_Exactly2FatJets", "passCut4_Semi_TriggerPlateau && Pass_Exactly2FatJets")
      .Define("passCut6_Semi_Higgs",           "passCut5_Semi_Exactly2FatJets && Pass_HiggsCandidateExists")
      .Define("passCut7_Semi_Vfat",            "passCut6_Semi_Higgs && Pass_VfatCandidateExists")
      .Define("passCut8_Semi_OverlapRemoval",  "passCut7_Semi_Vfat && passHEMVeto")  // Now includes HEM veto
      .Define("passCut9_Semi_4AK4",            "passCut8_Semi_OverlapRemoval && Pass_AtLeast4AK4Jets")
      .Define("passCut10_Semi_VqqPair",        "passCut9_Semi_4AK4 && Pass_VqqPair")
      .Define("passCut11_Semi_VBS",            "passCut10_Semi_VqqPair && Pass_VBSjets")
      .Define("passCut12_Semi_xwqq",           "passCut11_Semi_VBS && Pass_xwqq")
      .Define("passCut13_Semi_xbb",            "passCut12_Semi_xwqq && Pass_xbb")
      // Backward compatibility aliases, just to ignore changing names in command line while running
      // the anlaysis a thousand time for debugging -_-
      .Define("passCut10",                     "passCut13_Semi_xbb")
      .Define("passCut9",                      "passCut13_Semi_xbb");
    
    auto count_initial = df.Count();
    auto count_cut1 = df_out.Filter("passCut1_Semi_EventFilters").Count();
    auto count_cut2 = df_out.Filter("passCut2_Semi_InitialJetCut").Count();
    auto count_cut3 = df_out.Filter("passCut3_Semi_Triggers").Count();
    auto count_cut4 = df_out.Filter("passCut4_Semi_TriggerPlateau").Count();
    auto count_cut5 = df_out.Filter("passCut5_Semi_Exactly2FatJets").Count();
    auto count_cut6 = df_out.Filter("passCut6_Semi_Higgs").Count();
    auto count_cut7 = df_out.Filter("passCut7_Semi_Vfat").Count();
    auto count_cut8 = df_out.Filter("passCut8_Semi_OverlapRemoval").Count();
    auto count_cut9 = df_out.Filter("passCut9_Semi_4AK4").Count();
    auto count_cut10 = df_out.Filter("passCut10_Semi_VqqPair").Count();
    auto count_cut11 = df_out.Filter("passCut11_Semi_VBS").Count();
    auto count_cut12 = df_out.Filter("passCut12_Semi_xwqq").Count();
    auto count_cut13 = df_out.Filter("passCut13_Semi_xbb").Count();
    
    // debug counts
    auto count_vqq_not_vbs = df_out.Filter("passCut10_Semi_VqqPair && !Pass_VBSjets").Count();

    // Print cutflow
    std::cout << "===== Cutflow Report =====" << std::endl;
    std::cout << "Initial events: " << *count_initial << std::endl;
    std::cout << "After Event Filters: " << *count_cut1 << std::endl;
    std::cout << "After Initial Jet Cut: " << *count_cut2 << std::endl;
    std::cout << "After Triggers: " << *count_cut3 << std::endl;
    std::cout << "After Trigger Plateau: " << *count_cut4 << std::endl;
    std::cout << "After Exactly 2 Fat Jets: " << *count_cut5 << std::endl;
    std::cout << "After Higgs Selection: " << *count_cut6 << std::endl;
    std::cout << "After V Fat Selection: " << *count_cut7 << std::endl;
    std::cout << "After Overlap Removal: " << *count_cut8 << std::endl;
    std::cout << "After 4+ AK4 Jets: " << *count_cut9 << std::endl;
    std::cout << "After V→qq Selection: " << *count_cut10 << std::endl;
    std::cout << "After VBS Jets Selection: " << *count_cut11 << std::endl;
    std::cout << "After xwqq > 0.6: " << *count_cut12 << std::endl;
    std::cout << "After xbb > 0.8: " << *count_cut13 << " (Final)" << std::endl;
    
    std::cout << "\n=== Additional Debug Info ===" << std::endl;
    std::cout << "Events passing V→qq but NOT VBS: " << *count_vqq_not_vbs << std::endl;
    
    return df_out;
}

} // namespace Run2AllHad2FJ
