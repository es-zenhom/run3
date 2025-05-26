#include "commonSelections.h"
#include "selections_Run2AllHad3FJ.h"
#include "ABCDNet_Run2AllHad3FJ.h"

namespace Run2AllHad3FJ {
    RNode triggerSelections(RNode df)
    {
        std::cout << " -> Run Run2AllHad3FJ::triggerSelections()" << std::endl;
        auto df_triggers = df.Define("Pass_Triggers",
                                    "((is2018 && HLT_PFHT1050 == true) || "
                                    "(is2017 && HLT_PFHT1050 == true) || "
                                    "(is2016 && (HLT_PFHT800 == true || HLT_PFHT900 == true)) )");
        return df_triggers;
    }

    /*
    *   AK8 Jets selection
    */
    RNode bosonsReconstruction(RNode df)
    {
        std::cout << " -> Run Run2AllHad3FJ::bosonsReconstruction()" << std::endl;

        // Select H, V1, and V2 candidates and find their overlap
        auto df_out = selectHiggsCandidate(df);
        df_out = selectVCandidates(df_out);
        // NOTE: in current all-had this cut is applied after the fatjet selection, while 1-lep applies it before
        df_out = df_out.Define("HV1DeltaR", "ROOT::VecOps::DeltaR(Higgs_eta,  V1_eta, Higgs_phi, V1_phi)")
                    .Define("HV2DeltaR", "ROOT::VecOps::DeltaR(Higgs_eta, V2_eta, Higgs_phi, V2_phi)")
                    .Define("V1V2DeltaR", "ROOT::VecOps::DeltaR(V1_eta,  V2_eta, V1_phi, V2_phi)");

        return df_out;
    }

    RNode selectHiggsCandidate(RNode df)
    {
        std::cout << " -> Run Run2AllHad3FJ::selectHiggsCandidate()" << std::endl;
        auto df_out = df.Define("HiggsScoreIdx", "goodAK8Jets_HbbScore.size() != 0 ? ArgMax(goodAK8Jets_HbbScore) : 999")
                          .Define("HiggsScore", "HiggsScoreIdx != 999 ? goodAK8Jets_HbbScore[HiggsScoreIdx] : -999")
                          .Define("Pass_BoostedHiggsCandidateExists", "HiggsScore > 0")
                          .Define("Higgs_pt", "Pass_BoostedHiggsCandidateExists ? goodAK8Jets_pt[HiggsScoreIdx] : -999")
                          .Define("Higgs_eta", "Pass_BoostedHiggsCandidateExists ? goodAK8Jets_eta[HiggsScoreIdx] : -999")
                          .Define("Higgs_phi", "Pass_BoostedHiggsCandidateExists ? goodAK8Jets_phi[HiggsScoreIdx] : -999")
                          .Define("Higgs_mass", "Pass_BoostedHiggsCandidateExists ? goodAK8Jets_mass[HiggsScoreIdx] : -999")
                          .Define("Higgs_msoftdrop", "Pass_BoostedHiggsCandidateExists ? goodAK8Jets_msoftdrop[HiggsScoreIdx] : -999")
                          .Define("Higgs_particleNet_mass", "Pass_BoostedHiggsCandidateExists ? goodAK8Jets_particleNet_mass[HiggsScoreIdx] : -999");

        return df_out;
    }

    auto fFindV1CandidateIdx(const ROOT::RVecF &WqqScore, size_t &HiggsIdx)
    {
        float v1_score{-999};
        size_t v1_idx{999};
        for (size_t i = 0; i < WqqScore.size(); i++)
        {
            if (i == HiggsIdx)
            {
                continue;
            }
            if (WqqScore.at(i) > v1_score)
            {
                v1_score = WqqScore.at(i);
                v1_idx = i;
            }
        }
        return v1_idx;
    }

    auto fFindV2CandidateIdx(const ROOT::RVecF &WqqScore, size_t &HiggsIdx, size_t &V1Idx)
    {
        float v2_score{-999};
        size_t v2_idx{999};
        for (size_t i = 0; i < WqqScore.size(); i++)
        {
            if (i == HiggsIdx)
            {
                continue;
            }
            if (i == V1Idx)
            {
                continue;
            }
            if (WqqScore.at(i) > v2_score)
            {
                v2_score = WqqScore.at(i);
                v2_idx = i;
            }
        }
        return v2_idx;
    }

    RNode selectVCandidates(RNode df)
    {
        std::cout << " -> Run Run2AllHad3FJ::selectVCandidates()" << std::endl;
        // Leading V (V1)
        auto df_out = df.Define("V1ScoreIdx", fFindV1CandidateIdx, {"goodAK8Jets_WqqScore", "HiggsScoreIdx"})
                          .Define("V1Score", "V1ScoreIdx != 999 ? goodAK8Jets_WqqScore[V1ScoreIdx] : -999")
                          .Define("Pass_BoostedV1CandidateExists", "V1Score > 0")
                          .Define("V1_pt", "Pass_BoostedV1CandidateExists ? goodAK8Jets_pt[V1ScoreIdx] : -999")
                          .Define("V1_eta", "Pass_BoostedV1CandidateExists ? goodAK8Jets_eta[V1ScoreIdx] : -999")
                          .Define("V1_phi", "Pass_BoostedV1CandidateExists ? goodAK8Jets_phi[V1ScoreIdx] : -999")
                          .Define("V1_mass", "Pass_BoostedV1CandidateExists ? goodAK8Jets_mass[V1ScoreIdx] : -999")
                          .Define("V1_msoftdrop", "Pass_BoostedV1CandidateExists ? goodAK8Jets_msoftdrop[V1ScoreIdx] : -999")
                          .Define("V1_particleNet_mass", "Pass_BoostedV1CandidateExists ? goodAK8Jets_particleNet_mass[V1ScoreIdx] : -999")
                          // Trailing V (V2)
                          .Define("V2ScoreIdx", fFindV2CandidateIdx, {"goodAK8Jets_WqqScore", "HiggsScoreIdx", "V1ScoreIdx"})
                          .Define("V2Score", "V2ScoreIdx != 999 ? goodAK8Jets_WqqScore[V2ScoreIdx] : -999")
                          .Define("Pass_BoostedV2CandidateExists", "V2Score > 0")
                          .Define("V2_pt", "Pass_BoostedV2CandidateExists ? goodAK8Jets_pt[V2ScoreIdx] : -999")
                          .Define("V2_eta", "Pass_BoostedV2CandidateExists ? goodAK8Jets_eta[V2ScoreIdx] : -999")
                          .Define("V2_phi", "Pass_BoostedV2CandidateExists ? goodAK8Jets_phi[V2ScoreIdx] : -999")
                          .Define("V2_mass", "Pass_BoostedV2CandidateExists ? goodAK8Jets_mass[V2ScoreIdx] : -999")
                          .Define("V2_msoftdrop", "Pass_BoostedV2CandidateExists ? goodAK8Jets_msoftdrop[V2ScoreIdx] : -999")
                          .Define("V2_particleNet_mass", "Pass_BoostedV2CandidateExists ? goodAK8Jets_particleNet_mass[V2ScoreIdx] : -999");
        return df_out;
    }

    /*
    *   AK4 and VBS Jets selection
    */

    RNode AK4JetsSelections(RNode df)
    {
        std::cout << " -> Run Run2AllHad3FJ::AK4JetsSelections()" << std::endl;
        auto df_out = df.Define("AK4HDeltaR", VfDeltaR, {"Jet_eta", "Jet_phi", "Higgs_eta", "Higgs_phi"})
                        .Define("AK4V1DeltaR", VfDeltaR, {"Jet_eta", "Jet_phi", "V1_eta", "V1_phi"})
                        .Define("AK4V2DeltaR", VfDeltaR, {"Jet_eta", "Jet_phi", "V2_eta", "V2_phi"})
                        .Redefine("goodAK4Jets", "goodAK4Jets && "
                                                "AK4HDeltaR >= 0.8 && "
                                                "AK4V1DeltaR >= 0.8 && "
                                                "AK4V2DeltaR >= 0.8")
                        .Define("goodAK4FromTightBJet", "goodAK4Jets && (Jet_btagDeepFlavB > ak4tightBjetScore)");

        return df_out;
    }

    RNode VBSJetsReconstruction(RNode df)
    {
        std::cout << " -> Run Run2AllHad3FJ::VBSJetsSelections()" << std::endl;
        auto df_out = df.Redefine("goodVBSJets", "goodVBSJets && "
                                                 "CorrJet_pt >= 30 && " // ! tighter pT cut than pre-selection
                                                 "AK4HDeltaR >= 0.8 && "
                                                 "AK4V1DeltaR >= 0.8 && "
                                                 "AK4V2DeltaR >= 0.8")
                          .Redefine("goodVBSJets_pt", "CorrJet_pt[goodVBSJets]")
                          .Redefine("goodVBSJets_eta", "Jet_eta[goodVBSJets]")
                          .Redefine("goodVBSJets_phi", "Jet_phi[goodVBSJets]")
                          .Redefine("goodVBSJets_mass", "CorrJet_mass[goodVBSJets]");

        // Reconstruct VBS jets using pair with largest momenta in forward/backward regions
        // .Define("VBSjetidxs", VBS_MaxE, {"goodVBSJets_pt", "goodVBSJets_eta", "goodVBSJets_phi", "goodVBSJets_mass"})
        // .Define("VBSjet1pt", "goodVBSJets_pt[VBSjetidxs[0]]")
        // .Define("VBSjet1eta", "goodVBSJets_eta[VBSjetidxs[0]]")
        // .Define("VBSjet1phi", "goodVBSJets_phi[VBSjetidxs[0]]")
        // .Define("VBSjet1mass", "goodVBSJets_mass[VBSjetidxs[0]]")
        // .Define("VBSjet2pt", "goodVBSJets_pt[VBSjetidxs[1]]")
        // .Define("VBSjet2eta", "goodVBSJets_eta[VBSjetidxs[1]]")
        // .Define("VBSjet2phi", "goodVBSJets_phi[VBSjetidxs[1]]")
        // .Define("VBSjet2mass", "goodVBSJets_mass[VBSjetidxs[1]]")
        // .Define("VBSptjj", "VBSjet1pt + VBSjet2pt")
        // .Define("VBSdetajj", "abs(VBSjet1eta - VBSjet2eta)")
        // .Define("VBSMjj", fInvariantMass, {"VBSjet1pt", "VBSjet1eta", "VBSjet1phi", "VBSjet1mass", "VBSjet2pt", "VBSjet2eta", "VBSjet2phi", "VBSjet2mass"});

        // Reconstruct VBS jets using pair with largest abs(dEta)
        df_out = df_out.Define("Pass_TwoVBSJets", "Sum(goodVBSJets) >= 2") // FIXME could filter here to avoid handling the case of <2 goodVBSJets
                     .Define("goodVBSJets_idxcomb", "Pass_TwoVBSJets ? Combinations(goodVBSJets_eta, 2) : ROOT::VecOps::RVec<ROOT::VecOps::RVec<size_t>>{}") // combinations of unique pairs
                     .Define("goodVBSJets_dEta", "Pass_TwoVBSJets ? Take(goodVBSJets_eta, goodVBSJets_idxcomb[0]) - Take(goodVBSJets_eta, goodVBSJets_idxcomb[1])  : ROOT::VecOps::RVec<float>{}")
                     .Define("idxcomb_maxdEta", "Pass_TwoVBSJets ? ArgMax(abs(goodVBSJets_dEta))  : -999")
                     .Define("idx_vbsj1", "Pass_TwoVBSJets ? goodVBSJets_idxcomb[0][idxcomb_maxdEta]  : -999")
                     .Define("idx_vbsj2", "Pass_TwoVBSJets ? goodVBSJets_idxcomb[1][idxcomb_maxdEta]  : -999")
                     .Define("vbsj1_etaSel_pt", "Pass_TwoVBSJets ? goodVBSJets_pt[idx_vbsj1]  : -999")
                     .Define("vbsj2_etaSel_pt", "Pass_TwoVBSJets ? goodVBSJets_pt[idx_vbsj2]  : -999")
                     .Define("vbsj1_etaSel_eta", "Pass_TwoVBSJets ? goodVBSJets_eta[idx_vbsj1]  : -999")
                     .Define("vbsj2_etaSel_eta", "Pass_TwoVBSJets ? goodVBSJets_eta[idx_vbsj2]  : -999")
                     .Define("vbsj1_etaSel_phi", "Pass_TwoVBSJets ? goodVBSJets_phi[idx_vbsj1]  : -999")
                     .Define("vbsj2_etaSel_phi", "Pass_TwoVBSJets ? goodVBSJets_phi[idx_vbsj2]  : -999")
                     .Define("vbsj1_etaSel_mass", "Pass_TwoVBSJets ? goodVBSJets_mass[idx_vbsj1]  : -999")
                     .Define("vbsj2_etaSel_mass", "Pass_TwoVBSJets ? goodVBSJets_mass[idx_vbsj2]  : -999")
                     .Define("VBSdetajj", "abs(vbsj1_etaSel_eta - vbsj2_etaSel_eta)")
                     .Define("VBSMjj", fInvariantMass, {"vbsj1_etaSel_pt", "vbsj1_etaSel_eta", "vbsj1_etaSel_phi", "vbsj1_etaSel_mass", "vbsj2_etaSel_pt", "vbsj2_etaSel_eta", "vbsj2_etaSel_phi", "vbsj2_etaSel_mass"});

        return df_out;
    }

    RNode eventObjectsReconstruction(RNode df)
    {
        std::cout << " -> Run Run2AllHad3FJ::eventObjectSelections()" << std::endl;
        auto df_out = df.Define("puIDJets_pt", "CorrJet_pt[goodVBSJets || goodAK4Jets]") // used to calculate jet-level pileup weight
                        .Define("puIDJets_eta", "Jet_eta[goodVBSJets || goodAK4Jets]")
                        .Define("MET", "CorrMET_pt")
                        .Define("ht_ak4", "Sum(goodAK4Jets_pt)")
                        .Define("ht_ak8", "Sum(goodAK8Jets_pt)")
                        .Define("ST", "Higgs_pt + V1_pt + V2_pt + CorrMET_pt");

        return df_out;
    }

    RNode runABCDNet(RNode df)
    {
        std::cout << " -> Run Run2AllHad3FJ::runABCDNet()" << std::endl;

        const float WP_SR = 0.97; // Working point for all-hadronic signal region
        auto df_out = df.Define("ABCDNetScore", ABCDNet_Run2AllHad3FJ::evaluate_DNN, {"Higgs_pt", "Higgs_eta", "Higgs_phi", "Higgs_particleNet_mass", "V1_pt", "V1_eta", "V1_phi", "V1_particleNet_mass", "V2_pt", "V2_eta", "V2_phi", "V2_particleNet_mass", "VBSMjj"})
                        .Define("Pass_ABCDNetScoreSR", "ABCDNetScore > " + std::to_string(WP_SR))
                        .Define("Pass_VBSdetajj", "VBSdetajj > 5.0")
                        .Define("Pass_RegionAcut", "Pass_ABCDNetScoreSR && Pass_VBSdetajj")
                        .Define("Pass_RegionBcut", "Pass_ABCDNetScoreSR && !Pass_VBSdetajj")
                        .Define("Pass_RegionCcut", "!Pass_ABCDNetScoreSR && Pass_VBSdetajj")
                        .Define("Pass_RegionDcut", "!Pass_ABCDNetScoreSR && !Pass_VBSdetajj");

        return df_out;
    }

    RNode runAnalysis(RNode df)
    {
        std::cout << " -> Run Run2AllHad3FJ::runAnalysis()" << std::endl;
        
        auto df_out = runCommonSelections(df);
        df_out = triggerSelections(df_out);
        df_out = bosonsReconstruction(df_out);
        df_out = AK4JetsSelections(df_out);
        df_out = VBSJetsReconstruction(df_out);
        df_out = eventObjectsReconstruction(df_out);
        df_out = runABCDNet(df_out);

        df_out = df_out.Define("Pass_AtLeast3AK8Jets", "Sum(goodAK8Jets) >=3")
                     .Define("Pass_NoTightBak4Jet", "Sum(goodAK4FromTightBJet) == 0")
                     .Define("Pass_LeadAK8JetPtAbove550", "goodAK8Jets_pt[ptSortedGoodAK8Jets[0]] > 550")
                     .Define("Pass_3BosonCandidatesExist", "Pass_BoostedHiggsCandidateExists && Pass_BoostedV1CandidateExists && Pass_BoostedV2CandidateExists")
                     .Define("Pass_NoAK8JetsOverlap", "(HV1DeltaR>0.8) && (HV2DeltaR>0.8) && (V1V2DeltaR>0.8)")
                     // Loose signal region before ABCDNet
                     .Define("Pass_LooseSR", "(HiggsScore>0.5) && (V1Score>0.3) && (V2Score>0.3)")
                     // Tight signal recion after ABCDNet
                     .Define("Pass_TightSR", "(HiggsScore>0.8) && (V1Score>0.8) && (V2Score>0.7)");

        df_out = df_out.Define("passCut1", "Pass_EventFilters")
                     .Define("passCut2", "passCut1 && Pass_Triggers")
                     .Define("passCut3", "passCut2 && Pass_AtLeast3AK8Jets")
                     .Define("passCut4", "passCut3 && Pass_LeadAK8JetPtAbove550") 
                     .Define("passCut5", "passCut4 && Pass_3BosonCandidatesExist")
                     .Define("passCut6", "passCut5 && Pass_NoAK8JetsOverlap")
                     .Define("passCut7", "passCut6 && Pass_TwoVBSJets")
                     .Define("passCut8", "passCut7 && Pass_NoTightBak4Jet")
                     .Define("passCut9", "passCut8 && Pass_LooseSR")
                     .Define("passCut10", "passCut8 && Pass_TightSR")
                     .Define("RegionA", "passCut10 && Pass_RegionAcut")
                     .Define("RegionB", "passCut10 && Pass_RegionBcut")
                     .Define("RegionC", "passCut10 && Pass_RegionCcut")
                     .Define("RegionD", "passCut10 && Pass_RegionDcut");
        return df_out;
    }

} // Run2AllHad3FJ