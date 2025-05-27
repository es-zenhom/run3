#include "commonSelections.h"

RNode runCommonSelections(RNode df)
{
    auto df_out = flagSelections(df);
    df_out = AK8JetsPreselection(df_out);
    df_out = AK4JetsPreselection(df_out);
    df_out = VBSJetsPreselection(df_out);

    return df_out;
}

RNode flagSelections(RNode df)
{
    std::cout << " -> Run commonSelections::flagSelections()" << std::endl;
    auto df_out = df.Define("Pass_EventFilters",
                              "Flag_goodVertices &&"
                              "Flag_HBHENoiseFilter &&"
                              "Flag_HBHENoiseIsoFilter &&"
                              "Flag_EcalDeadCellTriggerPrimitiveFilter &&"
                              "Flag_BadPFMuonFilter &&"
                              "Flag_BadPFMuonDzFilter &&"
                              "Flag_hfNoisyHitsFilter &&"
                              "Flag_eeBadScFilter &&"
                              "( (is2016) || Flag_ecalBadCalibFilter) &&"           // apply only to 2017, 2018
                              "( (!isData) || Flag_globalSuperTightHalo2016Filter)" // apply only to data
    );
    return df_out;
}

/*
 *   AK8 Jets selection
 */
RNode AK8JetsPreselection(RNode df)
{
    std::cout << " -> Run commonSelections::AK8JetsPreselection()" << std::endl;
    // Select good AK8 jets
    auto df_out = df.Define("goodAK8Jets",
                          "CorrFatJet_pt > 300 && "
                          "abs(FatJet_eta) < 2.5 && "
                          "CorrFatJet_mass > 50 && "
                          "FatJet_msoftdrop > 40 && "
                          "FatJet_jetId > 0")
                   .Define("FatJet_HbbScore", "FatJet_particleNetMD_Xbb / (FatJet_particleNetMD_Xbb + FatJet_particleNetMD_QCD)")
                    .Define("FatJet_WqqScore", "(FatJet_particleNetMD_Xcc + FatJet_particleNetMD_Xqq) / (FatJet_particleNetMD_Xcc + FatJet_particleNetMD_Xqq + FatJet_particleNetMD_QCD)")
                    .Define("goodAK8Jets_pt", "CorrFatJet_pt[goodAK8Jets]")
                    .Define("goodAK8Jets_eta", "FatJet_eta[goodAK8Jets]")
                    .Define("goodAK8Jets_phi", "FatJet_phi[goodAK8Jets]")
                    .Define("goodAK8Jets_mass", "CorrFatJet_mass[goodAK8Jets]")
                    .Define("goodAK8Jets_msoftdrop", "FatJet_msoftdrop[goodAK8Jets]")
                    .Define("goodAK8Jets_particleNet_mass", "FatJet_particleNet_mass[goodAK8Jets]")
                    .Define("goodAK8Jets_HbbScore", "FatJet_HbbScore[goodAK8Jets]")
                    .Define("goodAK8Jets_WqqScore", "FatJet_WqqScore[goodAK8Jets]")
                    .Define("goodAK8Jets_nConstituents", "FatJet_nConstituents[goodAK8Jets]")
                    .Define("ht_goodAK8Jets", "Sum(goodAK8Jets_pt)")
                    .Define("n_goodAK8Jets", "Sum(goodAK8Jets)")
                    .Define("ptSortedGoodAK8Jets", "Argsort(-goodAK8Jets_pt)"); 

    return df_out;
}


/*
 *   AK4 and VBS Jets selection
 */
RNode AK4JetsPreselection(RNode df)
{
    std::cout << " -> Run commonSelections::AK4JetsPreselection()" << std::endl;
    auto df_out = df.Define("ak4tightBjetScore", tightDFBtagWP, {"sample_year"})
                      .Define("ak4mediumBjetScore", mediumDFBtagWP, {"sample_year"})
                      .Define("ak4looseBjetScore", looseDFBtagWP, {"sample_year"})
                      .Define("Jet_isTightBTag", "Jet_btagDeepFlavB > ak4tightBjetScore")
                      .Define("Jet_isMediumBTag", "Jet_btagDeepFlavB > ak4mediumBjetScore")
                      .Define("Jet_isLooseBTag", "Jet_btagDeepFlavB > ak4looseBjetScore")
                      .Define("Jet_minDrFromAnyGoodAK8Jet", dRfromClosestJet, {"Jet_eta", "Jet_phi", "goodAK8Jets_eta", "goodAK8Jets_phi"})
                      .Define("goodAK4Jets", "CorrJet_pt >= 20 && "
                                           "((is2016 && Jet_jetId >= 1) || (!is2016 && Jet_jetId >= 2)) && "
                                           "(CorrJet_pt >= 50 || (CorrJet_pt < 50 && Jet_puId != 0))"
                                           )
                      .Define("goodAK4Jets_pt", "CorrJet_pt[goodAK4Jets]")
                      .Define("goodAK4Jets_eta", "Jet_eta[goodAK4Jets]")
                      .Define("goodAK4Jets_phi", "Jet_phi[goodAK4Jets]")
                      .Define("goodAK4Jets_mass", "CorrJet_mass[goodAK4Jets]")
                      .Define("goodAK4Jets_isTightBTag", "Jet_isTightBTag[goodAK4Jets]")
                      .Define("goodAK4Jets_isMediumBTag", "Jet_isMediumBTag[goodAK4Jets]")
                      .Define("goodAK4Jets_isLooseBTag", "Jet_isLooseBTag[goodAK4Jets]")
                      .Define("ht_goodAK4Jets", "Sum(CorrJet_pt[goodAK4Jets])")
                      .Define("n_goodAK4Jets", "Sum(goodAK4Jets)")
                      .Define("ptSortedGoodAK4Jets", "Argsort(-CorrJet_pt)") //checkme
                      .Define("goodAK4Jets_minDrFromAnyGoodAK8Jet", dRfromClosestJet, {"goodAK4Jets_eta", "goodAK4Jets_phi", "goodAK8Jets_eta", "goodAK8Jets_phi"})
                      .Define("goodAK4Jets_passAK8OverlapRemoval", "goodAK4Jets_minDrFromAnyGoodAK8Jet>0.8")
                      .Define("n_goodAK4JetsWithAK8OverlapRemoval", "Sum(goodAK4Jets_passAK8OverlapRemoval)"); 

    return df_out;
}

RNode VBSJetsPreselection(RNode df)
{
    std::cout << " -> Run commonSelections::VBSJetsPreselection()" << std::endl;
    auto df_out = df.Define("goodVBSJets", "CorrJet_pt >= 30 && "
                                           "abs(Jet_eta) < 4.7 && "
                                           "((is2016 && Jet_jetId >= 1) || (!is2016 && Jet_jetId >= 2)) && "
                                           "(CorrJet_pt >= 50 || (CorrJet_pt < 50 && Jet_puId != 0))")
                      .Define("goodVBSJets_pt", "CorrJet_pt[goodVBSJets]")
                      .Define("goodVBSJets_eta", "Jet_eta[goodVBSJets]")
                      .Define("goodVBSJets_phi", "Jet_phi[goodVBSJets]")
                      .Define("goodVBSJets_mass", "CorrJet_mass[goodVBSJets]");
    return df_out;
}