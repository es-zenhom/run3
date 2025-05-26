#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDFHelpers.hxx"
#include <filesystem>

#include "weights.h"
#include "corrections.h"
#include "utils.h"

#include "commonSelections.h"
// ADD HERE ANALYZER
#include "selections_Run2AllHad3FJ.h"

#include "argparser.hpp"

struct MyArgs : public argparse::Args {
    std::string &spec = kwarg("i,input", "spec.json path");
    bool &cutflow = flag("cutflow", "print cutflow");
    bool &VERB = flag("v", "Verbose mode");
    bool &JEC = flag("jec", "JEC");
    bool &JER = flag("jer", "JER");
    bool &JMS = flag("jms", "JMS");
    bool &JMR = flag("jmr", "JMR");
    bool &METUnclustered = flag("met", "MET unclustered");
    int &nthreads = kwarg("n,nthreads", "number of threads").set_default(1);
    std::string &ana = kwarg("a,ana", "Tag of analyzer to use for event selection");
    std::string &output = kwarg("o,output", "output root file name without extension").set_default("");
    std::string &output_dir = kwarg("outdir", "output directory").set_default("./output/");
    std::string &variation = kwarg("var", "variation").set_default("nominal");
    std::string &JERvariation = kwarg("jervar", "JER variation").set_default("nominal");
    std::string &JEC_type = kwarg("jec_type", "JEC type").set_default("");
    std::string &cut = kwarg("cut", "cut on final snapshot").set_default("passCut9");
};


RNode runAnalysis(RNode df, MyArgs args, bool isSignal) {
    if (args.VERB) { std::cout << " -> Running "<< args.ana << "::runAnalysis()" << std::endl; }

    if (args.ana == "Run2AllHad3FJ") {
        return Run2AllHad3FJ::runAnalysis(df);
    }
    else{
        std::cerr << "Did not recognize analysis namespace: " << args.ana  << std::endl;
        std::exit(EXIT_FAILURE);
    }

}

RNode runDataAnalysis(RNode df, MyArgs args) {
    if (args.VERB) { std::cout << " -> Running runDataAnalysis()" << std::endl; }
    auto df_out = defineCorrectedCols(df);
    //df_out = applyPreSelections(df_out);
    df_out = runAnalysis(df_out, args, false);
    df_out = applyDataWeights(df_out);

    return df_out;
}

RNode runMCAnalysis(RNode df, MyArgs args, bool isSignal) {
    if (args.VERB) { std::cout << " -> Running runMCAnalysis()" << std::endl; }

    // corrections
    auto df_out = defineCorrectedCols(df);

    // apply pre preselection corrections
    if (args.JER){
        if (args.VERB) { std::cout << " -> Running JER corrections" << std::endl; }
        df_out = JetEnergyResolution(cset_jerc_2016preVFP, cset_jerc_2016postVFP, cset_jerc_2017, cset_jerc_2018, cset_jer_smear, df_out, args.JERvariation);
    }
    else if (args.METUnclustered) {
        if (args.VERB) { std::cout << " -> Running METUnclustered corrections" << std::endl; }
        df_out = METUnclusteredCorrections(df_out, args.variation);
    }
    else if (args.JEC){
        if (args.VERB) { std::cout << " -> Running JEC corrections" << std::endl; }
        df_out = JetEnergyCorrection(cset_jerc_2016preVFP, cset_jerc_2016postVFP, cset_jerc_2017, cset_jerc_2018, df_out, args.JEC_type, args.variation);
    }

    // run pre-selection
    if (args.VERB) { std::cout << " -> Apply pre-selection" << std::endl; }
    df_out = runAnalysis(df_out, args, isSignal);

    // apply post pre-selection corrections FIXME: check if this nees to be run before abcdNet
    if (args.JMS) {
        if (args.VERB) { std::cout << " -> Running JMS corrections" << std::endl; }
        df_out = JMS_Corrections(cset_jms, df_out, args.variation);
    }
    else if (args.JMR) {
        if (args.VERB) { std::cout << " -> Running JMR corrections" << std::endl; }
        df_out = JMR_Corrections(cset_jmr, df_out, args.variation);
    }

    if (args.VERB) { std::cout << " -> Apply MC weights" << std::endl; }
    df_out = applyMCWeights(df_out);

    return df_out;
}

void create_directory_if_not_exists(const std::string &dir, MyArgs args) {
    std::filesystem::path directory_path(dir);
    // Check if the directory exists
    if (std::filesystem::exists(directory_path)) {
        if (args.VERB) { std::cerr << "Directory already exists: " << dir << std::endl; }
    }
    // Try to create the directory and any missing parent directories
    else if (std::filesystem::create_directories(directory_path)) {
        std::cout << "Directory created: " << dir << std::endl;
    }
    else {
        std::cerr << "Failed to create directory: " << dir << std::endl;
        std::exit(EXIT_FAILURE); 
    }
}

int main(int argc, char** argv) {
    auto args = argparse::parse<MyArgs>(argc, argv);
    std::string input_spec = args.spec;
    std::string output_file = args.output;
    std::string output_dir = args.output_dir;
    create_directory_if_not_exists(output_dir, args);

    if (args.nthreads>1) {
        ROOT::EnableImplicitMT(args.nthreads);
    }
    ROOT::RDataFrame df_in = ROOT::RDF::Experimental::FromSpec(input_spec);
    ROOT::RDF::Experimental::AddProgressBar(df_in);

    // define metadata
    auto df = defineMetadata(df_in);
    

    // run analysis
    bool isData = false;
    bool isSignal = false;
    if (input_spec.find("data") != std::string::npos) {
        isData = true;
        if (output_file.empty()) {
            output_file = "data";
        }
    }
   else if (input_spec.find("bkg") != std::string::npos) {
        if (output_file.empty()) {
            output_file = "bkg";
        }
    }
    else if (input_spec.find("sig") != std::string::npos) {
        if (output_file.empty()) {
            output_file = "sig";
        }
        isSignal = true;
    }           
    else {
        std::cerr << "Could not guess output name from spec name, file must contain sig, bkg or data" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    ROOT::RDF::RNode df_final = (isData)
                                    ? runDataAnalysis(df, args)
                                    : runMCAnalysis(df, args, isSignal);

    if (args.cutflow)
    {
        if (args.ana == "Run2AllHad3FJ") {
            if (args.VERB) { std::cout << " -> Produce cutflow for predefined cuts" << std::endl; }
            std::vector<std::string> cutflow_cuts = {
                "Pass_AtLeast3AK8Jets", "Pass_LeadAK8JetPtAbove550", "Pass_NoTightBak4Jet",
                "Pass_3BosonCandidatesExist", "Pass_NoAK8JetsOverlap", "Pass_TwoVBSJets", "Pass_LooseSR", "Pass_TightSR", 
                "Pass_RegionAcut", "Pass_RegionBcut", "Pass_RegionCcut", "Pass_RegionDcut",
                "passCut1", "passCut2", "passCut3", "passCut4", "passCut5", "passCut6", "passCut7", 
                "passCut8", "passCut9", "passCut10", "RegionA", "RegionB", "RegionC", "RegionD"};
            auto cutflow = Cutflow(df_final, cutflow_cuts);
            cutflow.Print(output_dir + "/" + output_file + "_cutflow.txt");
         }
    }

    if (!args.cut.empty()){
        if (args.VERB) { std::cout << " -> Stored to file events that pass cut:" << args.cut << std::endl; }
        df_final = df_final.Filter(args.cut);
    }

    saveSnapshot(df_final, output_dir, output_file, isData);

    return 0;
}