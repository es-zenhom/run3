# VBS-VVH Analysis framework (hadronic)
The code in this repository is based on the work of aashayarora for the VBS-VVH 1-lepton channel: https://github.com/aashayarora/vbs-1lep/tree/main. 

## Table of Contents

- [Development workflow with github](#development-workflow-with-github)
- [Running the preselection](#running-the-preselection)
  - [Generate the config file](#generate-the-config-file)
  - [Compile and run](#compile-and-run)
  - [Run the full analysis, including all variations ](#run-the-full-analysis-including-all-variations)
  - [Add a new analyzer ](#add-a-new-analyzer)
- [Produce the datacards](#produce-the-datacards)


## Development workflow with github
Since multiple people will be working on the same codebase, always make changes on your development branch and submit a pull request to merge them into the main branch. 

1. Start by cloning the repository:
```bash
git clone git@github.com:mmfsz/VBS-VVH-analysis.git analysis
cd analysis
```

2. Create your own development branch:
    To ensure the main branch remains stable, always create a new branch for your work after fetching the latest changes:
    ```bash
    git pull origin main
    ```
    Create a new branch for your feature or bug fix:
    ```bash
    git checkout -b dev-<your-git-username>-<branche-name>
    ```
    Replace `<branch-name>` with a descriptive name for your branch, e.g., `dev-<your-git-username>/new-feature`.

3. Make your changes and commit. Use meaningful commit messages to describe your work:
    ```bash
    git add .
    git commit -m "Add a descriptive message about the changes"
    ```

4. Push your branch to the remote repository:
    ```bash
    git push origin <your-branch-name>
    ```

5. Submit a pull request:
  - Go to the repository's page on GitHub.
  - Click Pull Requests tab.
  - Select your branch and compare it with the main branch.
  - Provide a clear title and description of your changes.
  - Submit the pull request.

## Running the preselection
This setup assumes that you are running on UAF. 

```bash
git clone git@github.com:mmfsz/VBS-VVH-analysis.git analysis
cd analysis/preselection/.
source setup.sh
```

### Generate the config file
Check the following variables in `etc/makeConfig.py` are pointing correctly to your signal, background, and data inputs: i) `skims_base_dir`, ii) `skims_job_name`, iii) `file_regex`. 

Create the configuration file `input/input_sig.json` to run over your signal samples:
```bash
python etc/makeConfig.py --categories sig --output input/input_sig.json
```

### Compile and run
```bash
cd analysis/preselection/.
source setup.sh
make clean
make -j4
./bin/runAnalysis -n 4 -i input/{INPUT_CONFIG_NAME}.json -o {OUTPUT_FILE_NAME} --cutflow --outdir {OUTPUT_DIR_NAME} --ana {ANALYZER_TAG} --cut {CUT_NAME}
```

- `INPUT_CONFIG_NAME`: the name of the input configuration file generated in the previous step. For automatic sample type identification, add "sig", "bkg", or "data" in the name. 
- `OUTPUT_FILE_NAME` and `OUTPUT_DIR_NAME`: This code will generate ntuples at the path `{OUTPUT_DIR_NAME}/{OUTPUT_FILE_NAME}.root`. Make sure your output directory is in `/ceph/cms/store/user/$USER/` or in `/data/userdata/`. If you leave `OUTPUT_FILE_NAME` unset, the code will automatically set the output name to `{sig,bkg,data}.root` according to which string is present in the name of the input json config. 
- `ANALYZER_TAG`: the name of the analyzer you want to run to produce the ntuples, e.g. `--ana Run2AllHad3FJ` (see below)
- `CUT_NAME`: the name of the cut you want to apply to filter the events before storing them to the output ntuples. This is optional and depends on what cuts you have defined in the code. 
- `--cutflow`: Optionally store the cutflow at `{OUTPUT_DIR_NAME}/{OUTPUT_FILE_NAME}_cutflow.txt`

### Run the full analysis, including all variations 
Inside `runFullAnalysis.sh`, check that the configuration files `config_file_{sig,bkg,data}` are set correctly. 
```bash
cd analysis/preselection/.
source setup.sh
make clean
make -j8
sh runFullAnalysis.sh {ANALYZER_TAG} {OUTPUT_DIR_TAG}
```
The output files will be stored at `/ceph/cms/store/user/$USER/analysisNtuples/output_{ANALYZER_TAG}_{OUTPUT_DIR_TAG}/`.



### Add a new analyzer 
Common event filters and object selections are stored in `commonSelections.cpp`. 

Every channel should define it's own specific selections built on top of this. The idea is to keep the selection of "good objects" harmonized between the channels. In this example, we will add the new channel `Run2AllHad3FJ`. 

1. Create the files `selections_Run2AllHad3FJ.{cpp,h}`. Define all the content in these files inside the namespace `Run2AllHad3FJ` to avoid name collisions. Here you can refine the object selection, perform object reconstruction, and define cuts. You can also apply filters if you wish. The only requirement is to have a function called `runAnalysis()` that will call the common selections, plus all the local functions to run your analysis:

    ```
    // Make sure these are included at the top of your file
    #include "commonSelections.h"
    #include "selections_Run2AllHad3FJ.h"
    #include "ABCDNet_Run2AllHad3FJ.h"

    //...

    RNode runAnalysis(RNode df) {
        auto df_out = runCommonSelections(df); // Run the common selections defined in commonSelections.cpp
        df_out = triggerSelections(df_out);
        df_out = bosonsReconstruction(df_out);

        //...

        return df_out;
    }        
    ```

    Note that this analyzer includes the channel's ABCDNet information via `ABCDNet_Run2AllHad3FJ.h`. You can use this as an example for your own network. To export the DNN to C++ you can use this function: `https://github.com/cmstas/vbs/blob/main/abcdnet/scripts/export_dnn.py`.  


2. In `main.cpp`, add the call to `Run2AllHad3FJ::runAnalysis()` inside the function `runAnalysis()`:

    ```
    RNode runAnalysis(RNode df, MyArgs args, bool isSignal) {
        if (args.ana == "Run2AllHad3FJ") {
            return Run2AllHad3FJ::runAnalysis(df);
        }
        else if (args.ana == "MY_OTHER_CHANNEL") {
            return MY_OTHER_CHANNEL::runAnalysis(df);
        }
        else{
            std::cerr << "Did not recognize analysis namespace: " << args.ana  << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    ```
    Once compiled, you can now run your analyzer by including the flag  `--ana Run2AllHad3FJ`.

3. Optionally, to be able to print a cutflow add an option for your analyzer in the function `main()` inside `main.cpp`:
    ```
    // Inside main.cpp
    if (args.cutflow)
        {
            if (args.ana == "Run2AllHad3FJ") {
                std::vector<std::string> cutflow_cuts = {
                    "Pass_AtLeast3AK8Jets", "Pass_LeadAK8JetPtAbove550", ...};
                auto cutflow = Cutflow(df_final, cutflow_cuts);
                cutflow.Print(output_dir + "/" + output_file + "_cutflow.txt");
            }
        }
    ```


## Produce the datacards
```
python3 make_datacard.py --indir  /ceph/cms/store/user/$USER/analysisNtuples/output_{ANALYZER_TAG}_{OUTPUT_DIR_TAG}/ --chtag {CHANNEL_TAG} 
```

- `indir`: directory with input files `sig.root`, `data.root`, and `sig_{variation}.root`
- `chtag`: channel tag required to have unique names for rate parameters and some uncertainties, e.g. `AllHad3FJ`
- `lhe_tag`: c2v or c2wc2z
- `year`: optional parameter to produce datacards for one single year

The output datacards will be stored at `analysis/datacard/datacards/output_{ANALYZER_TAG}_{OUTPUT_DIR_TAG}/`.
