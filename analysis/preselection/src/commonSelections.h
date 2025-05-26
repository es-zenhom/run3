#ifndef commonSelections_H
#define commonSelections_H

#include "ROOT/RDataFrame.hxx"

#include "utils.h"
#include "corrections.h"

using RNode = ROOT::RDF::RNode;

RNode runCommonSelections(RNode df);
RNode flagSelections(RNode df);
RNode AK8JetsPreselection(RNode df_);
RNode AK4JetsPreselection(RNode df);
RNode VBSJetsPreselection(RNode df);

#endif // commonSelections_H