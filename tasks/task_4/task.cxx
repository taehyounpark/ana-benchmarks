#include "ana/analysis.h"
#include "ana/aggregate.h"

#include "TCanvas.h"

#include "rootana/Tree.h"
#include "rootana/Histogram.h"

template <typename T>
using RVec = ROOT::RVec<T>;
using RVecUI = RVec<unsigned int>;
using RVecI = RVec<int>;
using RVecF = RVec<float>;
using RVecD = RVec<double>;

using cut = ana::selection::cut;
using weight = ana::selection::weight;

int main() {
  
  ana::multithread::enable();
  auto ds = ana::analysis<Tree>({"Run2012B_SingleMu.root"}, "Events");

  auto met = ds.read<float>("MET_pt");
  auto jets_pt = ds.read<RVecF>("Jet_pt");
  auto jets_eta = ds.read<RVecF>("Jet_eta");

  auto njets_pt40 = ds.define([](RVecF const& jets_pt){return jets_pt[jets_pt > 40.0].size();})(jets_pt);
  auto cut_2jets = ds.filter<cut>("2jets")(njets_pt40 > ds.constant(2));

  auto met_hist = ds.book<Histogram<1,float>>("met",100,0,200).fill(met).at(cut_2jets);

  TCanvas c;
  met_hist->Draw();
  c.SaveAs("task_4.pdf");

}