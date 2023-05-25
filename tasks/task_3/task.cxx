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
  auto jets_pt = ds.read<RVecF>("Jet_pt");
  auto jets_eta = ds.read<RVecF>("Jet_eta");
  auto jets_pt_sel = jets_pt[jets_eta < ds.constant(1.0)];
  auto all = ds.filter<cut>("all")(ds.constant(true));
  auto jets_pt_hist = ds.book<Histogram<1,RVecF>>("jets_pt",100,0,200).fill(jets_pt).at(all);

  TCanvas c;
  jets_pt_hist->Draw();
  c.SaveAs("task_2.pdf");

}