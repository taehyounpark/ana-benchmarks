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

struct Jets : public ana::column::aggregate<Jets(RVecF pt, RVecD eta, RVecD phi, RVecD m)>
{
};

int main() {
  
  ana::multithread::enable();
  auto ds = ana::analysis<Tree>({"Run2012B_SingleMu.root"}, "Events");
  auto met = ds.read<float>("MET_pt");
  auto all = ds.filter<cut>("all")(ds.constant(true));
  auto met_hist = ds.book<Histogram<1,float>>("met",100,0,200).fill(met).at(all);

  TCanvas c;
  met_hist->Draw();
  c.SaveAs("task_1.pdf");

}