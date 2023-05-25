#include "ana/analysis.h"
#include "ana/aggregate.h"

#include "TCanvas.h"

#include "rootana/Tree.h"
#include "rootana/Histogram.h"

template <typename T> using RVec = ROOT::RVec<T>;
using RVecUI = RVec<unsigned int>;
using RVecI = RVec<int>;
using RVecF = RVec<float>;
using RVecD = RVec<double>;

using cut = ana::selection::cut;
using weight = ana::selection::weight;

void task(int n) {
  ana::multithread::enable(n);
  auto ds = ana::analysis<Tree>({"Run2012B_SingleMu.root"}, "Events");
  ds.limit_entries(1000000000);
  auto jets_pt = ds.read<RVecF>("Jet_pt");
  auto jets_eta = ds.read<RVecF>("Jet_eta");
  auto jets_pt_sel = jets_pt[jets_eta > ds.constant(-1.0) && jets_eta < ds.constant(1.0)];
  auto all = ds.filter<cut>("all")(ds.constant(true));
  auto jets_pt_hist = ds.book<Histogram<1,RVecF>>("jets_pt",45,15,60).fill(jets_pt_sel).at(all);
  TCanvas c;
  jets_pt_hist->Draw();
  c.SaveAs("task_3.pdf");
}

int main(int argc, char **argv) {
  int nthreads = 0;
  if (argc==2) { nthreads=strtol(argv[1], nullptr, 0); }
  auto tic = std::chrono::steady_clock::now();
  task(nthreads);
  auto toc = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = toc-tic;
  std::cout << "used threads = " << ana::multithread::concurrency() << ", elapsed time = " << elapsed_seconds.count() << "s" << std::endl;
}