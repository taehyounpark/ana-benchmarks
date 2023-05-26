#include <cstdlib>
#include <chrono>
#include "ana/analysis.h"

#include "TCanvas.h"

#include "rootana/Tree.h"
#include "rootana/Hist.h"

template <typename T>
using RVec = ROOT::RVec<T>;
using RVecUI = RVec<unsigned int>;
using RVecI = RVec<int>;
using RVecF = RVec<float>;
using RVecD = RVec<double>;

using cut = ana::selection::cut;
using weight = ana::selection::weight;

void task(int n) {
  ana::multithread::enable(n);
  auto ds = ana::analysis<Tree>({"Run2012B_SingleMu.root"}, "Events");
  auto met = ds.read<float>("MET_pt");
  auto all = ds.filter<cut>("all")(ds.constant(true));
  auto met_hist = ds.book<Hist<1,float>>("met",100,0,200).fill(met).at(all);
  TCanvas c;
  met_hist->Draw();
  c.SaveAs("task_1.pdf");
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