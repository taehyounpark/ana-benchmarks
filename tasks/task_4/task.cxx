#include "ana/analysis.h"

#include <ROOT/RVec.hxx>
#include "TCanvas.h"

#include "rootana/Tree.h"
#include "rootana/Hist.h"

template <typename T>
using Vec = ROOT::RVec<T>;
using VecUI = Vec<unsigned int>;
using VecI = Vec<int>;
using VecF = Vec<float>;
using VecD = Vec<double>;

using cut = ana::selection::cut;
using weight = ana::selection::weight;

void task(int n) {
  ana::multithread::enable(n);
  auto df = ana::dataflow<Tree>({"Run2012B_SingleMu.root"}, "Events");
  auto met = df.read<float>("MET_pt");
  auto jets_pt = df.read<VecF>("Jet_pt");
  auto jets_eta = df.read<VecF>("Jet_eta");
  auto njets_pt40 = df.define([](VecF const& jets_pt){return jets_pt[jets_pt > 40.0].size();})(jets_pt);
  auto cut_2jets = df.filter<cut>("2jets")(njets_pt40 >= df.constant(2));
  auto met_hist = df.book<Hist<1,float>>("met",100,0,200).fill(met).at(cut_2jets);
  TCanvas c;
  met_hist->Draw();
  c.SaveAs("task_4.pdf");
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