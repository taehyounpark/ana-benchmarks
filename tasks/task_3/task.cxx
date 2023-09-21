#include <array>

#include "ana/analogical.h"

#include <ROOT/RVec.hxx>
#include "TCanvas.h"

#include "AnalysisPlugins/Tree.h"
#include "AnalysisPlugins/Hist.h"

using cut = ana::selection::cut;
using weight = ana::selection::weight;

template <typename T> using Vec = ROOT::RVec<T>;
using VecUI = Vec<unsigned int>;
using VecI = Vec<int>;
using VecF = Vec<float>;
using VecD = Vec<double>;

void task(int n) {
  ana::multithread::enable(n);
  auto df = ana::dataflow<Tree>({"Run2012B_SingleMu.root"}, "Events");
  auto n_jet = df.read<unsigned int>("nJet");
  auto jets_pt = df.read<VecF>("Jet_pt");
  auto jets_eta = df.read<VecF>("Jet_eta");
  auto jets_phi = df.read<VecF>("Jet_phi");
  auto jets_m = df.read<VecF>("Jet_mass");
  auto jets_pt_sel = jets_pt[jets_eta > df.constant(-1.0) && jets_eta < df.constant(1.0)];
  auto all = df.filter<cut>("all")(df.constant(true));
  auto jets_pt_hist = df.book<Hist<1,VecF>>("jets_pt",45,15,60).fill(jets_pt_sel).at(all);
  TCanvas c;
  jets_pt_hist->Draw();
  c.SaveAs("task_3.png");
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