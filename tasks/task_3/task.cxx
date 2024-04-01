#include "HepQuery/Tree.h"
#include "HepQuery/Hist.h"

#include "TCanvas.h"
#include "Math/Vector4D.h"
#include "ROOT/RVec.hxx"

using XYZTVector = ROOT::Math::XYZTVector;
using PtEtaPhiMVector = ROOT::Math::PtEtaPhiMVector;

template <typename T> using Vec = ROOT::RVec<T>;
using VecUI = Vec<unsigned int>;
using VecI = Vec<int>;
using VecF = Vec<float>;
using VecD = Vec<double>;

#include "queryosity.h"

using dataflow = queryosity::dataflow;
namespace multithread = queryosity::multithread;
namespace dataset = queryosity::dataset;
namespace column = queryosity::column;
namespace query = queryosity::query;
namespace systematic = queryosity::systematic;

#include <chrono>
#include <functional>
#include <algorithm>
#include <cstdlib>
void task(int n) {
  dataflow df(multithread::enable(n));
  auto tree_files = std::vector<std::string>{"Run2012B_SingleMu.root"};
  std::string tree_name = "Events";
  auto ds = df.load(dataset::input<HepQ::Tree>(tree_files,tree_name));
  auto n_jet = df.read<unsigned int>("nJet");
  auto jets_pt = df.read<VecF>("Jet_pt");
  auto jets_eta = df.read<VecF>("Jet_eta");
  auto jets_phi = df.read<VecF>("Jet_phi");
  auto jets_m = df.read<VecF>("Jet_mass");
  auto jets_pt_sel = jets_pt[jets_eta > df.constant(-1.0) && jets_eta < df.constant(1.0)];
  auto all = df.filter<cut>("all")(df.constant(true));
  auto jets_pt_hist = df.book<HepQ::Hist<1,VecF>>("jets_pt",45,15,60).fill(jets_pt_sel).at(all);
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
  std::cout << "used threads = " << multithread::concurrency() << ", elapsed time = " << elapsed_seconds.count() << "s" << std::endl;
}