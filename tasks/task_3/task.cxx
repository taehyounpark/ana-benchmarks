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
  std::vector<std::string> tree_files{"Run2012B_SingleMu.root"};
  std::string tree_name = "Events";
  auto ds = df.load(dataset::input<HepQ::Tree>(tree_files,tree_name));
  auto jets_pt = ds.read(dataset::column<VecF>("Jet_pt"));
  auto jets_pt_sel = df.define(column::expression([](VecF const& pts){ return pts[ROOT::VecOps::abs(pts) < 1.0];}))(jets_pt);
  auto cut_jet = df.filter(column::expression([](VecF const& pts){return pts.size();}))(jets_pt_sel);
  auto jets_pt_hist = df.get(query::output<HepQ::Hist<1,VecF>>("jets_pt",45,15,60)).fill(jets_pt_sel).at(cut_jet);
  TCanvas c;
  jets_pt_hist->Draw();
  c.SaveAs("task_3.png");
}

int main(int argc, char **argv) {
  int nthreads = 0;
  if (argc == 2) {
    nthreads = strtol(argv[1], nullptr, 0);
  }
  auto tic = std::chrono::steady_clock::now();
  task(nthreads);
  auto toc = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = toc - tic;
  std::cout << "used threads = " << nthreads
            << ", elapsed time = " << elapsed_seconds.count() << "s"
            << std::endl;
}