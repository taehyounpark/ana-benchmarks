#include "HepQuery/Hist.h"
#include "HepQuery/Tree.h"

#include "Math/Vector4D.h"
#include "ROOT/RVec.hxx"
#include "TCanvas.h"

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
namespace selection = queryosity::selection;
namespace query = queryosity::query;
namespace systematic = queryosity::systematic;

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <functional>

class DimuInvMassComb : public column::definition<VecF(
                            VecF pt, VecF eta, VecF phi, VecF m, VecF q)> {
public:
  virtual VecF evaluate(column::observable<VecF> pt,
                        column::observable<VecF> eta,
                        column::observable<VecF> phi,
                        column::observable<VecF> m,
                        column::observable<VecF> q) const override {
    VecF masses;
    const auto c = ROOT::VecOps::Combinations(*pt, 2);
    for (auto i = 0u; i < c[0].size(); i++) {
      const auto i1 = c[0][i];
      const auto i2 = c[1][i];
      if (q->at(i1) == q->at(i2))
        continue;
      const PtEtaPhiMVector p1(pt->at(i1), eta->at(i1), phi->at(i1), m->at(i1));
      const PtEtaPhiMVector p2(pt->at(i2), eta->at(i2), phi->at(i2), m->at(i2));
      masses.push_back((p1 + p2).mass());
    }
    return masses;
  }
};

void task(int n) {
  dataflow df(multithread::enable(n));
  std::vector<std::string> tree_files{"Run2012B_SingleMu.root"};
  std::string tree_name = "Events";
  auto ds = df.load(dataset::input<HepQ::Tree>(tree_files, tree_name));
  auto met = ds.read(dataset::column<float>("MET_pt"));
  auto nmuons = ds.read(dataset::column<unsigned int>("nMuon"));
  auto muons_pt = ds.read(dataset::column<VecF>("Muon_pt"));
  auto muons_eta = ds.read(dataset::column<VecF>("Muon_eta"));
  auto muons_phi = ds.read(dataset::column<VecF>("Muon_phi"));
  auto muons_m = ds.read(dataset::column<VecF>("Muon_mass"));
  auto muons_q = ds.read(dataset::column<VecI>("Muon_charge"));
  auto dimuons_m = df.define(column::definition<DimuInvMassComb>())(
      muons_pt, muons_eta, muons_phi, muons_m, muons_q);

  auto cut_dimuon = df.filter(
      column::expression([](unsigned int n) { return n >= 2; }))(nmuons);
  auto cut_dimuon_os_60m120 =
      cut_dimuon.filter(column::expression([](VecF const &masses) {
        return Sum(masses > 60 && masses < 120) > 0;
      }))(dimuons_m);

  auto met_hist =
      df.get(query::output<HepQ::Hist<1, float>>("met", 100, 0, 200))
          .fill(met)
          .at(cut_dimuon_os_60m120);
  TCanvas c;
  met_hist->Draw();
  c.SaveAs("task_5.png");
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