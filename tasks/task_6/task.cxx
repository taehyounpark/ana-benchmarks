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
namespace query = queryosity::query;
namespace systematic = queryosity::systematic;

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <functional>

class TopTriJet : public column::definition<Vec<std::size_t>(Vec<XYZTVector>)> {
public:
  TopTriJet(float top_mass) : m_top_mass(top_mass){};
  virtual Vec<std::size_t>
  evaluate(column::observable<Vec<XYZTVector>> jets_p4) const override {
    constexpr std::size_t n = 3;
    float distance = 1e9;
    std::size_t idx1 = 0, idx2 = 1, idx3 = 2;
    for (std::size_t i = 0; i <= jets_p4->size() - n; ++i) {
      auto p1 = (*jets_p4)[i];
      for (std::size_t j = i + 1; j <= jets_p4->size() - n + 1; ++j) {
        auto p2 = (*jets_p4)[j];
        for (std::size_t k = j + 1; k <= jets_p4->size() - n + 2; ++k) {
          auto p3 = (*jets_p4)[k];
          const auto candidate_mass = (p1 + p2 + p3).mass();
          const auto candidate_distance = std::abs(candidate_mass - m_top_mass);
          if (candidate_distance < distance) {
            distance = candidate_distance;
            idx1 = i;
            idx2 = j;
            idx3 = k;
          }
        }
      }
    }
    return {idx1, idx2, idx3};
  }

protected:
  const float m_top_mass;
};

auto get_trijet_pt = [](Vec<float> const &pt, Vec<float> const &eta,
                        Vec<float> const &phi, Vec<float> const &mass,
                        Vec<std::size_t> const &idx) -> float {
  const auto p1 = ROOT::Math::PtEtaPhiMVector(pt[idx[0]], eta[idx[0]],
                                              phi[idx[0]], mass[idx[0]]);
  const auto p2 = ROOT::Math::PtEtaPhiMVector(pt[idx[1]], eta[idx[1]],
                                              phi[idx[1]], mass[idx[1]]);
  const auto p3 = ROOT::Math::PtEtaPhiMVector(pt[idx[2]], eta[idx[2]],
                                              phi[idx[2]], mass[idx[2]]);
  return (p1 + p2 + p3).pt();
};

auto get_trijet_maxval = [](Vec<float> const &vals,
                            Vec<std::size_t> const &idx) -> float {
  return Max(Take(vals, idx));
};

void task(int n) {
  dataflow df(multithread::enable(n));

  auto tree_files = std::vector<std::string>{"Run2012B_SingleMu.root"};
  std::string tree_name = "Events";
  auto ds = df.load(dataset::input<HepQ::Tree>(tree_files, tree_name));

  auto njets = ds.read(dataset::column<unsigned int>("nJet"));
  auto jets_pt = ds.read(dataset::column<VecF>("Jet_pt"));
  auto jets_eta = ds.read(dataset::column<VecF>("Jet_eta"));
  auto jets_phi = ds.read(dataset::column<VecF>("Jet_phi"));
  auto jets_m = ds.read(dataset::column<VecF>("Jet_mass"));
  auto jets_btag = ds.read(dataset::column<VecF>("Jet_btag"));

  auto cut_3jets = df.filter(
      column::expression([](unsigned int njets) { return njets >= 3; }))(njets);

  auto jets_p4 = df.define(column::expression(
      [](Vec<float> pt, Vec<float> eta, Vec<float> phi, Vec<float> m) {
        return ROOT::VecOps::Construct<XYZTVector>(
            ROOT::VecOps::Construct<PtEtaPhiMVector>(pt, eta, phi, m));
      }))(jets_pt, jets_eta, jets_phi, jets_m);
  auto top_trijet = df.define(column::definition<TopTriJet>(172.5))(jets_p4);
  auto trijet_pt = df.define(column::expression(get_trijet_pt))(
      jets_pt, jets_eta, jets_phi, jets_m, top_trijet);
  auto trijet_maxbtag =
      df.define(column::expression(get_trijet_maxval))(jets_btag, top_trijet);

  auto trijet_pt_hist =
      df.get(query::output<HepQ::Hist<1, float>>("trijet_pt", 100, 15, 40))
          .fill(trijet_pt)
          .at(cut_3jets);
  auto trijet_maxbtag_hist =
      df.get(query::output<HepQ::Hist<1, float>>("trijet_maxbtag", 100, 0, 1))
          .fill(trijet_maxbtag)
          .at(cut_3jets);

  TCanvas c;
  c.Divide(2, 1);
  c.cd(1);
  trijet_pt_hist->Draw();
  c.cd(2);
  trijet_maxbtag_hist->Draw();
  c.SaveAs("task_6.png");
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