#include "HepQuery/Hist.h"
#include "HepQuery/Tree.h"

#include "Math/Vector4D.h"
#include "ROOT/RVec.hxx"
#include "TCanvas.h"

using XYZTVector = ROOT::Math::XYZTVector;
using PtEtaPhiMVector = ROOT::Math::PtEtaPhiMVector;
constexpr static unsigned int PLACEHOLDER_VALUE = 99999;

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

auto concatf = [](VecF const &a, VecF const &b) { return Concatenate(a, b); };
auto concati = [](const ROOT::RVec<int> &a, const ROOT::RVec<int> &b) {
  return Concatenate(a, b);
};

auto transverse_mass = [](VecF const &Lepton_pt, VecF const &Lepton_phi,
                          float MET_pt, float MET_phi, unsigned int idx) {
  return sqrt(2.0 * Lepton_pt[idx] * MET_pt *
              (1.0 - cos(ROOT::VecOps::DeltaPhi(MET_phi, Lepton_phi[idx]))));
};

auto lepton_flavour = [](unsigned int nMuon, unsigned int nElectron) {
  return Concatenate(ROOT::RVec<int>(nMuon, 0), ROOT::RVec<int>(nElectron, 1));
};

unsigned int additional_lepton_idx(Vec<float> const &pt, Vec<float> const &eta,
                                   Vec<float> const &phi,
                                   Vec<float> const &mass,
                                   Vec<int> const &charge,
                                   Vec<int> const &flavour) {
  const auto c = Combinations(pt, 2);
  float best_mass = PLACEHOLDER_VALUE;
  unsigned int best_i1 = PLACEHOLDER_VALUE;
  unsigned int best_i2 = PLACEHOLDER_VALUE;
  const auto z_mass = 91.2;
  const auto make_p4 = [&](std::size_t idx) {
    return ROOT::Math::PtEtaPhiMVector(pt[idx], eta[idx], phi[idx], mass[idx]);
  };

  for (auto i = 0u; i < c[0].size(); i++) {
    const auto i1 = c[0][i];
    const auto i2 = c[1][i];
    if (charge[i1] == charge[i2])
      continue;
    if (flavour[i1] != flavour[i2])
      continue;
    const auto tmp_mass = (make_p4(i1) + make_p4(i2)).mass();
    if (std::abs(tmp_mass - z_mass) < std::abs(best_mass - z_mass)) {
      best_mass = tmp_mass;
      best_i1 = i1;
      best_i2 = i2;
    }
  }

  if (best_i1 == PLACEHOLDER_VALUE)
    return PLACEHOLDER_VALUE;

  float max_pt = -999;
  unsigned int lep_idx = PLACEHOLDER_VALUE;
  for (auto i = 0u; i < pt.size(); i++) {
    if (i != best_i1 && i != best_i2 && pt[i] > max_pt) {
      max_pt = pt[i];
      lep_idx = i;
    }
  }

  return lep_idx;
}

void task(int n) {

  dataflow df(multithread::enable(n));
  auto tree_files = std::vector<std::string>{"Run2012B_SingleMu.root"};
  std::string tree_name = "Events";
  auto ds = df.load(dataset::input<HepQ::Tree>(tree_files, tree_name));

  auto n_muon = ds.read(dataset::column<unsigned int>("nMuon"));
  auto mus_pt = ds.read(dataset::column<VecF>("Muon_pt"));
  auto mus_eta = ds.read(dataset::column<VecF>("Muon_eta"));
  auto mus_phi = ds.read(dataset::column<VecF>("Muon_phi"));
  auto mus_m = ds.read(dataset::column<VecF>("Muon_mass"));
  auto mus_q = ds.read(dataset::column<VecI>("Muon_charge"));

  auto n_elec = ds.read(dataset::column<unsigned int>("nElectron"));
  auto els_pt = ds.read(dataset::column<VecF>("Electron_pt"));
  auto els_eta = ds.read(dataset::column<VecF>("Electron_eta"));
  auto els_phi = ds.read(dataset::column<VecF>("Electron_phi"));
  auto els_m = ds.read(dataset::column<VecF>("Electron_mass"));
  auto els_q = ds.read(dataset::column<VecI>("Electron_charge"));

  auto met_pt = ds.read(dataset::column<float>("MET_pt"));
  auto met_phi = ds.read(dataset::column<float>("MET_phi"));

  auto leps_pt = df.define(column::expression(concatf))(mus_pt, els_pt);
  auto leps_eta = df.define(column::expression(concatf))(mus_eta, els_eta);
  auto leps_phi = df.define(column::expression(concatf))(mus_phi, els_phi);
  auto leps_m = df.define(column::expression(concatf))(mus_m, els_m);
  auto leps_q = df.define(column::expression(concati))(mus_q, els_q);
  auto leps_type =
      df.define(column::expression(lepton_flavour))(n_muon, n_elec);

  auto add_lep_idx = df.define(column::expression(additional_lepton_idx))(
      leps_pt, leps_eta, leps_phi, leps_m, leps_q, leps_type);
  auto mt = df.define(column::expression(transverse_mass))(
      leps_pt, leps_phi, met_pt, met_phi, add_lep_idx);

  auto three = df.define(column::constant<unsigned int>(3));
  auto cut_3l_sfos = df.filter((n_muon + n_elec) >= three)
                         .filter(column::expression([](unsigned int idx) {
                           return (idx != PLACEHOLDER_VALUE);
                         }))(add_lep_idx);

  auto mt_hist = df.get(query::output<HepQ::Hist<1, float>>("mt", 100, 0, 200))
                     .fill(mt)
                     .at(cut_3l_sfos);

  TCanvas c;
  mt_hist->Draw();
  c.SaveAs("task_8.png");
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