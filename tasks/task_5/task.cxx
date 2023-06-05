#include "ana/analysis.h"

#include "TCanvas.h"
#include "ROOT/RVec.hxx"
#include "Math/Vector4D.h"

#include "rootana/Tree.h"
#include "rootana/Hist.h"

template <typename T>
using Vec = ROOT::RVec<T>;
using VecUI = Vec<unsigned int>;
using VecI = Vec<int>;
using VecF = Vec<float>;
using VecD = Vec<double>;
using FourVector = ROOT::Math::PtEtaPhiMVector;

using cut = ana::selection::cut;
using weight = ana::selection::weight;

class dimuon_invariant_masses : public ana::column::definition<VecF(VecF pt, VecF eta, VecF phi, VecF m, VecF q)>
{
public:
  virtual VecF evaluate(ana::observable<VecF> pt, ana::observable<VecF> eta, ana::observable<VecF> phi, ana::observable<VecF> m, ana::observable<VecF> q) const override
  {
    VecF masses;
    const auto c = ROOT::VecOps::Combinations(*pt, 2);
    for (auto i = 0u; i < c[0].size(); i++) {
      const auto i1 = c[0][i];
      const auto i2 = c[1][i];
      if (q->at(i1) == q->at(i2)) continue;
      const FourVector p1(pt->at(i1), eta->at(i1), phi->at(i1), m->at(i1));
      const FourVector p2(pt->at(i2), eta->at(i2), phi->at(i2), m->at(i2));
      masses.push_back((p1 + p2).mass());
    }
    return masses;
  }
};

void task(int n) {
  ana::multithread::enable(n);
  auto df = ana::dataflow<Tree>({"Run2012B_SingleMu.root"}, "Events");
  auto met = df.read<float>("MET_pt");
  auto nmuons = df.read<unsigned int>("nMuon");
  auto muons_pt = df.read<VecF>("Muon_pt");
  auto muons_eta = df.read<VecF>("Muon_eta");
  auto muons_phi = df.read<VecF>("Muon_phi");
  auto muons_m = df.read<VecF>("Muon_mass");
  auto muons_q = df.read<VecI>("Muon_charge");
  auto dimuons_m = df.define<dimuon_invariant_masses>()(muons_pt,muons_eta,muons_phi,muons_m,muons_q);

  // require 2 muons beforehand to ensure combinatorics can work
  auto cut_dimuon = df.filter<cut>("dimuon")(nmuons >= df.constant<unsigned int>(2));
  // do the combinatorics
  auto cut_dimuon_os_60m120 = cut_dimuon.filter<cut>("dimuon_os_60m120",[](VecF const& dimuons_m){return Sum(dimuons_m > 60 && dimuons_m < 120) > 0;})(dimuons_m);

  auto met_hist = df.book<Hist<1,float>>("met",100,0,200).fill(met).at(cut_dimuon_os_60m120);
  TCanvas c;
  met_hist->Draw();
  c.SaveAs("task_5.png");
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