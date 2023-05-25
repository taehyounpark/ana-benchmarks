#include "ana/analysis.h"
#include "ana/aggregate.h"

#include "TCanvas.h"
#include "ROOT/RVec.hxx"
#include "Math/Vector4D.h"

#include "rootana/Tree.h"
#include "rootana/Histogram.h"

template <typename T>
using RVec = ROOT::RVec<T>;
using RVecUI = RVec<unsigned int>;
using RVecI = RVec<int>;
using RVecF = RVec<float>;
using RVecD = RVec<double>;
using FourVector = ROOT::Math::PtEtaPhiMVector;

using cut = ana::selection::cut;
using weight = ana::selection::weight;

class dimuon_invariant_masses : public ana::column::definition<RVecF(RVecF pt, RVecF eta, RVecF phi, RVecF m, RVecF q)>
{
public:
  virtual RVecF evaluate(ana::observable<RVecF> pt, ana::observable<RVecF> eta, ana::observable<RVecF> phi, ana::observable<RVecF> m, ana::observable<RVecF> q) const override
  {
    ROOT::RVec<float> masses;
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

int main() {
  
  ana::multithread::enable();
  auto ds = ana::analysis<Tree>({"Run2012B_SingleMu.root"}, "Events");

  auto met = ds.read<float>("MET_pt");

  auto nmuons = ds.read<int>("nMuon");
  auto muons_pt = ds.read<RVecF>("Muon_pt");
  auto muons_eta = ds.read<RVecF>("Muon_eta");
  auto muons_phi = ds.read<RVecF>("Muon_phi");
  auto muons_m = ds.read<RVecF>("Muon_mass");
  auto muons_q = ds.read<RVecF>("Muon_charge");

  auto dimuons_m = ds.define<dimuon_invariant_masses>()(muons_pt,muons_eta,muons_phi,muons_m,muons_q);
  auto cut_dimuon_os_60m120 = ds.filter<cut>("dimuon_os_60m120",[](RVecF const& dimuons_m){return Sum(dimuons_m > 60 && dimuons_m < 120) > 0;})(dimuons_m);

  auto met_hist = ds.book<Histogram<1,float>>("met",100,0,200).fill(met).at(cut_dimuon_os_60m120);

  TCanvas c;
  // met_hist->Draw();
  c.SaveAs("task_4.pdf");

}