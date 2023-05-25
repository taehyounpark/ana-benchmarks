#include <cstdlib>
#include <chrono>
#include <algorithm>
#include <functional>

#include "ana/analysis.h"
#include "ana/aggregate.h"
#include "ana/definition.h"

#include "TCanvas.h"
#include "ROOT/RVec.hxx"

#include "rootana/Tree.h"
#include "rootana/Histogram.h"

using cut = ana::selection::cut;
using weight = ana::selection::weight;

template <typename T> using Vec = ROOT::RVec<T>;
using VecUI = Vec<unsigned int>;
using VecI = Vec<int>;
using VecF = Vec<float>;
using VecD = Vec<double>;


class JetLepDRReq : public ana::column::definition<VecI(VecF,VecF,VecF,VecF,VecF)>
{
public:
  JetLepDRReq(float minDR, float pt2Min) : m_minDR(minDR), m_pt2Min(pt2Min) {}
  virtual ~JetLepDRReq() = default;
  virtual VecI evaluate(ana::observable<VecF> eta1, ana::observable<VecF> phi1, ana::observable<VecF> pt2, ana::observable<VecF> eta2, ana::observable<VecF> phi2) const override
  {
    VecI mask(eta1->size(), 1);
    if (eta2->size() == 0) {
      return mask;
    }

    const auto ptcut = (*pt2) > m_pt2Min;
    const auto eta2_ptcut = (*eta2)[ptcut];
    const auto phi2_ptcut = (*phi2)[ptcut];
    if (eta2_ptcut.size() == 0) {
      return mask;
    }

    const auto c = ROOT::VecOps::Combinations(*eta1, eta2_ptcut);
    for (auto i = 0u; i < c[0].size(); i++) {
      const auto i1 = c[0][i];
      const auto i2 = c[1][i];
      const auto dr = ROOT::VecOps::DeltaR((*eta1)[i1], eta2_ptcut[i2], (*phi1)[i1], phi2_ptcut[i2]);
      if (dr < m_minDR) mask[i1] = 0;
    }
    return mask;
  }
protected:
  const float m_minDR;
  const float m_pt2Min;
};

void task(int n) {

  ana::multithread::enable(n);
  auto ds = ana::analysis<Tree>({"Run2012B_SingleMu.root"}, "Events");

  auto n_jet = ds.read<unsigned int>("nJet");
  auto jets_pt = ds.read<VecF>("Jet_pt");
  auto jets_eta = ds.read<VecF>("Jet_eta");
  auto jets_phi = ds.read<VecF>("Jet_phi");
  auto jets_m = ds.read<VecF>("Jet_mass");

  auto n_muon = ds.read<unsigned int>("nMuon");
  auto mus_pt = ds.read<VecF>("Muon_pt");
  auto mus_eta = ds.read<VecF>("Muon_eta");
  auto mus_phi = ds.read<VecF>("Muon_phi");

  auto n_elec = ds.read<unsigned int>("nElectron");
  auto els_pt = ds.read<VecF>("Electron_pt");
  auto els_eta = ds.read<VecF>("Electron_eta");
  auto els_phi = ds.read<VecF>("Electron_phi");

  auto jets_ptcut = ds.define([](VecF const& pts){return pts>30;})(jets_pt);
  auto jets_mudr  = ds.define<JetLepDRReq>(0.4,10.0)(jets_eta, jets_phi, mus_pt, mus_eta, mus_phi);
  auto jets_eldr  = ds.define<JetLepDRReq>(0.4,10.0)(jets_eta, jets_phi, els_pt, els_eta, els_phi);
  auto goodjet_mask    = jets_ptcut && jets_mudr && jets_eldr;
  auto goodjet_sumpt = ds.define([](VecI const& good, VecF const& pt) { return Sum(pt[good]); })(goodjet_mask, jets_pt);

  auto cut_1jet = ds.filter<cut>("1jet")(n_jet >= ds.constant(1));
  auto cut_goodjet = cut_1jet.filter<cut>("goodjet",[](VecI const& goodjet){return Sum(goodjet);})(goodjet_mask);

  auto goodjet_sumpt_hist = ds.book<Histogram<1,float>>("goodjet_sumpt",185,15,200).fill(goodjet_sumpt).at(cut_goodjet);

  TCanvas c;
  goodjet_sumpt_hist->Draw();
  c.SaveAs("task_7.pdf");
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