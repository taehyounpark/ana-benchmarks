#include <cstdlib>
#include <chrono>
#include <algorithm>
#include <functional>

#include "ana/analysis.h"
#include "ana/aggregate.h"

#include "TCanvas.h"
#include "TLorentzVector.h"
#include "ROOT/RVec.hxx"

#include "rootana/Tree.h"
#include "rootana/Histogram.h"

template <typename T> using Vec = ROOT::RVec<T>;
using VecUI = Vec<unsigned int>;
using VecI = Vec<int>;
using VecF = Vec<float>;
using VecD = Vec<double>;
using P4 = TLorentzVector;

class TriJetLoc : public ana::column::definition<Vec<std::size_t>(Vec<P4>)>
{
  virtual Vec<std::size_t> evaluate(ana::observable<Vec<P4>> jets_p4) const override
  {
    constexpr std::size_t n = 3;
    float distance = 1e9;
    const auto top_mass = 172.5;
    std::size_t idx1 = 0, idx2 = 1, idx3 = 2;
    for (std::size_t i = 0; i <= jets_p4->size() - n; i++) {
      auto p1 = (*jets_p4)[i];
      for (std::size_t j = i + 1; j <= jets_p4->size() - n + 1; j++) {
        auto p2 = (*jets_p4)[j];
        for (std::size_t k = j + 1; k <= jets_p4->size() - n + 2; k++) {
          auto p3 = (*jets_p4)[k];
          const auto candidate_mass = (p1 + p2 + p3).M();
          const auto candidate_distance = std::abs(candidate_mass - top_mass);
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
};

float get_trijet_pt(Vec<P4> const& p4s, Vec<std::size_t> const& idx)
{
  return (p4s[0] + p4s[1] + p4s[2]).Pt();
}

float get_trijet_maxval(Vec<float> const& vals, Vec<std::size_t> const& idx)
{
  auto trijet_vals = std::vector<float>(3,0);
  trijet_vals[0] = vals[idx[0]];
  trijet_vals[1] = vals[idx[1]];
  trijet_vals[2] = vals[idx[2]];
  return *std::max_element(trijet_vals.begin(), trijet_vals.end());
}

using cut = ana::selection::cut;
using weight = ana::selection::weight;

void task(int n) {

  ana::multithread::enable(n);
  auto ds = ana::analysis<Tree>({"Run2012B_SingleMu.root"}, "Events");

  auto njets = ds.read<unsigned int>("nJet");
  auto jets_pt = ds.read<VecF>("Jet_pt");
  auto jets_eta = ds.read<VecF>("Jet_eta");
  auto jets_phi = ds.read<VecF>("Jet_phi");
  auto jets_m = ds.read<VecF>("Jet_mass");
  auto jets_btag = ds.read<VecF>("Jet_btag");

  auto cut_3jets = ds.filter<cut>("3jets")(njets >= ds.constant(3));

  auto jets_p4 = ds.define([](VecF const& pts, VecF const& etas, VecF const& phis, VecF const& ms){
    Vec<P4> p4s(pts.size(), P4());
    for (size_t i=0 ; i<pts.size(); ++i) {
      p4s[i].SetPtEtaPhiM(pts[i],etas[i],phis[i],ms[i]);
    }
    return p4s;
  })(jets_pt, jets_eta, jets_phi, jets_m);
  auto trijet_loc = ds.define<TriJetLoc>()(jets_p4);
  auto trijet_pt = ds.define(std::function(get_trijet_pt))(jets_p4, trijet_loc);
  auto trijet_maxbtag = ds.define(std::function(get_trijet_maxval))(jets_btag, trijet_loc);

  auto trijet_pt_hist = ds.book<Histogram<1,float>>("trijet_pt",25,15,40).fill(trijet_pt).at(cut_3jets);
  auto trijet_maxbtag_hist = ds.book<Histogram<1,float>>("trijet_maxbtag",50,0,1).fill(trijet_maxbtag).at(cut_3jets);
  
  TCanvas c;
  c.Divide(2,1);
  c.cd(1);
  trijet_pt_hist->Draw();
  c.cd(2);
  trijet_maxbtag_hist->Draw();
  c.SaveAs("task_6.pdf");
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