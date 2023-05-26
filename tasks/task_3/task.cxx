#include <array>

#include "ana/analysis.h"
#include "ana/representation.h"

#include "TCanvas.h"
#include "ROOT/RVec.hxx"

#include "rootana/Tree.h"
#include "rootana/Hist.h"

using cut = ana::selection::cut;
using weight = ana::selection::weight;

template <typename T> using Vec = ROOT::RVec<T>;
using VecUI = Vec<unsigned int>;
using VecI = Vec<int>;
using VecF = Vec<float>;
using VecD = Vec<double>;

template <typename T> using obs = typename ana::observable<T>;
enum class JetProperty { pt, eta, phi, m };
class Jet : public ana::column::representation<Jet(float,float,float,float)>
{
public:
  float pt() const { return this->value<JetProperty::pt>(); }
};

template<std::size_t... I>
auto make_jet_tuple(ana::analysis<Tree>& ds, 
  ana::analysis<Tree>::lazy<Tree::Branch<VecF>> const& jets_pt,
  ana::analysis<Tree>::lazy<Tree::Branch<VecF>> const& jets_eta,
  ana::analysis<Tree>::lazy<Tree::Branch<VecF>> const& jets_phi,
  ana::analysis<Tree>::lazy<Tree::Branch<VecF>> const& jets_m
) {
  return std::make_tuple(
    std::invoke(
      [&](int i){
        auto idx = ds.constant(i);
        return ds.define<Jet>()(jets_pt[ds.constant(i)], jets_eta[ds.constant(i)], jets_phi[ds.constant(i)], jets_m[ds.constant(i)]);
      },
      I
    )...
  );
}


void task(int n) {
  ana::multithread::enable(n);
  auto ds = ana::analysis<Tree>({"Run2012B_SingleMu.root"}, "Events");
  auto n_jet = ds.read<unsigned int>("nJet");
  auto jets_pt = ds.read<VecF>("Jet_pt");
  auto jets_eta = ds.read<VecF>("Jet_eta");
  auto jets_phi = ds.read<VecF>("Jet_phi");
  auto jets_m = ds.read<VecF>("Jet_mass");

  auto jets_pt_sel = jets_pt[jets_eta > ds.constant(-1.0) && jets_eta < ds.constant(1.0)];

  // make a tuple of up to 10 jets
  auto jet_tuple = make_jet_tuple<0,1,2,3,4,5,5,6,7,8,9>(ds,jets_pt,jets_eta,jets_phi,jets_m);

  // reminder -- each "Jet" is simply a "representation", so nothing is being computed unless we make it to,
  // i.e. doesn't matter if an event only has e.g. 6 jets, as long as we know not to ask for the 7th!
  // which means until we can move even "non-existent" jets as much as we want for the time being.
  using J = Jet const&;
  auto make_jet_list = [](J j1, J j2, J j3, J j4, J j5, J j6, J j7, J j8, J j9, J j10){
    Vec<Jet> jets;
    jets.push_back(j1);
    jets.push_back(j2);
    jets.push_back(j3);
    jets.push_back(j4);
    jets.push_back(j5);
    jets.push_back(j6);
    jets.push_back(j7);
    jets.push_back(j8);
    jets.push_back(j9);
    jets.push_back(j10);
    return jets;
  };
  auto jet_list = ds.define(make_jet_list)(
    std::get<0>(jet_tuple), 
    std::get<1>(jet_tuple), 
    std::get<2>(jet_tuple), 
    std::get<3>(jet_tuple), 
    std::get<4>(jet_tuple), 
    std::get<5>(jet_tuple), 
    std::get<6>(jet_tuple), 
    std::get<7>(jet_tuple), 
    std::get<8>(jet_tuple), 
    std::get<9>(jet_tuple)
  );

  // first njet argument is important:
  // need to know beforehand when to stop!
  auto jet_eta_selection = [](int njet, Vec<Jet> const& jets) {
    Vec<Jet> jets_sel;
    for (size_t ijet=0 ; ijet<njet ; ++ijet) {
      if (std::abs(jets[ijet].value<JetProperty::eta>()) < 1.0) {
        jets_sel.push_back(jets[ijet]);
      }
    }
    return jets_sel;
  };

  // now we have a Vec<Jet> whose eta is < 1.0
  auto jet_list_sel = ds.define(jet_eta_selection)(n_jet, jet_list);

  auto incl = ds.filter<cut>("all")(ds.constant(true));
  auto none = ds.filter<cut>("none")(ds.constant(false));
  auto all = incl || none;

  auto jets_pt_hist = ds.book<Hist<1,VecF>>("jets_pt",45,15,60).fill(jets_pt_sel).at(all);

  TCanvas c;
  jets_pt_hist->Draw();
  c.SaveAs("task_3.pdf");
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