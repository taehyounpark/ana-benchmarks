#include "AnaQuery/Hist.h"
#include "AnaQuery/Tree.h"

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

#include <queryosity.hpp>

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
void task(int n) {
  dataflow df(multithread::enable(n));
  std::vector<std::string> tree_files{"Run2012B_SingleMu.root"};
  std::string tree_name = "Events";
  auto ds = df.load(dataset::input<AnaQ::Tree>(tree_files, tree_name));
  auto met = ds.read(dataset::column<float>("MET_pt"));
  auto all = df.filter(column::constant(true));
  auto met_hist =
      df.get(query::output<AnaQ::Hist<1, float>>("met", 100, 0, 200))
          .fill(met)
          .at(all);
  TCanvas c;
  met_hist->Draw();
  c.SaveAs("task_1.png");
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