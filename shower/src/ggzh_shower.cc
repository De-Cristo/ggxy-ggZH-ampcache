// ggzh-shower: shower POWHEG ggxy_ggZH events (LHE) with Pythia 8 and analyse them with Rivet
// in the same process (no HepMC file on disk).
//
//   ggzh-shower --lhe <file.lhe> --out <name> [--cmnd <file.cmnd>]... [--set "Key = value"]...
//               [--nevents N] [--seed S] [--analysis <Rivet analysis[:options]>]... [--no-beam-check]
//
// The settings come from the command files (cards/*.cmnd, read in order) and --set lines.
// With POWHEG:veto > 0 the POWHEG hooks (Pythia8Plugins/PowhegHooks.h) veto shower emissions
// harder than the POWHEG one (NLO+PS); otherwise the events are showered as they are (LO+PS).
// The LHE weights (POWHEG <rwgt> block) and Pythia's shower-variation weights reach Rivet
// through Pythia's weight container and the HepMC3 converter.  Output: <name>.yoda.
// Pythia's own Rivet bridge (Pythia8Plugins/Pythia8Rivet.h of 8.315) does not compile with
// Rivet 4.1, so the HepMC3 event is handed to Rivet::AnalysisHandler directly.
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/PowhegHooks.h"
#include "Pythia8Plugins/HepMC3.h"
#include "Rivet/AnalysisHandler.hh"

#include <iostream>
#include <string>
#include <vector>

using namespace Pythia8;

static void usage() {
  std::cerr << "usage: ggzh-shower --lhe <file.lhe> --out <name> [--cmnd <file>]... [--set \"Key = value\"]...\n"
               "                   [--nevents N] [--seed S] [--analysis <name[:options]>]... [--no-beam-check]\n";
}

int main(int argc, char* argv[]) {
  std::string lhe, out;
  std::vector<std::string> cmnds, sets, analyses;
  long nEvents = -1;
  long seed = -1;
  bool beamCheck = true;     // off for tests without parton level (no beam remnants in the record)
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    auto next = [&](const char* what) -> std::string {
      if (i + 1 >= argc) { std::cerr << "missing value for " << what << "\n"; usage(); std::exit(2); }
      return argv[++i];
    };
    if (a == "--lhe") lhe = next("--lhe");
    else if (a == "--out") out = next("--out");
    else if (a == "--cmnd") cmnds.push_back(next("--cmnd"));
    else if (a == "--set") sets.push_back(next("--set"));
    else if (a == "--analysis") analyses.push_back(next("--analysis"));
    else if (a == "--nevents") nEvents = std::stol(next("--nevents"));
    else if (a == "--seed") seed = std::stol(next("--seed"));
    else if (a == "--no-beam-check") beamCheck = false;
    else { std::cerr << "unknown argument " << a << "\n"; usage(); return 2; }
  }
  if (lhe.empty() || out.empty()) { usage(); return 2; }
  if (analyses.empty()) analyses = {"SimplifiedTemplateCrossSections:HiggsProdMode=GG2ZH", "GGZH_NLOPS"};

  Pythia pythia;
  for (const auto& c : cmnds)
    if (!pythia.readFile(c)) { std::cerr << "cannot read " << c << "\n"; return 3; }
  pythia.readString("Beams:frameType = 4");
  pythia.readString("Beams:LHEF = " + lhe);
  if (seed >= 0) {
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = " + std::to_string(seed % 900000000));
  }
  for (const auto& s : sets)
    if (!pythia.readString(s)) { std::cerr << "bad setting " << s << "\n"; return 3; }

  // NLO+PS: veto shower emissions harder than the POWHEG emission
  shared_ptr<PowhegHooks> powhegHooks;
  const bool veto = pythia.settings.mode("POWHEG:veto") > 0;
  if (veto) {
    powhegHooks = make_shared<PowhegHooks>();
    pythia.setUserHooksPtr((UserHooksPtr)powhegHooks);
  }
  if (!pythia.init()) { std::cerr << "Pythia initialisation failed\n"; return 4; }

  Pythia8ToHepMC toHepMC;
  Rivet::AnalysisHandler rivet;
  rivet.setCheckBeams(beamCheck);
  for (const auto& a : analyses) rivet.addAnalysis(a);
  bool rivetInit = false;

  long iEvent = 0, nFail = 0;
  unsigned long nISRveto = 0, nFSRveto = 0;
  const int nAllowed = pythia.settings.mode("Main:timesAllowErrors");
  while (nEvents < 0 || iEvent < nEvents) {
    if (!pythia.next()) {
      if (pythia.info.atEndOfFile()) break;
      if (++nFail > nAllowed) { std::cerr << "too many failed events\n"; return 5; }
      continue;
    }
    toHepMC.fillNextEvent(pythia);
    if (!rivetInit) { rivet.init(toHepMC.event()); rivetInit = true; }
    rivet.analyze(toHepMC.event());
    if (veto) { nISRveto += powhegHooks->getNISRveto(); nFSRveto += powhegHooks->getNFSRveto(); }
    ++iEvent;
  }
  pythia.stat();
  if (rivetInit) { rivet.finalize(); rivet.writeData(out + ".yoda"); }
  std::cout << "ggzh-shower: " << iEvent << " events, " << nFail << " failed; POWHEG veto "
            << (veto ? "on" : "off") << ", ISR vetoes " << nISRveto << ", FSR vetoes " << nFSRveto
            << "; weights per event " << (rivetInit ? toHepMC.event().weights().size() : 0) << "\n";
  return 0;
}
