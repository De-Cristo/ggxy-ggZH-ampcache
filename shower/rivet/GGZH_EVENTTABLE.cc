// -*- C++ -*-
// GGZH_EVENTTABLE: one text row per event of gg -> Z(l+l-)H after the shower, for the LO+PS -> NLO+PS
// reweighting study (rwt/).  Objects exactly as in GGZH_NLOPS (shower/rivet/GGZH_NLOPS.cc):
// Z and H are the last copies of the bosons in the event record; jets are anti-kT R = 0.4 on all
// stable particles except the descendants of the Higgs and the Z boson, pT > 30 GeV.
//
// Columns (header line first, then one line per event):
//   ievt                          event counter (0, 1, ...), for the even/odd split
//   w:nominal w:1001 ... w:1009   the nominal weight (Rivet's default weight index) and the nine scale
//                                 weights (named by LHE id or MUR<r>_MUF<f>), when present; w:2000 = NNPDF3.1
//                                 central (PDF weight list of the package), when present
//   ptz yz etaz phiz mz           last-copy Z (mz = m_ll)
//   pth yh etah phih              last-copy H
//   mzh ptzh yzh drzh             Z + H system; drzh with pseudorapidity
//   njets30 ptj1                  jets pT > 30 GeV, no rapidity cut; ptj1 = -1 without a jet
//   njets30eta47 ptj1eta47        the same with |eta| < 4.7
//   cell                          Stage 1.3 cell 3*ptbin + min(njets30, 2), pT(Z) edges 75/150/250/400/600;
//                                 -1 if |y_H| > 2.5 (FWDH); as GGZH_NLOPS's stxs13_fine_own
// Events without a Z or H in the record are written with nan kinematics (counted in the log).
// Option OUT=<file> (default ggzh-events.csv).
#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include <cstdio>
#include <limits>

namespace Rivet {

  class GGZH_EVENTTABLE : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(GGZH_EVENTTABLE);

    void init() {
      declare(FinalState(), "FS");
      const std::string out = getOption("OUT", "ggzh-events.csv");
      _f = std::fopen(out.c_str(), "w");
      if (_f == nullptr) throw Error("GGZH_EVENTTABLE: cannot open " + out);
    }

    Particle lastCopy(const Particle& p) {
      Particle cur = p;
      bool found = true;
      while (found) {
        found = false;
        for (const Particle& c : cur.children())
          if (c.pid() == cur.pid()) { cur = c; found = true; break; }
      }
      return cur;
    }

    bool descendsFrom(const Particle& p, const Particle& anc) {
      HepMC3::ConstGenVertexPtr v = p.genParticle()->production_vertex();
      if (v == nullptr) return false;
      for (const auto& a : HepMCUtils::particles(v, Relatives::ANCESTORS))
        if (a == anc.genParticle()) return true;
      return false;
    }

    void header() {
      const std::vector<std::string>& names = handler().weightNames();
      const size_t inom = handler().defaultWeightIndex();
      _keep.push_back(inom);
      // the nine scale weights: LHE ids 1001..1009 (POWHEG rwl list), which Pythia/Rivet may name
      // MUR<r>_MUF<f> instead; written as w:1001 ... w:1009 in either case
      static const char* ids[]  = {"1001", "1002", "1003", "1004", "1005", "1006", "1007", "1008", "1009"};
      static const char* murf[] = {"MUR1.0_MUF1.0", "MUR1.0_MUF2.0", "MUR1.0_MUF0.5", "MUR2.0_MUF1.0", "MUR2.0_MUF2.0",
                                   "MUR2.0_MUF0.5", "MUR0.5_MUF1.0", "MUR0.5_MUF2.0", "MUR0.5_MUF0.5"};
      for (size_t k = 0; k < 9; ++k)
        for (size_t i = 0; i < names.size(); ++i) {
          if (i == inom) continue;
          const std::string& n = names[i];
          const std::string s(ids[k]);
          const bool byId = n.size() >= s.size() && n.compare(n.size() - s.size(), s.size(), s) == 0;
          if (byId || n == murf[k]) { _keep.push_back(i); _keepName.push_back(ids[k]); break; }
        }
      // the NNPDF3.1 central member (LHE id 2000 in the package's weight lists), when present: lets a
      // PDF4LHC21 sample be used with NNPDF3.1 (the CMS PDF) by taking w:2000 as the nominal weight
      for (size_t i = 0; i < names.size(); ++i) {
        const std::string& n = names[i];
        if (i != inom && n.size() >= 4 && n.compare(n.size() - 4, 4, "2000") == 0) { _keep.push_back(i); _keepName.push_back("2000"); break; }
      }
      std::fprintf(_f, "ievt");
      std::fprintf(_f, ",w:nominal");
      for (const std::string& n : _keepName) std::fprintf(_f, ",w:%s", n.c_str());
      std::fprintf(_f, ",ptz,yz,etaz,phiz,mz,pth,yh,etah,phih,mzh,ptzh,yzh,drzh,njets30,ptj1,njets30eta47,ptj1eta47,cell\n");
      MSG_INFO("weights: " << names.size() << " in the event, " << _keep.size() << " written (nominal index " << inom << ")");
      for (size_t k = 0; k < _keepName.size(); ++k) MSG_INFO("  w:" << _keepName[k] << " = " << names[_keep[k + 1]]);
      _header = true;
    }

    void analyze(const Event& event) {
      if (!_header) header();
      const std::valarray<double> w = event.weights();
      std::fprintf(_f, "%zu", _n++);
      for (size_t i : _keep) std::fprintf(_f, ",%.8g", i < w.size() ? w[i] : std::numeric_limits<double>::quiet_NaN());

      Particle Z, H;
      bool hasZ = false, hasH = false;
      for (HepMC3::ConstGenParticlePtr gp : HepMCUtils::particles(event.genEvent())) {
        if (!hasZ && gp->pdg_id() == 23) { Z = lastCopy(Particle(gp)); hasZ = true; }
        if (!hasH && gp->pdg_id() == 25) { H = lastCopy(Particle(gp)); hasH = true; }
        if (hasZ && hasH) break;
      }
      if (!hasZ || !hasH) {
        ++_nMissing;
        std::fprintf(_f, ",nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,-1,-1,-1,-1,-1\n");
        return;
      }
      Particles hadrons;
      for (const Particle& p : apply<FinalState>(event, "FS").particles()) {
        if (descendsFrom(p, H) || descendsFrom(p, Z)) continue;
        hadrons += p;
      }
      FinalState fsTemp;
      FastJets fj(fsTemp, JetAlg::ANTIKT, 0.4);
      fj.calc(hadrons);
      const Jets jets = fj.jetsByPt(Cuts::pT > 30.0);
      int n47 = 0; double pt47 = -1;
      for (const Jet& j : jets) if (j.abseta() < 4.7) { if (n47 == 0) pt47 = j.pT(); ++n47; }

      const FourMomentum pz = Z.momentum(), ph = H.momentum(), pzh = pz + ph;
      const double drzh = deltaR(pz.eta(), pz.phi(), ph.eta(), ph.phi());
      int cell = -1;
      if (ph.absrap() <= 2.5) {
        const double pt = pz.pT();
        const int ptbin = pt < 75 ? 0 : pt < 150 ? 1 : pt < 250 ? 2 : pt < 400 ? 3 : pt < 600 ? 4 : 5;
        cell = 3 * ptbin + std::min<int>(jets.size(), 2);
      }
      std::fprintf(_f, ",%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%zu,%.6g,%d,%.6g,%d\n",
                   pz.pT(), pz.rapidity(), pz.eta(), pz.phi(), pz.mass(), ph.pT(), ph.rapidity(), ph.eta(), ph.phi(),
                   pzh.mass(), pzh.pT(), pzh.rapidity(), drzh, jets.size(), jets.empty() ? -1.0 : jets[0].pT(), n47, pt47, cell);
    }

    void finalize() {
      if (_nMissing) MSG_WARNING(_nMissing << " events without a Z or Higgs boson in the record");
      // Rivet 4 calls finalize once per weight stream, after all events: close the file the first time
      if (_f) { std::fclose(_f); _f = nullptr; MSG_INFO(_n << " events written"); }
    }

  private:
    std::FILE* _f = nullptr;
    std::vector<size_t> _keep;
    std::vector<std::string> _keepName;
    bool _header = false;
    size_t _n = 0, _nMissing = 0;
  };

  RIVET_DECLARE_PLUGIN(GGZH_EVENTTABLE);

}
