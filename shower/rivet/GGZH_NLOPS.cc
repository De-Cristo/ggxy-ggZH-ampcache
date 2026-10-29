// -*- C++ -*-
// GGZH_NLOPS: distributions for gg -> Z(l+l-)H at NLO+PS / LO+PS, companion of the LHC Higgs WG
// routine SimplifiedTemplateCrossSections (STXS categories).
//
// Z and H are the last copies of the bosons in the event record.  Jets: anti-kT R = 0.4 on all
// stable particles except the descendants of the Higgs boson and of the Z boson (as the STXS
// routine does for leptonic V decays), p_T > 30 GeV, no rapidity cut.
// Normalisation: sum of weights over the number of events (POWHEG's convention), for every
// weight of the event (Rivet multi-weights), so the histograms are cross sections per bin.
// "stxs13_fine_own": the 18 Stage 1.3 fine cells computed here (p_T(Z) [0,75,150,250,400,600,inf)
// x 0J/1J/>=2J, index 3*ptbin + jetbin; 18 = |y_H| > 2.5), to cross-check the STXS routine.
// "sigma": one bin with every event, i.e. the inclusive cross section.
#include "Rivet/Analysis.hh"
#include <cstdio>
#include <algorithm>
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  class GGZH_NLOPS : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(GGZH_NLOPS);

    void init() {
      declare(FinalState(), "FS");
      const std::string check = getOption("CHECKOUT", "");
      if (!check.empty()) {
        _check = std::fopen(check.c_str(), "w");
        if (!_check) throw Error("cannot open GGZH_NLOPS CHECKOUT");
        std::fprintf(_check, "ptz,pth,mzh,ptzh,njets30,ptj1,etaj1,phij1,ptj2,etaj2,phij2,ptl1,etal1,phil1,ptl2,etal2,phil2\n");
      }
      auto edges = [](double lo, double hi, double step) {
        std::vector<double> e;
        for (double x = lo; x <= hi + 1e-9; x += step) e.push_back(x);
        return e;
      };
      book(_h["ptz"], "ptz", edges(0., 1000., 25.));
      book(_h["ptzwide"], "ptzwide", edges(0., 4000., 100.));
      book(_h["pth"], "pth", edges(0., 1000., 25.));
      book(_h["mzh"], "mzh", edges(150., 1250., 25.));
      book(_h["ptzh"], "ptzh", edges(0., 300., 10.));
      book(_h["yzh"], "yzh", 20, -3., 3.);
      book(_h["njets"], "njets", std::vector<double>{-0.5, 0.5, 1.5, 2.5, 3.5});
      book(_h["ptj1"], "ptj1", edges(0., 1000., 25.));
      book(_h["stxs13_fine_own"], "stxs13_fine_own", 19, -0.5, 18.5);
      book(_h["sigma"], "sigma", 1, 0., 1.);
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

    void analyze(const Event& event) {
      _h["sigma"]->fill(0.5);
      Particle Z, H;
      bool hasZ = false, hasH = false;
      for (HepMC3::ConstGenParticlePtr gp : HepMCUtils::particles(event.genEvent())) {
        if (!hasZ && gp->pdg_id() == 23) { Z = lastCopy(Particle(gp)); hasZ = true; }
        if (!hasH && gp->pdg_id() == 25) { H = lastCopy(Particle(gp)); hasH = true; }
        if (hasZ && hasH) break;
      }
      if (!hasZ || !hasH) { ++_nMissing; return; }

      Particles hadrons;
      for (const Particle& p : apply<FinalState>(event, "FS").particles()) {
        if (descendsFrom(p, H) || descendsFrom(p, Z)) continue;
        hadrons += p;
      }
      FinalState fsTemp;
      FastJets fj(fsTemp, JetAlg::ANTIKT, 0.4);
      fj.calc(hadrons);
      const Jets jets = fj.jetsByPt(Cuts::pT > 30.0);

      const FourMomentum pz = Z.momentum(), ph = H.momentum(), pzh = pz + ph;
      if (_check) {
        Jets accepted;
        for (const Jet& j : jets) if (j.abseta() < 4.7) accepted.push_back(j);
        std::fprintf(_check, "%.12g,%.12g,%.12g,%.12g,%zu,%.12g", pz.pT(), ph.pT(), pzh.mass(), pzh.pT(),
                     jets.size(), accepted.empty() ? -1.0 : accepted[0].pT());
        for (size_t k = 0; k < 2; ++k) {
          if (k) std::fprintf(_check, ",%.12g", accepted.size() > k ? accepted[k].pT() : -1.0);
          std::fprintf(_check, ",%.12g,%.12g", accepted.size() > k ? accepted[k].eta() : -1.0,
                       accepted.size() > k ? accepted[k].phi() : -1.0);
        }
        Particles children;
        for (const Particle& child : Z.children())
          if (child.abspid() == 11 || child.abspid() == 13 || child.abspid() == 15) children.push_back(child);
        std::sort(children.begin(), children.end(), [](const Particle& a, const Particle& b) { return a.pT() > b.pT(); });
        for (size_t k = 0; k < 2; ++k)
          std::fprintf(_check, ",%.12g,%.12g,%.12g", children.size() > k ? children[k].pT() : -1.0,
                       children.size() > k ? children[k].eta() : -1.0, children.size() > k ? children[k].phi() : -1.0);
        std::fprintf(_check, "\n");
      }
      _h["ptz"]->fill(pz.pT());
      _h["ptzwide"]->fill(pz.pT());
      _h["pth"]->fill(ph.pT());
      _h["mzh"]->fill(pzh.mass());
      _h["ptzh"]->fill(pzh.pT());
      _h["yzh"]->fill(pzh.rapidity());
      _h["njets"]->fill(std::min<size_t>(jets.size(), 3));
      if (!jets.empty()) _h["ptj1"]->fill(jets[0].pT());

      int cell = 18;
      if (ph.absrap() <= 2.5) {
        const double pt = pz.pT();
        const int ptbin = pt < 75 ? 0 : pt < 150 ? 1 : pt < 250 ? 2 : pt < 400 ? 3 : pt < 600 ? 4 : 5;
        const int jetbin = std::min<int>(jets.size(), 2);
        cell = 3 * ptbin + jetbin;
      }
      _h["stxs13_fine_own"]->fill(cell);
    }

    void finalize() {
      if (_check) { std::fclose(_check); _check = nullptr; }
      if (_nMissing) MSG_WARNING(_nMissing << " events without a Z or Higgs boson in the record");
      const double sf = 1.0 / numEvents();
      for (auto& h : _h) scale(h.second, sf);
    }

  private:
    std::map<std::string, Histo1DPtr> _h;
    size_t _nMissing = 0;
    std::FILE* _check = nullptr;
  };

  RIVET_DECLARE_PLUGIN(GGZH_NLOPS);

}
