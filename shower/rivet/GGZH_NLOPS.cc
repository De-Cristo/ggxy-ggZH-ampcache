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
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  class GGZH_NLOPS : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(GGZH_NLOPS);

    void init() {
      declare(FinalState(), "FS");
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
      if (_nMissing) MSG_WARNING(_nMissing << " events without a Z or Higgs boson in the record");
      const double sf = 1.0 / numEvents();
      for (auto& h : _h) scale(h.second, sf);
    }

  private:
    std::map<std::string, Histo1DPtr> _h;
    size_t _nMissing = 0;
  };

  RIVET_DECLARE_PLUGIN(GGZH_NLOPS);

}
