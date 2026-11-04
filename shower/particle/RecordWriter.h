// Passive record writer. No random draws, no mutation of Pythia event/process.
#pragma once
#include "Pythia8/LesHouches.h"
#include "fastjet/ClusterSequence.hh"
#include <fstream>
#include <iomanip>
#include <set>
#include <cmath>
using namespace Pythia8;
struct CountingLHEF : LHAupLHEF {
  long ordinal = -1;
  CountingLHEF(Info* info, const char* file): LHAupLHEF(info,file,nullptr,true,false) {}
  bool setEvent(int id=0) override {
    bool ok=LHAupLHEF::setEvent(id); if(ok) ++ordinal; return ok;
  }
};
class RecordWriter {
  std::ofstream o;
  void ints(const std::vector<int>& v) {
    o << '['; for(size_t k=0;k<v.size();++k) {if(k)o<<',';o<<v[k];} o<<']';
  }
  void record(const Event& ev) {
    o << '[';
    for(int i=0;i<ev.size();++i) {
      if(i)o<<',';const auto& p=ev[i];
      o << '[' << p.id()<<','<<p.status()<<','<<p.px()<<','<<p.py()<<','<<p.pz()<<','<<p.e()<<','<<p.m()<<','
        <<p.col()<<','<<p.acol()<<','<<(p.isFinal()?1:0)<<',';
      ints(p.motherList());o<<',';ints(p.daughterList());o<<']';
    } o<<']';
  }
  bool desc(const Event& e,int i,int a) {
    if(i==a)return true;
    std::vector<int> todo=e[i].motherList();std::set<int> seen;
    while(!todo.empty()) {
      int m=todo.back();todo.pop_back();if(m==a)return true;
      if(m<=0 || m>=e.size() || !seen.insert(m).second)continue;
      auto ms=e[m].motherList();todo.insert(todo.end(),ms.begin(),ms.end());
    }return false;
  }
  void jets(const std::vector<fastjet::PseudoJet>& particles,fastjet::JetAlgorithm algo,double R) {
    o<<"{\"jets\":[";
    if(!particles.empty()) {
      fastjet::ClusterSequence cs(particles,fastjet::JetDefinition(algo,R));
      auto js=fastjet::sorted_by_pt(cs.inclusive_jets(10.));
      for(size_t j=0;j<js.size();++j) {
        if(j)o<<',';auto& p=js[j];o<<'['<<p.px()<<','<<p.py()<<','<<p.pz()<<','<<p.e()<<",[";
        auto cc=p.constituents();for(size_t k=0;k<cc.size();++k) {if(k)o<<',';o<<cc[k].user_index();}o<<"]]";
      }o<<"],\"dmerge\":[";
      for(int k=0;k<8;++k) {if(k)o<<',';o<<(k<int(particles.size())?cs.exclusive_dmerge(k):-1.);}
    }else o<<"],\"dmerge\":[-1,-1,-1,-1,-1,-1,-1,-1";
    o<<"]}";
  }
public:
  RecordWriter(const std::string& path):o(path) {o<<std::setprecision(17);if(!o)throw std::runtime_error("record output");}
  void write(Pythia& p,long ordinal,long accepted) {
    const auto& e=p.event;int z=-1,h=-1;
    for(int i=1;i<e.size();++i) {
      if(e[i].id()!=23 && e[i].id()!=25)continue;
      bool last=true;for(int d:e[i].daughterList())if(d>0 && d<e.size() && e[d].id()==e[i].id())last=false;
      if(last) {if(e[i].id()==23 && z<0)z=i; if(e[i].id()==25 && h<0)h=i;}
    }
    std::vector<fastjet::PseudoJet> stable,partons;
    std::vector<int> ref;
    for(int i=1;i<e.size();++i) {
      if((z>=0&&desc(e,i,z)) || (h>=0&&desc(e,i,h)))continue;
      const auto& q=e[i];fastjet::PseudoJet v(q.px(),q.py(),q.pz(),q.e());v.set_user_index(i);
      if(q.isFinal())stable.push_back(v);
      int id=std::abs(q.id()),st=std::abs(q.status());
      if((id<=6 || id==21)&&st>=71&&st<=79&&q.pT()>=10.) {
        ref.push_back(i);bool copied=false;
        for(int d:q.daughterList())if(d>0&&d<e.size()) {
          int di=std::abs(e[d].id()),ds=std::abs(e[d].status());
          if((di<=6||di==21)&&ds>=71&&ds<=79)copied=true;
        }
        if(!copied)partons.push_back(v);
      }
    }
    o<<"{\"lhe_ordinal\":"<<ordinal<<",\"accepted\":"<<accepted<<",\"z\":"<<z<<",\"h\":"<<h
     <<",\"pythia_incoming\":["<<p.info.x1()<<','<<p.info.x2()<<','<<p.info.id1()<<','<<p.info.id2()<<']'
     <<",\"event\":";record(e);o<<",\"process\":";record(p.process);
    o<<",\"radiation\":";jets(stable,fastjet::antikt_algorithm,.4);
    o<<",\"partons\":";jets(partons,fastjet::kt_algorithm,1.);
    o<<",\"parton_indices\":[";for(size_t k=0;k<partons.size();++k){if(k)o<<',';o<<partons[k].user_index();}
    o<<"],\"status71_79_indices\":";ints(ref);
    // Jet replacement uses anti-kT radiation jets as kT clustering inputs.
    std::vector<fastjet::PseudoJet> replacement;
    if(!stable.empty()) {
      fastjet::ClusterSequence cs(stable,fastjet::JetDefinition(fastjet::antikt_algorithm,.4));
      replacement=fastjet::sorted_by_pt(cs.inclusive_jets(10.));
      for(size_t j=0;j<replacement.size();++j){auto v=replacement[j];replacement[j]=fastjet::PseudoJet(v.px(),v.py(),v.pz(),v.e());replacement[j].set_user_index(j);}
    }
    o<<",\"replacement_kt\":";jets(replacement,fastjet::kt_algorithm,1.);o<<"}\n";
    if(!o)throw std::runtime_error("particle writer failed");
  }
};
