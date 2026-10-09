// Fixed-order LO / NLO gg -> l- l+ H (one lepton generation) with the ggxy Monte Carlo.
// Derived from ggxy-release examples/ggzh-nlo/nlo-ggllh.cpp (commit 9edd1d5), the program
// behind arXiv:2603.15762, Table 3 and Fig. 6.  Changes: PDF set, number of points and seed
// from the command line (default NNPDF3.1, LHAPDF 325300), and five extra histograms
// (pT(Z), pT(H) and m_ZH in 25 GeV bins, pT(ZH) in 10 GeV bins, wide pT(Z)).  Physics
// settings unchanged: 13.6 TeV, mu_R = mu_F = m_ZH/2, on-shell m_t = 172.5 GeV,
// 15 < m_ll < 150 GeV, Z radiation off the external quark lines on.
//   fo-ggllh <lo|nlo> <seed> <outdir/> [np] [np_opt] [pdf] [sqrt(s) in GeV, default 13600]
// Histograms are written in nb per unit of the observable, for ONE lepton flavour.
#include <iostream>
#include <iomanip>
#include <string>
#include <chrono>
#include <complex>
#include <algorithm>

#include "cs_dipoles/cs_dipoles.h"
#include "tools/lorentz_vec.h"
#include "tools/params.h"
#include "tools/histo.h"
#include "tools/mc_gen.h"
#include "mc_config/mc_ggllh.h"


//function to define renormalization, factorization and MSbar scales
void set_scale(const std::vector<int>& iproc, const std::vector<lorentz_vec>& p, const params& pars, std::vector<double>& muR, std::vector<double>& muF, double& mut){
   lorentz_vec pzh(0.0,0.0,0.0,0.0);
   for(unsigned i=0;i<iproc.size();i++){
      if(std::abs(iproc[i]) == 11 or std::abs(iproc[i])==12 or iproc[i] == 25) pzh+=p[i];
   }
   muR[0]=pzh.mass()/2.0;
   muF[0]=muR[0];

   muR[1]=pars.mtscheme == 1 ? pars.mtmt : pars.mt;
   muF[1]=muR[1];

   //set MSbar scale
   if(pars.mtscheme == 1) mut = pars.mtmt; //mt(mt)
}


//phase-space cuts
bool event_cut(const std::vector<int>& iproc, const std::vector<lorentz_vec>& plab, const params& pars){
   //can be used to impose cuts on the final states

   //e.g. a jet-veto
   //jet only present in real corrections pj=-phh
   //const double ptveto = 50.0;
   //if(ptj > ptveto) return false;

   const double mllmin =  15.0;
   const double mllmax = 150.0;

   lorentz_vec pz(0.0,0.0,0.0,0.0);
   for(unsigned i=0;i<iproc.size();i++){
      if(std::abs(iproc[i]) == 11 or std::abs(iproc[i]) == 12) pz+=plab[i];
   }

   const double mll = pz.mass();
   if(mll < mllmin or mll > mllmax) return false;

   return true;
}


//function to init histograms
void init_histograms(histo& hlist){
   hlist.add("cross",1,0.0,1.0);      //total cross section
   hlist.add("ptz",20,0.0,600.0);     //pT of Z boson
   hlist.add("pth",20,0.0,600.0);     //pT of Higgs boson
   hlist.add("ptzh",20,0.0,600.0);    //pT of ZH system
   hlist.add("mzh",20,250.0,1000.0);  //mass of HH system
   hlist.add("yz",20,-3.0,3.0);       //rapidity of Z boson
   hlist.add("yh",20,-3.0,3.0);       //rapidity of Higgs boson
   hlist.add("yzh",20,-3.0,3.0);      //rapitidy of ZH system

   //histogram of decay products
   hlist.add("ptlm",20,0.0,400.0);        //pt of negative charged lepton
   hlist.add("ptlp",20,0.0,400.0);        //pt of positive charged lepton
   hlist.add("mll",40,40.0,140.0);        //mass of ll(Z) system
   hlist.add("ylm",20,-2.5,2.5);          //rapidity of negative charged lepton
   hlist.add("ylp",20,-2.5,2.5);          //rapidity of positive charged lepton
   hlist.add("phill",20,0.0,M_PI);        //phi angle between both charged leptons
   hlist.add("cosll",20,-1.0,1.0);        //cos\theta between both charged leptons
   hlist.add("rll",20,0.0,4.0);           //\Delta R between both charged leptons
   hlist.add("mll2",40,15.0,150.0);       //mass of ll(Z) system

   //additional binnings (25 GeV bins contain the STXS pT(V) boundaries 75, 150, 250, 400, 600 GeV)
   hlist.add("ptz25",40,0.0,1000.0);      //17 pT of Z boson
   hlist.add("pth25",40,0.0,1000.0);      //18 pT of Higgs boson
   hlist.add("ptzh10",30,0.0,300.0);      //19 pT of ZH system
   hlist.add("mzh25",44,150.0,1250.0);    //20 mass of ZH system
   hlist.add("ptzwide",40,0.0,4000.0);    //21 pT of Z boson, full range (100 GeV bins)
}


//function to fill histograms
void fill_histograms(histo& hlist, const std::vector<int>& iproc, const std::vector<lorentz_vec>& plab, const params& pars){
   lorentz_vec ph,plm,plp;
   double ptz,pth;

   lorentz_vec pz = lorentz_vec(0.0,0.0,0.0,0.0);

   bool ilepton = false;
   for(unsigned i=0;i<iproc.size();i++){
      if(std::abs(iproc[i]) == 11 or std::abs(iproc[i]) == 12) pz += plab[i];
      //if(iproc[i] == 11 or iproc[i] == 12){
      if(iproc[i] == 11){
         ilepton = true;
         plm = plab[i];
      //} else if(iproc[i] == -11 or iproc[i] == -12){
      } else if(iproc[i] == -11){
         ilepton = true;
         plp = plab[i];
      } else if(iproc[i] == 25){
         ph = plab[i];
      }
   }

   lorentz_vec pzh = pz + ph;

   //further observables can easily be added, see src/tools/lorentz_vec.h for more built-in functions of the lorentz_vec class

   hlist.fill( 0,0.5);                     //total cross section
   hlist.fill( 1,pz.pt());                 //pT of Z boson
   hlist.fill( 2,ph.pt());                 //pT of Higgs boson
   hlist.fill( 3,pzh.pt());                //pT of ZH system
   hlist.fill( 4,pzh.mass());              //mass of ZH system
   hlist.fill( 5,pz.rap());                //rapidity of Z boson
   hlist.fill( 6,ph.rap());                //rapidity of soft Higgs
   hlist.fill( 7,pzh.rap());               //rapitidy of HH system
   hlist.fill(17,pz.pt());
   hlist.fill(18,ph.pt());
   hlist.fill(19,pzh.pt());
   hlist.fill(20,pzh.mass());
   hlist.fill(21,pz.pt());

   if(not ilepton) return;

   //extra observables for charged leptons
   hlist.fill( 8,plm.pt());                //pt of negative charged lepton
   hlist.fill( 9,plp.pt());                //pt of positive charged lepton
   hlist.fill(10,(plm+plp).mass());        //mass of ll(Z) system
   hlist.fill(11,plm.rap());               //rapidity of negative charged lepton
   hlist.fill(12,plp.rap());               //rapidity of positive charged lepton
   hlist.fill(13,get_phi(plm,plp));        //phi angle between both charged leptons
   hlist.fill(14,get_cos(plm,plp));        //cos\theta between both charged leptons
   hlist.fill(15,get_deltaR(plm,plp));     //\Delta R between both charged leptons
   hlist.fill(16,(plm+plp).mass());        //mass of ll(Z) system
}



int main(int argc, char* argv[])
{
   std::cout << std::setprecision(17);
   std::cerr << std::setprecision(17);

   int seed=1001;                                   //random number seed

   //setup: https://arxiv.org/pdf/2508.09905
   double ss=13600.0;                               //center-of-mass energy
   double mt=172.5;                                 //top quark mass
   double mz=91.1876;                               //Z boson mass
   double wz=2.4955;                                //Z boson decay width
   double mh=125.0;                                 //Higgs boson mass
   double mw=80.385;                                //W boson mass
   double gfermi=1.16637e-5;                        //fermi constant
   double kappa_hzz=1.0;                            //rescaling of HZZ vertex
   std::string pdf_name="NNPDF31_nnlo_as_0118_mc_hessian_pdfas";  //PDF set (LHAPDF 325300)
   unsigned mtscheme=0;                             //top-quark mass scheme: 0 - on-shell, 1 - MSbar
   int iopt_ext_Z = 1;                              //0: deactivate, 1: activate Z radiation from external quark lines
   std::string decaymode = "ll";                    //ll: l-l+, vv: vl vl~

   int np    = 1000000;         //number of generated phase-space points for MC integration
   int np_opt=  150000;         //number of generated phase-space points for phase-space optimization


   //mode = "lo"     -> run  LO calculation
   //mode = "nlo"    -> run NLO calculation
   //the real subtracted cross section and the Born like contributions can also be generated separately
   //mode = "RS"     -> run real subtracted part (gg->gHH + qg->qHH + gq->qHH + qqx->gHH)
   //mode = "LOVDIP" -> run Born+Virtual+integrated dipoles
   std::string mode="nlo";
   std::string path="";
   if(argc >= 2) mode=std::string(argv[1]);
   if(argc >= 3) seed=atoi(argv[2]);
   if(argc >= 4) path=std::string(argv[3]);
   if(argc >= 5) np=atoi(argv[4]);
   if(argc >= 6) np_opt=atoi(argv[5]);
   if(argc >= 7) pdf_name=std::string(argv[6]);
   if(argc >= 8) ss=atof(argv[7]);                  //centre-of-mass energy [GeV]
   std::cout<<"fo-ggllh: mode "<<mode<<" seed "<<seed<<" np "<<np<<" np_opt "<<np_opt<<" pdf "<<pdf_name<<" sqrt(s) "<<ss<<std::endl;

   std::string outfile=path+"histo_"+mode+"_"+std::to_string(seed)+".txt";  //output file for histograms

   //init mc generator
   mc_gen gen=mc_gen(seed,ss);
   configure_mc_ggllh(mode,gen,decaymode);

   //set parameters
   gen.set_mass_top(mt);
   gen.set_mass_zboson(mz,wz);
   gen.set_mass_higgs(mh);
   gen.set_mass_wboson(mw);
   gen.set_gfermi(gfermi);
   gen.set_kappa_hzz(kappa_hzz);
   gen.set_mtscheme(mtscheme);
   gen.set_iopt_ext_Z(iopt_ext_Z);


   //init names of scales, gridtype, PDF set
   gen.set_scales_pdfs({"mzh2","mt"},1,pdf_name);
   gen.set_scalefunc(set_scale);


   //set functions to init and fill histograms
   gen.set_hinit(init_histograms);
   gen.set_hfill(fill_histograms);


   //phase-space cuts
   gen.set_event_cut(event_cut);


   //initialization finished, no parameters should be changed afterwards
   gen.finish_init();


   std::chrono::time_point<std::chrono::steady_clock> start;
   std::chrono::duration<double> timeEval1(0);
   start = std::chrono::steady_clock::now();

   //phase-space optimization
   gen.integrate(np_opt);
   gen.reset_weights();
   gen.set_phase_optim(false);
   //integration
   std::vector<double> res=gen.integrate(np);
   gen.output(outfile);



   std::cout<<mode<<std::endl;

   std::cout<<"res_ggxy = "<<res[0]<<" +- "<<res[1]<<std::endl;
   timeEval1 += (std::chrono::steady_clock::now() - start);
   std::cout<<"Eval time: "<<timeEval1.count()<< std::endl;


  return 0;
}

//#####################################################################

