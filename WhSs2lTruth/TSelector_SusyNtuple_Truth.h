#ifndef SusyAna_TSelector_SusyNtuple_Truth_h
#define SusyAna_TSelector_SusyNtuple_Truth_h

//////////////////////////////////////////////////////////
// General script to implement basic selection with all //
// signal region cut methods.                           //
//////////////////////////////////////////////////////////


#include "SusyNtuple/SusyNtTruthAna.h"
#include "SusyNtuple/SusyNtTools.h"
#include "SusyNtuple/SusyNtAna.h"
#include "SusyNtuple/DilTrigLogic.h"
#include "SusyNtuple/SusyDefs.h"
#include "SusyNtuple/WhTruthExtractor.h"

// Root Packages
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include <fstream>


using Susy::Lepton;
using Susy::Jet;
using Susy::Met;
struct FourMom {
    double px, py, pz, E;
    int q;
    bool isMu, isEl, isJet;
    FourMom() : px(0), py(0), pz(0), E(0), q(0), isMu(false), isEl(false), isJet(false) {}
#ifndef __CINT__
// cint is not able to parse 'complex' code; see
// http://root.cern.ch/drupal/content/interacting-shared-libraries-rootcint
    FourMom& set4mom(const Lepton &l) { px=l.Px(); py=l.Py(); pz=l.Pz(); E=l.E(); q=l.q; return *this; }
    FourMom& set4mom(const Jet &j)    { px=j.Px(); py=j.Py(); pz=j.Pz(); E=j.E(); return *this; }
    FourMom& setMu(const Lepton &l) { isMu=true; isEl = isJet = false; return set4mom(l); }
    FourMom& setEl(const Lepton &l) { isEl=true; isMu = isJet = false; return set4mom(l); }
    FourMom& setJet(const Jet &j)   { isJet=true; isMu = isEl = false; return set4mom(j); }
    FourMom& setMet(const Met &m)   { isJet=isMu=isEl=false; px=m.lv().Px(); py=m.lv().Py(); E=m.lv().E(); return *this; }
#endif
};

struct EventParameters {
    double weight;
    unsigned int eventNumber;
    unsigned int runNumber;
    bool ismc;
    EventParameters() : weight(0), eventNumber(0), runNumber(0), ismc(0) {}
#ifndef __CINT__
    EventParameters& setWeight(const double &w) { weight=w; return *this; }
    EventParameters& setEvent(const unsigned int &e) { eventNumber=e; return *this; }
    EventParameters& setRun(const unsigned int &r) { runNumber=r; return *this; }
    EventParameters& setIsmc(const bool &im) { ismc=im; return *this; }
#endif
};

class TSelector_SusyNtuple_Truth : public SusyNtTruthAna
{

  public:
    
//   Event* event_test;
//   EventParameters   pars; 
//   TBranch        *b_pars; 
  TBranch* m_eventParameters_b;
  EventParameters* m_eventParameters;
  
//   int truthParticles_pdgId;
//   TBranch *truthParticles_pdgId;
  
  unsigned int eventNumber;
  TBranch *b_pars_eventNumber;    
    
  TSelector_SusyNtuple_Truth();
  virtual ~TSelector_SusyNtuple_Truth(){};
  
  // Output Text File
  ofstream out;
  // Begin is called before looping on entries
  virtual void    Begin(TTree *tree);
  virtual void    SlaveBegin(TTree *tree);
  virtual Bool_t  Process(Long64_t entry);

  virtual void    Terminate();
  virtual void    SlaveTerminate();
		     
    // Cut methods
    /* bool CheckRealLeptons(const ElectronVector& signal_electrons, MuonVector& signal_muons); */
    /* bool CheckChargeFlipElectrons(const ElectronVector& signal_electrons); */
    float getBTagWeight(const Event* evt, BTagSys SysSettingBTag);
    float recalcMetRel(TLorentzVector metLV, TLorentzVector l1, TLorentzVector l2, const TruthJetVector jets, bool useForward);
    float calcMT2(TLorentzVector metlv, TLorentzVector l0, TLorentzVector l1);
    float calcMT2J(TLorentzVector metlv, TLorentzVector l0, TLorentzVector l1, TLorentzVector j0, TLorentzVector j1);
    bool checkLeptonPt(const LeptonVector& leptons);
    float mTWW(TLorentzVector _ll, TLorentzVector _nu, bool MvvTrue);
    float calcHT(TLorentzVector l1, TLorentzVector l2, TLorentzVector met, TruthJetVector &signalJets);
    float calcMeff(TLorentzVector met, const JetVector &signalJets);
    float calcMZTauTau_coll(const TLorentzVector &signal_lep_0, const TLorentzVector &signal_lep_1, const TLorentzVector &met);
    float calcMZTauTau_mmc(TLorentzVector lep1, TLorentzVector lep2, int tau0_decay_type, int tau1_decay_type);
    float calcMt(TLorentzVector _l, TLorentzVector _nu);
    float calcSumMv1(const JetVector &signalJets);
    
    void defineHistos();
    void defineHistos_sysUncert();
    void writeHistos();
    void writeHistos_sysUncert();
//     void addHistos();
    void fillHistos_EE(int cutnumber, float weight);
    void fillHistos_MM(int cutnumber, float weight);
    void fillHistos_EM(int cutnumber, float weight);
    
    void calcJet_variables(TLorentzVector met_TLV, TruthJetVector signalJets);
    void calc_EE_variables(TruthParticle* el0, TruthParticle* el1, TLorentzVector el0_TLV, TLorentzVector el1_TLV, TruthJetVector signalJets, TLorentzVector signalJet0_TLV, TLorentzVector signalJet1_TLV, TLorentzVector met_TLV);
    void calc_MM_variables(TruthParticle* mu0, TruthParticle* mu1, TLorentzVector mu0_TLV, TLorentzVector mu1_TLV, TruthJetVector signalJets, TLorentzVector signalJet0_TLV, TLorentzVector signalJet1_TLV, TLorentzVector met_TLV);
    void calc_EM_variables(TruthParticle* mu, TruthParticle* el, TLorentzVector mu_TLV, TLorentzVector el_TLV, TruthJetVector signalJets, TLorentzVector signalJet0_TLV, TLorentzVector signalJet1_TLV, TLorentzVector met_TLV);
    
    void fillHistos_EE_SRSS1(float cut_EE, float weight_ALL_EE);    
    void fillHistos_MM_SRSS1(float cut_MM, float weight_ALL_MM);    
    void fillHistos_EM_SRSS1(float cut_EM, float weight_ALL_EM);
    
    
    float calc_D0(bool unbiased, const Lepton* lep);
    static bool compareElecMomentum (Electron* e0, Electron* e1);
    static bool compareMuonMomentum (Muon* mu0, Muon* mu1);
//     ElectronVector getSoftElectrons(SusyNtObject* susyNt, SusyNtSys sys);
//     ElectronVector getOverlapElectrons(SusyNtObject* susyNt, SusyNtSys sys);
//     MuonVector getSoftMuons(SusyNtObject* susyNt, SusyNtSys sys);
//     MuonVector getOverlapMuons(SusyNtObject* susyNt, SusyNtSys sys);
    bool isCMSJet(const Susy::Jet* jet);
    int numberOfCMSJets(const JetVector& jets);    
    vector<TLorentzVector> overlapRemoval(vector<TLorentzVector> objects_type1, vector<TLorentzVector> indices_2, double dr, bool sameType, bool removeSoft) ;
    bool doEventCleaning_andFillHistos(TruthJetVector baseJets, const TruthMet* met, TruthParticleVector baseMuons, TruthParticleVector baseElectrons, int flag, float weight_ALL_EE, float weight_ALL_MM, float weight_ALL_EM, SusyNtSys SysSetting, bool n0150BugFix);
    bool initTupleMaker(const std::string &outFilename, const std::string &treename);
    bool initFile(const std::string &outFilename);
    bool initTree(const std::string &treename);
    bool initTreeBranches();
    bool fillTupleMaker(const double weight, const unsigned int run, const unsigned int event, const bool isMc, const Susy::Lepton &l0, const Susy::Lepton &l1, const Susy::Met &met, const LeptonVector &otherLeptons, const JetVector &jets);
    LeptonVector getAnyElOrMu(SusyNtObject &susyNt, const Lepton *l0, const Lepton *l1);
    bool closeTupleMaker();
//     SusyNtSys returnSysUncertType(int isys);
//     int returnSysUncertNumber(SusyNtSys UncertType);
//     int returnSysUncertBTagNumber(BTagSys UncertType);
    char *GetID(Int_t type);
    float getLeptonSF(const LeptonVector& leptons, int SFUncertType);
    void overlapProcedure(TruthParticleVector& elecs, TruthParticleVector& muons, TruthParticleVector& taus, TruthJetVector& jets);
    void e_e_overlap(TruthParticleVector& elecs, float minDr);
  
    // e-j overlap
    void e_j_overlap(TruthParticleVector& elecs, TruthJetVector& jets, float minDr, 
                     bool removeJets = true);
  
    // m-j overlap
    void m_j_overlap(TruthParticleVector& muons, TruthJetVector jets, float minDr);

    // e-m overlap 
    void e_m_overlap(TruthParticleVector& elecs, TruthParticleVector& muons, float minDr);
  
    // m-m overlap
    void m_m_overlap(TruthParticleVector& muons, float minDr);

    // t-e overlap
    void t_e_overlap(TruthParticleVector& taus, TruthParticleVector& elecs, float minDr);

    // t-m overlap
    void t_m_overlap(TruthParticleVector& taus, TruthParticleVector& muons, float minDr);

    // t-j overlap
    void t_j_overlap(TruthParticleVector& taus, TruthJetVector& jets, float minDr, 
                     bool removeJets = true);
    
    TruthJetVector getTruthJetsL20   (TruthJetVector& truthBaseJets);
    TruthJetVector getTruthJetsF30   (TruthJetVector& truthBaseJets);
    TruthJetVector getTruthJetsB20   (TruthJetVector& truthBaseJets);
    TruthJetVector getTruthJetsSignal   (TruthJetVector& truthBaseJets);

    // Selection region
    void setSelection(std::string s) { m_sel = s; }

    void print_counters() const;

    // debug check
    bool debugEvent();
    bool makeNTuple;
    bool run_on_SusyNtuple;
    bool calcSysUncert;

    int nSignalJets;
    float MET;
    float pTj0;
    float pTj1;
    float eta_j0;
    float eta_j1;
    float mjj;
    float DeltaPhijj;
    float pTjj;
    float DeltaPhiMETj0;
    float DeltaPhiMETj1;
    float DeltaRjj;
    float DeltaEtajj;
    float DeltaYjj;
    float DeltaPhiMETjj;
    int NBJets;
    int NCJets;
    int NFJets;
    float meff;
    
    float mZTT_coll;
    float mZTT_mmc; 
    
    //#####################################
    float pTl0_EE;
    float pTl1_EE;
    float etal0_EE;
    float etal1_EE;
    float DeltaR_EE;
    float pTll_EE;
    float Mll_EE;
    float METrel_EE;
    float MET_EE;
    float HT_EE;
    float mTWW_EE;
    float mT_EE;
    float mTmin_EE;
    float mTmax_EE;
    float mTl0MET_EE;
    float mTl1MET_EE;
    float mMET_EE;
    float DeltaPhi_EE;
    float DeltaPhiMETl0_EE;
    float DeltaPhiMETl1_EE;
    float DeltaPhiMETll_EE;
    float mT2_EE;  
    float mT2J_EE;       
    float DeltaPhilljj_EE;
    float DeltaPhil0jj_EE;
    float DeltaPhil1jj_EE;
    float DeltaRlljj_EE;
    float Mljj_EE;
    float Mlj_EE;
    float DeltaEtall_EE;
    
    float D0_branch_l0_EE;
    float D0_branch_l1_EE;
    float D0err_branch_l0_EE;
    float D0err_branch_l1_EE;
    
    float sD0Signif_branch_l0_EE;
    float sD0Signif_branch_l1_EE;
  
    int Nleptons_ZcandImpact_EE;
    float etcone30lZcandImpact_EE;
    
    float mllZcandImpact_EE;      
    float mTllZcandImpact_EE;
    float pTlZcandImpact_EE;
    float etalZcandImpact_EE;
    int IClZcandImpact_EE;
    float ptcone30lZcandImpact_EE;
    float d0SiglZcandImpact_EE;
    float z0SinThetalZcandImpact_EE;
    
    bool ZcandLep_exists_EE;
    bool ZcandLep_passesPT_EE;
    bool ZcandLep_passesEta_EE;
    bool ZcandLep_passesPTcone_EE;
    bool ZcandLep_passesETcone_EE;
    bool ZcandLep_passesD0_EE; 
    bool ZcandLep_passesZ0_EE; 
    bool ZcandLep_PassesMedium_EE;
    bool ZcandLep_PassesTight_EE; 
    bool ZcandLep_PassesORAndMllCut_EE;
    bool ZcandLep_PassesPR_EE;
    
    //#####################################
    
    float pTl0_MM;
    float pTl1_MM;
    float etal0_MM;
    float etal1_MM;
    float DeltaR_MM;
    float pTll_MM;
    float Mll_MM;
    float METrel_MM;
    float MET_MM;
    float HT_MM;
    float mTWW_MM;
    float mT_MM;
    float mTmin_MM;
    float mTmax_MM;
    float mTl0MET_MM;
    float mTl1MET_MM;
    float mMET_MM;
    float DeltaPhi_MM;
    float DeltaPhiMETl0_MM;
    float DeltaPhiMETl1_MM;
    float DeltaPhiMETll_MM;
    float mT2_MM;  
    float mT2J_MM;       
    float DeltaPhilljj_MM;
    float DeltaPhil0jj_MM;
    float DeltaPhil1jj_MM;
    float DeltaRlljj_MM;
    float Mljj_MM;
    float Mlj_MM;
    float DeltaEtall_MM;
    
    float D0_branch_l0_MM;
    float D0_branch_l1_MM;
    float D0err_branch_l0_MM;
    float D0err_branch_l1_MM;
    
    float sD0Signif_branch_l0_MM;
    float sD0Signif_branch_l1_MM;
    
    int Nleptons_ZcandImpact_MM;
    float mllZcandImpact_MM;      
    float mTllZcandImpact_MM;
    float pTlZcandImpact_MM;
    float etalZcandImpact_MM;
    int IClZcandImpact_MM;
    float ptcone30lZcandImpact_MM;
    float etcone30lZcandImpact_MM;
    float d0SiglZcandImpact_MM;
    float z0SinThetalZcandImpact_MM;
    
    bool ZcandLep_exists_MM;
    bool ZcandLep_passesPT_MM;
    bool ZcandLep_passesEta_MM;
    bool ZcandLep_passesPTcone_MM;
    bool ZcandLep_passesETcone_MM;
    bool ZcandLep_passesD0_MM; 
    bool ZcandLep_passesZ0_MM; 
    bool ZcandLep_PassesMedium_MM;
    bool ZcandLep_PassesTight_MM; 
    bool ZcandLep_PassesORAndMllCut_MM;
    bool ZcandLep_PassesPR_MM;
    //#####################################
    
    float pTl0_EM;
    float pTl1_EM;
    float etal0_EM;
    float etal1_EM;
    float DeltaR_EM;
    float pTll_EM;
    float Mll_EM;
    float METrel_EM;
    float MET_EM;
    float HT_EM;
    float mTWW_EM;
    float mT_EM;
    float mTmin_EM;
    float mTmax_EM;
    float mTl0MET_EM;
    float mTl1MET_EM;
    float mMET_EM;
    float DeltaPhi_EM;
    float DeltaPhiMETl0_EM;
    float DeltaPhiMETl1_EM;
    float DeltaPhiMETll_EM;
    float mT2_EM;  
    float mT2J_EM;       
    float DeltaPhilljj_EM;
    float DeltaPhil0jj_EM;
    float DeltaPhil1jj_EM;
    float DeltaRlljj_EM;
    float Mljj_EM;
    float Mlj_EM;
    float DeltaEtall_EM;
    
    float D0_branch_l0_EM;
    float D0_branch_l1_EM;
    float D0err_branch_l0_EM;
    float D0err_branch_l1_EM;
    
    float sD0Signif_branch_l0_EM;
    float sD0Signif_branch_l1_EM;
    
 
    int Nleptons_ZcandImpact_mu_EM;
    int Nleptons_ZcandImpact_el_EM;
    
    float mllZcandImpact_mu_EM;      
    float mTllZcandImpact_mu_EM;
    float pTlZcandImpact_mu_EM;
    float etalZcandImpact_mu_EM;
    int IClZcandImpact_mu_EM;
    float ptcone30lZcandImpact_mu_EM;
    float etcone30lZcandImpact_mu_EM;
    float d0SiglZcandImpact_mu_EM;
    float z0SinThetalZcandImpact_mu_EM;    
    
    float mllZcandImpact_el_EM;      
    float mTllZcandImpact_el_EM;
    float pTlZcandImpact_el_EM;
    float etalZcandImpact_el_EM;
    int IClZcandImpact_el_EM;
    float ptcone30lZcandImpact_el_EM;
    float etcone30lZcandImpact_el_EM;   
    float d0SiglZcandImpact_el_EM;
    float z0SinThetalZcandImpact_el_EM;
    
    bool ZcandLep_exists_mu_EM;
    bool ZcandLep_passesPT_mu_EM;
    bool ZcandLep_passesEta_mu_EM;
    bool ZcandLep_passesPTcone_mu_EM;
    bool ZcandLep_passesETcone_mu_EM;
    bool ZcandLep_passesD0_mu_EM; 
    bool ZcandLep_passesZ0_mu_EM; 
    bool ZcandLep_PassesMedium_mu_EM;
    bool ZcandLep_PassesTight_mu_EM; 
    bool ZcandLep_PassesORAndMllCut_mu_EM;
    bool ZcandLep_PassesPR_mu_EM;
    
    bool ZcandLep_exists_el_EM;
    bool ZcandLep_passesPT_el_EM;
    bool ZcandLep_passesEta_el_EM;
    bool ZcandLep_passesPTcone_el_EM;
    bool ZcandLep_passesETcone_el_EM;
    bool ZcandLep_passesD0_el_EM; 
    bool ZcandLep_passesZ0_el_EM; 
    bool ZcandLep_PassesMedium_el_EM;
    bool ZcandLep_PassesTight_el_EM; 
    bool ZcandLep_PassesORAndMllCut_el_EM;
    bool ZcandLep_PassesPR_el_EM;
    
    //#####################################
    
    ClassDef(TSelector_SusyNtuple_Truth, 1);

  protected:
    std::string         m_sel;          // event selection string

    DilTrigLogic*      m_trigObjWithoutRU;      // My trigger logic class
    SUSY::CrossSectionDB*                       m_susyXsec;     // SUSY cross section database

    // Cut variables    
    bool                m_selectB;      // switch to select b-tagged jets
    bool                m_vetoB;        // switch to veto b-tagged jets
    bool                m_selectSFOS;   // switch to select SFOS pairs
    bool                m_vetoSFOS;     // switch to veto SFOS pairs

    bool                m_writeOut;     // switch to control output dump

    
    /* RecoTruthMatch                m_recoTruthMatch;       // Lepton truth matching tool */
    
    unsigned int mcid_of_first_entry;
    float sumw_from_histo;
    
    bool m_kIsData;
    bool calcFakeContribution; 
    int isys;
    SusyNtSys SysSetting;
    BTagSys SysSettingBTag;
    int qflipSysUncertType;
    int SFUncertType;
    bool n0150BugFix;
    
    /* enum LEP_TYPE{PR=0, CONV, HF, LF, TYPE_Undef}; */
    
    // Build a vector that is the difference between two vectors.
// Caveat1: computing the difference triggers copies of the vectors
// Caveat2: the result is not guaranteed to be sorted
// Caveat3: a-b != b-a
// Based on: http://stackoverflow.com/questions/14175858/c-subtract-vectors
  template <typename T>
  std::vector<T> subtract_vector(std::vector<T>& a, const std::vector<T>& b)
  {
    std::vector<T> difference;
    std::set_difference(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter( difference ));
    return difference;
  }
  
   TTree *tree_out;
   
   size_t counter_input;
   size_t counter_event_cleaning;
   size_t counter_EE, counter_EM, counter_MM;
   size_t counter_EE_SRSS1, counter_EE_SRSS2;
   size_t counter_EM_SRSS1, counter_EM_SRSS2;
   size_t counter_MM_SRSS1, counter_MM_SRSS2;
 private:
    TFile *file_;
    TTree *tree_;
    FourMom l0_, l1_, met_;
    std::vector<FourMom> jets_, lowptLepts_;
    EventParameters eventPars_;
    
};

#endif
