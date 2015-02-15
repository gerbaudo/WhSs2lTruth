#include "Mt2/mt2_bisect.h"
#include "WhSs2lTruth/TSelector_SusyNtuple_Truth.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TChainElement.h"

#include <iomanip>
#include <string>
//#include <strstream>
using namespace std;
using namespace Susy;
// using FourMom;

/*--------------------------------------------------------------------------------*/
// TSelector_SusyNtuple_Truth Constructor
/*--------------------------------------------------------------------------------*/
TSelector_SusyNtuple_Truth::TSelector_SusyNtuple_Truth() :
        m_selectB (false),
        m_vetoB (false),
        m_selectSFOS(false),
        m_vetoSFOS (false),
        m_writeOut (false),
        tree_out(0),
        counter_input(0),
        counter_event_cleaning(0),
        counter_EE(0),
        counter_EM(0),
        counter_MM(0),
        counter_EE_SRSS1(0), counter_EE_SRSS2(0),
        counter_EM_SRSS1(0), counter_EM_SRSS2(0),
        counter_MM_SRSS1(0), counter_MM_SRSS2(0),
        file_(0),
        tree_(0)
{
//   setAnaType(Ana_2LepWH);
  setDebug(false);

  if(m_writeOut) {
    out.open("event.dump");
  }
}

/*--------------------------------------------------------------------------------*/
// The Begin() function is called at the start of the query.
/*--------------------------------------------------------------------------------*/
void TSelector_SusyNtuple_Truth::Begin(TTree * /*tree*/)
{

  n0150BugFix = true;
  TString option = GetOption();
  SysSetting = NtSys_NOM;


}
void TSelector_SusyNtuple_Truth::SlaveBegin(TTree* /*tree*/)
{
    calcSysUncert = true;
    makeNTuple = false;
    defineHistos();
    if(calcSysUncert) defineHistos_sysUncert();
    SusyNtTruthAna::Begin(0);
    m_trigObjWithoutRU = new DilTrigLogic("Moriond",false/*No Reweight Utils!*/);

    cout << "initialize chargeFlip tool" << endl;
    if(makeNTuple) initTupleMaker("/data/etp3/jwittkow/analysis_SUSYTools_03_04_SusyNt_01_16/SusySel_Egamma_6_NEW.root", "SusySel");
}



Bool_t TSelector_SusyNtuple_Truth::Process(Long64_t entry)
{
    GetEntry(entry); //function from SusyNtAna.h
    clearTruthObjects(); //function from SusyNtAna.cxx
    m_chainEntry++;
    counter_input++;
    unsigned int mcid;
//   cout << "nt.evt()->mcChannel= " << nt.evt()->mcChannel << endl;
    if(nt.evt()->isMC){
        mcid = nt.evt()->mcChannel;
        calcFakeContribution = false;
    } else {
        mcid = 111111;
        calcFakeContribution = true;
    }
    m_kIsData = !nt.evt()->isMC;


    if(m_dbg || m_chainEntry%50000==0) {
        // float recalc_sumw = 0.;
//     SumwMapKey genKey(mcid, 0);
//     map<unsigned int, float>::const_iterator sumwMapIter = m_sumwMap.find(genKey);
//     if(sumwMapIter != m_sumwMap.end()) recalc_sumw = sumwMapIter->second;
        cout << "**** Processing entry " << setw(6) << m_chainEntry
             << " run " << setw(6) << nt.evt()->run
             << " event " << setw(7) << nt.evt()->event
             << " sumw= " << nt.evt()->sumw << "   ****" << endl;
    }

    // charge flip background contribution in SS channels: for the
    // e^pm e^pm and e^pm mu^pm channels, processes that are
    // opposite-sign in truth but where one electron has undergone a
    // “charge flip”. Contributions from WW, ttbar, Z/gamma* and
    // single top are via charge-flip
    // In previous analysis, it has been observed that the charge flip
    // rate in data is lower than that in the simulation by about
    // 20%. Because of this disagreement, the electron charge flip
    // rate is measured in data as a function of |eta| and combined
    // with the smaller dependence on pT taken from simulation.


    //select Leptons, Jets, ... automatically with SusyNt methods:
    clearTruthObjects();
//   selectTruthObjects();
    const Susy::TruthMet* truthMet = getTruthMet(&nt);
    TLorentzVector truthMet_TLV = truthMet->lv();
    TruthParticleVector truthElectrons = getPreTruthLeptons(&nt,11); //10 GeV, |eta|<2.47
    TruthParticleVector truthMuons     = getPreTruthLeptons(&nt,13);// 10 GeV, |eta|<2.4
    TruthParticleVector truthTaus      = getPreTruthLeptons(&nt,15);// 20 GeV, |eta|<2.5
    TruthJetVector truthJets      = getPreTruthJets   (&nt   ); // 20 GeV, |eta|<4.5
    overlapProcedure(truthElectrons, truthMuons, truthTaus, truthJets);
    removeSFOSPair(truthElectrons);
    removeSFOSPair(truthMuons    );

    TruthParticleVector signal_truthElectrons = truthElectrons; //leptons: same as baseline selection.
    TruthParticleVector signal_truthMuons = truthMuons;
    TruthParticleVector signal_truthTaus = truthTaus;
    TruthJetVector signal_L20truthJets = getTruthJetsL20(truthJets); //jets: pt>20 GeV, |eta|<=2.4, not matched to b quark
    TruthJetVector signal_F30truthJets = getTruthJetsF30(truthJets); //jets: pt>30 GeV, |eta|>2.4
    TruthJetVector signal_B20truthJets = getTruthJetsB20(truthJets); //jets: pt>20 GeV, |eta|<=2.4, matched to b quark
    TruthJetVector signal_truthJets = getTruthJetsSignal(truthJets); //jets: pt>20 GeV, |eta|<=2.4 or pt > 30 GeV and |eta| > 2.4

    calcJet_variables(truthMet_TLV, signal_truthJets);

    int flag = nt.evt()->cutFlags[SysSetting];
    float weight_ALL = 1.;
    float weight_ALL_EE = weight_ALL;
    float weight_ALL_MM = weight_ALL;
    float weight_ALL_EM = weight_ALL;


    //Anders: I'm skipping many of the event cleaning cuts. Some of
    //them are simply not applicable to MC truth - There are no cosmic
    //muons, and there are no mismeasured or misidentified
    //objects. The goodPrimaryVertex cut, for example, is based on the
    //number of tracks leading to a given vertex. The corresponding
    //number does not exist in truth. All you have is the (x, y, z)
    //position of each vertex. It's possible that you could simulate
    //the vx_nTracks variable if you looked up the official definition
    //of how tracks are assigned to vertices, but not without
    //significant extra work.
    //event cleaning doesn't need to be changed after agreement on cutflow:
    if(doEventCleaning_andFillHistos(truthJets, truthMet, truthMuons, truthElectrons,
                                     flag, weight_ALL_EE, weight_ALL_MM, weight_ALL_EM,
                                     SysSetting, n0150BugFix)) {
        counter_event_cleaning++;
        //Set TLorentzVector for jets:
        TruthJet* jet0_buff=NULL;
        TruthJet* jet1_buff=NULL;
        TLorentzVector signalJet0_buff_TLV, signalJet1_buff_TLV;
        //nSignalJets = signal_truthJets.size(); // computed in calcJet_variables
        if(nSignalJets >=1){
            jet0_buff = signal_truthJets.at(0);
            signalJet0_buff_TLV.SetPtEtaPhiE(jet0_buff->Pt(), jet0_buff->eta, jet0_buff->phi,
                                             jet0_buff->Pt()*cosh(jet0_buff->eta));
            signalJet0_buff_TLV.SetPtEtaPhiM(jet0_buff->Pt(), jet0_buff->eta, jet0_buff->phi, jet0_buff->m);
            if(nSignalJets >=2){
                jet1_buff = signal_truthJets.at(1);
                signalJet1_buff_TLV.SetPtEtaPhiE(jet1_buff->Pt(), jet1_buff->eta, jet1_buff->phi,
                                                 jet1_buff->Pt()*cosh(jet1_buff->eta));
                signalJet1_buff_TLV.SetPtEtaPhiM(jet1_buff->Pt(), jet1_buff->eta, jet1_buff->phi,
                                                 jet1_buff->m);

            }
        }
        //now use buffer to fill jets variables (otherwise code crashes):
        TruthJet* jet0;
        TruthJet* jet1;
        TLorentzVector signalJet0_TLV, signalJet1_TLV;
        jet0 = jet0_buff;
        signalJet0_TLV = signalJet0_buff_TLV;
        if(nSignalJets >=2){
            jet0 = (jet0_buff->Pt() > jet1_buff->Pt()) ? jet0_buff : jet1_buff;
            jet1 = (jet0_buff->Pt() > jet1_buff->Pt()) ? jet1_buff : jet0_buff;
            signalJet1_TLV = (jet0_buff->Pt() > jet1_buff->Pt()) ? signalJet1_buff_TLV : signalJet0_buff_TLV;
            signalJet0_TLV = (jet0_buff->Pt() > jet1_buff->Pt()) ? signalJet0_buff_TLV : signalJet1_buff_TLV;
        }
        bool useForwardJets = true;
        //make a copy of electrons and MET which will eventually be changed by ChargeFlip tool:

        float cutnumber;

        //////////////// EE ////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        if(truthElectrons.size() == 2){
            cutnumber = 15.; fillHistos_EE_SRSS1(cutnumber, weight_ALL_EE);
            if(signal_truthElectrons.size() == 2){
                counter_EE++;
                cutnumber = 16.; fillHistos_EE_SRSS1(cutnumber, weight_ALL_EE);
                if(signal_truthTaus.size()==0){
                    cutnumber = 17.; fillHistos_EE_SRSS1(cutnumber, weight_ALL_EE);
                    if(( 14 < signal_truthElectrons.at(0)->Pt() &&
                         14 < signal_truthElectrons.at(1)->Pt() ) ||
                       (25 < signal_truthElectrons.at(0)->Pt() &&
                        10 < signal_truthElectrons.at(1)->Pt() &&
                        signal_truthElectrons.at(1)->Pt() <= 14 )) {
                        cutnumber = 18.; fillHistos_EE_SRSS1(cutnumber, weight_ALL_EE);
                        cutnumber = 19.; fillHistos_EE_SRSS1(cutnumber, weight_ALL_EE);
                        cutnumber = 20.; fillHistos_EE_SRSS1(cutnumber, weight_ALL_EE);
                        // genuine SS events
                        if(signal_truthElectrons.at(0)->charge * signal_truthElectrons.at(1)->charge > 0.){
                            cutnumber = 21.; fillHistos_EE_SRSS1(cutnumber, weight_ALL_EE);
                            // TruthParticleVector electrons_SS;
                            TruthParticle* el0_SS;
                            TruthParticle* el1_SS;
                            el0_SS = (signal_truthElectrons.at(0)->Pt() > signal_truthElectrons.at(1)->Pt()) ?
                                signal_truthElectrons.at(0) : signal_truthElectrons.at(1);
                            el1_SS = (signal_truthElectrons.at(0)->Pt() > signal_truthElectrons.at(1)->Pt()) ?
                                signal_truthElectrons.at(1) : signal_truthElectrons.at(0);
                            TLorentzVector el0_SS_TLV(*el0_SS);
                            TLorentzVector el1_SS_TLV(*el1_SS);


                            calc_EE_variables(el0_SS, el1_SS, el0_SS_TLV, el1_SS_TLV,
                                              signal_truthJets, signalJet0_TLV, signalJet1_TLV, truthMet_TLV);

                            if((mllZcandImpact_EE > MZ+20. || mllZcandImpact_EE < MZ-20.)){
                                cutnumber = 22.;  fillHistos_EE_SRSS1(cutnumber, weight_ALL_EE);
                                if(signal_F30truthJets.size()==0){
                                    cutnumber = 23.;  fillHistos_EE_SRSS1(cutnumber, weight_ALL_EE);
                                    if(signal_B20truthJets.size()==0){
                                        cutnumber = 24.;  fillHistos_EE_SRSS1(cutnumber, weight_ALL_EE);
                                        ///////////////////////////////////////////////////////////////////
                                        //                1jet                //////////////
                                        ///////////////////////////////////////////////////////////////////
                                        if(signal_truthJets.size()==1){
                                            cutnumber = 25.;  fillHistos_EE_SRSS1(cutnumber, weight_ALL_EE);
                                            if(el0_SS->Pt() > 30. && el1_SS->Pt() > 20.){
                                                cutnumber=26.; fillHistos_EE_SRSS1(cutnumber, weight_ALL_EE);

                                                if((el0_SS_TLV + el1_SS_TLV).M() > MZ+10. ||
                                                   (el0_SS_TLV + el1_SS_TLV).M() < MZ-10.){
                                                    cutnumber = 27.;  fillHistos_EE_SRSS1(cutnumber, weight_ALL_EE); //ZVeto
                                                    cutnumber = 28.;  fillHistos_EE_SRSS1(cutnumber, weight_ALL_EE);
                                                    cutnumber = 29.;  fillHistos_EE_SRSS1(cutnumber, weight_ALL_EE);
                                                    if(Mlj_EE < 90.){
                                                        cutnumber = 30.;  fillHistos_EE_SRSS1(cutnumber, weight_ALL_EE);
                                                        if(HT_EE > 200.){
                                                            cutnumber = 31.;  fillHistos_EE_SRSS1(cutnumber, weight_ALL_EE);
                                                            // from MET, calculate METrel

                                                            float METrel_EE = recalcMetRel(truthMet_TLV, el0_SS_TLV, el1_SS_TLV, signal_truthJets, useForwardJets);
                                                            if(METrel_EE > 55.){
                                                                cutnumber = 32.;  fillHistos_EE_SRSS1(cutnumber, weight_ALL_EE);
                                                                cutnumber = 33.;  fillHistos_EE_SRSS1(cutnumber, weight_ALL_EE);
                                                                counter_EE_SRSS1++;
                                                            } // metrel>55
                                                        }//ht>200
                                                    }//mlj>90
                                                }// mZ veto
                                            } //pt1>30 pt2>20
                                        } // nj==1
                                        ///////////////////////////////////////////////////////////////////
                                        //                2,3jet                //////////////
                                        ///////////////////////////////////////////////////////////////////
                                        if(signal_truthJets.size() >=2 && signal_truthJets.size() <=3){
                                            cutnumber = 34.;  fillHistos_EE_SRSS1(cutnumber, weight_ALL_EE);
                                            if(el0_SS->Pt() > 30. && el1_SS->Pt() > 20.){

                                                cutnumber = 35.;  fillHistos_EE_SRSS1(cutnumber, weight_ALL_EE);
                                                if((el0_SS_TLV + el1_SS_TLV).M() > MZ+10. ||
                                                   (el0_SS_TLV + el1_SS_TLV).M() < MZ-10.){
                                                    cutnumber = 36.;  fillHistos_EE_SRSS1(cutnumber, weight_ALL_EE); //ZVeto
                                                    cutnumber = 37.;  fillHistos_EE_SRSS1(cutnumber, weight_ALL_EE);
                                                    if(mTmax_EE > 110.){
                                                        cutnumber = 38.;  fillHistos_EE_SRSS1(cutnumber, weight_ALL_EE);
                                                        if(Mljj_EE < 120.){
                                                            cutnumber = 39.;  fillHistos_EE_SRSS1(cutnumber, weight_ALL_EE);
                                                            cutnumber = 40.;  fillHistos_EE_SRSS1(cutnumber, weight_ALL_EE);
                                                            float METrel_EE = recalcMetRel(truthMet_TLV,
                                                                                           el0_SS_TLV,
                                                                                           el1_SS_TLV,
                                                                                           signal_truthJets,
                                                                                           useForwardJets);
                                                            if(METrel_EE>=30.){
                                                                cutnumber = 41.;  fillHistos_EE_SRSS1(cutnumber, weight_ALL_EE);
                                                                counter_EE_SRSS2++;
                                                            } // metrel>30
                                                        } // mljj<120
                                                    } // mtmax >110
                                                } //m_Z veto
                                            } // pt1>30 pt2>20
                                        }// nj >=2
                                    } // veto b
                                } // veto f
                            } // veto 3rd lep
                        } // samesign
                    } // trigger pt
                } // veto tau
            } // n_sig_ele==2
        } // n_ele==2

        //////////////////// MM ////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        if(truthMuons.size() == 2){
            cutnumber = 15.; fillHistos_MM_SRSS1(cutnumber, weight_ALL_MM);
            if(signal_truthMuons.size()==2){
                counter_MM++;
                cutnumber = 16.; fillHistos_MM_SRSS1(cutnumber, weight_ALL_MM);
                if(signal_truthTaus.size()==0){
                    cutnumber = 17.; fillHistos_MM_SRSS1(cutnumber, weight_ALL_MM);

                    if(( 18 < signal_truthMuons.at(0)->Pt() &&
                         18 < signal_truthMuons.at(1)->Pt() ) ||
                       (18 < signal_truthMuons.at(0)->Pt() &&
                        14 < signal_truthMuons.at(1)->Pt() &&
                        signal_truthMuons.at(1)->Pt() <= 18) ||
                       (18 < signal_truthMuons.at(0)->Pt() &&
                        8 < signal_truthMuons.at(1)->Pt() &&
                        signal_truthMuons.at(1)->Pt() <= 14 ) ||
                       (14 < signal_truthMuons.at(0)->Pt() &&
                        signal_truthMuons.at(0)->Pt() <= 18 &&
                        14 < signal_truthMuons.at(1)->Pt() ) ){
                        cutnumber = 18.; fillHistos_MM_SRSS1(cutnumber, weight_ALL_MM);
                        cutnumber = 19.; fillHistos_MM_SRSS1(cutnumber, weight_ALL_MM);
                        cutnumber = 20.; fillHistos_MM_SRSS1(cutnumber, weight_ALL_MM);
                        if(signal_truthMuons.at(0)->charge * signal_truthMuons.at(1)->charge > 0.){
                            cutnumber = 21.; fillHistos_MM_SRSS1(cutnumber, weight_ALL_MM);
                            TruthParticle* mu0;
                            TruthParticle* mu1;

                            mu0 = (signal_truthMuons.at(0)->Pt() > signal_truthMuons.at(1)->Pt()) ?
                                signal_truthMuons.at(0) : signal_truthMuons.at(1);
                            mu1 = (signal_truthMuons.at(0)->Pt() > signal_truthMuons.at(1)->Pt()) ?
                                signal_truthMuons.at(1) : signal_truthMuons.at(0);
                            TLorentzVector mu0_TLV(*mu0);
                            TLorentzVector mu1_TLV(*mu1);

                            calc_MM_variables(mu0, mu1, mu0_TLV, mu1_TLV,
                                              signal_truthJets, signalJet0_TLV, signalJet1_TLV, truthMet_TLV);

                            if((mllZcandImpact_MM > MZ+20. || mllZcandImpact_MM < MZ-20.)){
                                cutnumber = 22.;  fillHistos_MM_SRSS1(cutnumber, weight_ALL_MM);
                                if(signal_F30truthJets.size()==0){
                                    cutnumber = 23.;  fillHistos_MM_SRSS1(cutnumber, weight_ALL_MM);
                                    if(signal_B20truthJets.size()==0){
                                        cutnumber = 24.;  fillHistos_MM_SRSS1(cutnumber, weight_ALL_MM);
                                        ///////////////////////////////////////////////////////////////
                                        ///////////       1 jet /////////////////////////////
                                        ///////////////////////////////////////////////////////////////
                                        if(signal_truthJets.size() ==1){
                                            cutnumber = 25.;  fillHistos_MM_SRSS1(cutnumber, weight_ALL_MM);
                                            if(mu0->Pt() > 30 && mu1->Pt() > 20){
                                                cutnumber = 26.;  fillHistos_MM_SRSS1(cutnumber, weight_ALL_MM);
                                                cutnumber = 27.;  fillHistos_MM_SRSS1(cutnumber, weight_ALL_MM); //ZVeto
                                                if(DeltaEtall_MM < 1.5){
                                                    cutnumber = 28.;  fillHistos_MM_SRSS1(cutnumber, weight_ALL_MM);
                                                    if(mTmax_MM > 110.){
                                                        cutnumber = 29.;  fillHistos_MM_SRSS1(cutnumber, weight_ALL_MM);
                                                        if(Mlj_MM < 90.){
                                                            cutnumber = 30.;  fillHistos_MM_SRSS1(cutnumber, weight_ALL_MM);

                                                            if(HT_MM > 200.){
                                                                cutnumber = 31.;  fillHistos_MM_SRSS1(cutnumber, weight_ALL_MM);
                                                                counter_MM_SRSS1++;
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }//end ==1 jets
                                        ///////////////////////////////////////////////////////////////
                                        ///////////       2,3 jets /////////////////////////////
                                        ///////////////////////////////////////////////////////////////
                                        if(signal_truthJets.size() >=2 && signal_truthJets.size() <=3){
                                            cutnumber = 34.;  fillHistos_MM_SRSS1(cutnumber, weight_ALL_MM);
                                            if(mu0->Pt() > 30 && mu1->Pt() > 30){
                                                cutnumber = 35.;  fillHistos_MM_SRSS1(cutnumber, weight_ALL_MM);
                                                cutnumber = 36.;  fillHistos_MM_SRSS1(cutnumber, weight_ALL_MM); //ZVeto
                                                if(DeltaEtall_MM<1.5){
                                                    cutnumber = 37.;  fillHistos_MM_SRSS1(cutnumber, weight_ALL_MM);
                                                    cutnumber = 38.;  fillHistos_MM_SRSS1(cutnumber, weight_ALL_MM);
                                                    if(Mljj_MM < 120.){
                                                        cutnumber = 39.;  fillHistos_MM_SRSS1(cutnumber, weight_ALL_MM);
                                                        if(HT_MM >= 200.){
                                                            cutnumber = 40.;  fillHistos_MM_SRSS1(cutnumber, weight_ALL_MM);
                                                            counter_MM_SRSS2++;
                                                        }
                                                    }
                                                }
                                            }
                                        }//end >=2 jets
                                    }
                                }
                            }
                        }
                    }
                }
            }            
        }

        //////////////////////////// EM ////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        if(truthElectrons.size() == 1 && truthMuons.size() == 1){
            cutnumber = 15.; fillHistos_EM_SRSS1(cutnumber, weight_ALL_EM);
            if(signal_truthElectrons.size()==1 && signal_truthMuons.size()==1){
                counter_EM++;
                cutnumber = 16.; fillHistos_EM_SRSS1(cutnumber, weight_ALL_EM);
                if(signal_truthTaus.size()==0){
                    cutnumber = 17.; fillHistos_EM_SRSS1(cutnumber, weight_ALL_EM);
                    if( ( 14 < signal_truthElectrons.at(0)->Pt() &&
                          8 < signal_truthMuons.at(0)->Pt()) ||
                        ( 10 < signal_truthElectrons.at(0)->Pt() &&
                          signal_truthElectrons.at(0)->Pt() <= 14 &&
                          18 < signal_truthMuons.at(0)->Pt()) ){
                        cutnumber = 18.; fillHistos_EM_SRSS1(cutnumber, weight_ALL_EM);
                        cutnumber = 19.; fillHistos_EM_SRSS1(cutnumber, weight_ALL_EM);
                        cutnumber = 20.; fillHistos_EM_SRSS1(cutnumber, weight_ALL_EM);
                        if(signal_truthElectrons.at(0)->charge * signal_truthMuons.at(0)->charge > 0.){
                            cutnumber = 21.; fillHistos_EM_SRSS1(cutnumber, weight_ALL_EM);
                            TruthParticle* mu;
                            TruthParticle* el;
                            mu = signal_truthMuons.at(0);
                            el = signal_truthElectrons.at(0);
                            TLorentzVector mu_TLV(*mu);
                            TLorentzVector el_TLV(*el);
                            calc_EM_variables(mu, el, mu_TLV, el_TLV,
                                              signal_truthJets, signalJet0_TLV, signalJet1_TLV, truthMet_TLV);

                            if((mllZcandImpact_mu_EM > MZ+20. || mllZcandImpact_mu_EM < MZ-20.) &&
                               (mllZcandImpact_el_EM > MZ+20. || mllZcandImpact_el_EM < MZ-20.)){
                                cutnumber = 22.; fillHistos_EM_SRSS1(cutnumber, weight_ALL_EM);
                                if(signal_F30truthJets.size()==0){
                                    cutnumber = 23.; fillHistos_EM_SRSS1(cutnumber, weight_ALL_EM);
                                    if(signal_B20truthJets.size()==0){
                                        cutnumber = 24.; fillHistos_EM_SRSS1(cutnumber, weight_ALL_EM);
                                        /////////////////////////////////////////////////////////////////
                                        //             1 jet ///////////////////////////////////
                                        /////////////////////////////////////////////////////////////////
                                        if(nSignalJets ==1){
                                            cutnumber = 25.; fillHistos_EM_SRSS1(cutnumber, weight_ALL_EM);
                                            if(el->Pt() > 30. && mu->Pt() > 30.){
                                                cutnumber = 26.; fillHistos_EM_SRSS1(cutnumber, weight_ALL_EM);
                                                cutnumber = 27.; fillHistos_EM_SRSS1(cutnumber, weight_ALL_EM);                                                if(DeltaEtall_EM < 1.5){
                                                    cutnumber = 28.; fillHistos_EM_SRSS1(cutnumber, weight_ALL_EM);
                                                    if(mTmax_EM > 110){
                                                        cutnumber = 29.; fillHistos_EM_SRSS1(cutnumber, weight_ALL_EM);
                                                        if(Mlj_EM < 90.){
                                                            cutnumber = 30.; fillHistos_EM_SRSS1(cutnumber, weight_ALL_EM);
                                                            if(HT_EM > 200.){
                                                                cutnumber = 31.; fillHistos_EM_SRSS1(cutnumber, weight_ALL_EM);
                                                                cutnumber = 32.; fillHistos_EM_SRSS1(cutnumber, weight_ALL_EM);
                                                                cutnumber = 33.; fillHistos_EM_SRSS1(cutnumber, weight_ALL_EM);
                                                                counter_EM_SRSS1++;
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }//end ==1jet
                                        /////////////////////////////////////////////////////////////////
                                        //             2,3 jets ///////////////////////////////////
                                        /////////////////////////////////////////////////////////////////
                                        if(nSignalJets >=2 && nSignalJets <=3){
                                            cutnumber = 34.; fillHistos_EM_SRSS1(cutnumber, weight_ALL_EM);
                                            if(el->Pt() > 30. && mu->Pt() > 30.){
                                                cutnumber = 35.; fillHistos_EM_SRSS1(cutnumber, weight_ALL_EM);
                                                cutnumber = 36.; fillHistos_EM_SRSS1(cutnumber, weight_ALL_EM); //ZVeto
                                                if(DeltaEtall_EM < 1.5){
                                                    cutnumber = 37.; fillHistos_EM_SRSS1(cutnumber, weight_ALL_EM);
                                                    if(mTmax_EM > 110){
                                                        cutnumber = 38.; fillHistos_EM_SRSS1(cutnumber, weight_ALL_EM);
                                                        if(Mljj_EM < 120.){
                                                            cutnumber = 39.; fillHistos_EM_SRSS1(cutnumber, weight_ALL_EM);
                                                            if(HT_EM > 200.){
                                                                cutnumber = 41.; fillHistos_EM_SRSS1(cutnumber, weight_ALL_EM);
                                                                counter_EM_SRSS2++;
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }//end >=2 jets
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return kTRUE;
}



/*--------------------------------------------------------------------------------*/
// Analysis cuts
/*--------------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------------*/

// bool TSelector_SusyNtuple_Truth::CheckRealLeptons(const ElectronVector& signal_electrons, MuonVector& signal_muons){

//   for(uint i=0; i<signal_electrons.size(); i++){
//     Electron* signal_electron = signal_electrons.at(i);
//     if(signal_electron->isChargeFlip) return false;
//     if(signal_electron->truthType != RecoTruthMatch::PROMPT) return false;
//   }

//     for(uint i=0; i<signal_muons.size(); i++){
//       Muon* signal_muon = signal_muons.at(i);
//       if(signal_muon->truthType != RecoTruthMatch::PROMPT) return false;
//   }

//   return true;

// }

/*--------------------------------------------------------------------------------*/
// bool TSelector_SusyNtuple_Truth::CheckChargeFlipElectrons(const ElectronVector& signal_electrons){

//   for(uint i=0; i<signal_electrons.size(); i++){
// Electron* signal_electron = signal_electrons.at(i);
// if(signal_electron->isChargeFlip) return false;
// // check if signal electron has no charge flip
//   }
//   return true;

// }
/*--------------------------------------------------------------------------------*/
float TSelector_SusyNtuple_Truth::getBTagWeight(const Event* evt, BTagSys SysSettingBTag)
{

//   JetVector tempJets = getBTagSFJets2Lep(m_baseJets);
//   return bTagSF(evt, tempJets, nt.evt()->mcChannel, SysSettingBTag);
return 0.;

}

//   if(!nt->evt()->isMC) return 1;
//   if(!USE_BWEIGHT)     return 1;
//
//   JetVector  valJets;
//   valJets.clear();
//   for(uint i=0; i<jets->size(); ++i){
//     Jet* jet = jets->at(i);
//     if( jet->Pt() < JET_PT_L20_CUT        ) continue;
//     if( fabs(jet->detEta) > JET_ETA_CUT_2L ) continue;
//     valJets.push_back(jet);
//   }
//
//   if(valJets.size()==0) return 1;//safety.
//
//   //Get sys naming convention
//   uint _sys = DGSys_NOM;
//   if(iSys==DGSys_BJet_DN) _sys=BTag_BJet_DN;
//   if(iSys==DGSys_CJet_DN) _sys=BTag_CJet_DN;
//   if(iSys==DGSys_LJet_DN) _sys=BTag_LJet_DN;
//   if(iSys==DGSys_BJet_UP) _sys=BTag_BJet_UP;
//   if(iSys==DGSys_CJet_UP) _sys=BTag_CJet_UP;
//   if(iSys==DGSys_LJet_UP) _sys=BTag_LJet_UP;
//
//   return bTagSF(nt->evt(),valJets, nt->evt()->mcChannel, (BTagSys) _sys);

/*--------------------------------------------------------------------------------*/
bool TSelector_SusyNtuple_Truth::checkLeptonPt(const LeptonVector& leptons)
{
    //cout << "leptons[0]->Pt()= " << leptons[0]->Pt() << " leptons[1]->Pt()= " << leptons[1]->Pt() << endl;
    if(leptons[0]->Pt()>leptons[1]->Pt()){
        if(leptons[0]->Pt()>35. && leptons[1]->Pt()>20.) return true;
        else return false;
    } else {
        if(leptons[1]->Pt()>35. && leptons[0]->Pt()>20.) return true;
        else return false;
    }
    return true;
}
/*--------------------------------------------------------------------------------*/
// Anyes: https://svnweb.cern.ch/trac/atlasinst/browser/Institutes/UCIrvine/ataffard/SusyWeakProdAna/trunk/Root/PhysicsTools.cxx#L197
float TSelector_SusyNtuple_Truth::mTWW(TLorentzVector _ll, TLorentzVector _nu, bool MvvTrue)
{
    float dphi = acos(cos(_ll.Phi()-_nu.Phi()));
    float mll = _ll.M();
    float mvv=0;
    if(!MvvTrue) mvv=mll;
    float mT=0;
    mT = sqrt( pow(mll,2) + pow(mvv,2)
               + 2*( sqrt(pow(_ll.Pt(),2) + pow(mll,2)) * sqrt(pow(_nu.Pt(),2) + pow(mvv,2))
                     - _ll.Pt()*_nu.Pt()*cos(dphi) ) );    
    return mT;
}
/*--------------------------------------------------------------------------------*/
//Anyes: https://svnweb.cern.ch/trac/atlasinst/browser/Institutes/UCIrvine/ataffard/SusyWeakProdAna/trunk/Root/ToyNt.cxx#L404
float TSelector_SusyNtuple_Truth::calcHT(TLorentzVector l1, TLorentzVector l2,
                                         TLorentzVector met, TruthJetVector &signalJets)
{
    float HT = 0;
    HT += l1.Pt();
    HT += l2.Pt();
    for(uint i=0; i<signalJets.size(); i++){
        if(signalJets[i]->Pt() > 20) HT += signalJets[i]->Pt();
    }
    HT += met.E();
    return HT;
// ht = Meff(*leptons, *jets, met, JET_PT_CUT); //The function is defined in SusyNtTools.cxx & is the sum pT of leptons, jets and met
// where JET_PT_CUT = 20 (the default in the function uses pt=40)

// float SusyNtTools::Meff(const LeptonVector& leps, const JetVector& jets, const Met* met)
// {
// float meff = 0;
// for(uint i=0; i<leps.size(); i++) meff += leps[i]->Pt();
// for(uint i=0; i<jets.size(); i++){
// if(jets[i]->Pt() > 40) meff += jets[i]->Pt();
// }
// meff += met->Et;
// return meff;
// }
}

/*--------------------------------------------------------------------------------*/
float TSelector_SusyNtuple_Truth::calcMeff(TLorentzVector met, const JetVector &signalJets){
  //Anyes: https://svnweb.cern.ch/trac/atlasinst/browser/Institutes/UCIrvine/ataffard/SusyWeakProdAna/trunk/Root/ToyNt.cxx#L404
  //mEff : same as HT but without the lepton
    float HT = 0;
// HT += l1.Pt();
// HT += l2.Pt();
    for(uint i=0; i<signalJets.size(); i++){
        if(signalJets[i]->Pt() > 40) HT += signalJets[i]->Pt();
    }
    HT += met.E();
    return HT;
}
/*--------------------------------------------------------------------------------*/
float TSelector_SusyNtuple_Truth::calcMt(TLorentzVector _l, TLorentzVector _nu)
{
    float dphi = acos(cos(_l.Phi()-_nu.Phi()));
    return sqrt(2*_l.Pt()*_nu.Pt() * (1- cos(dphi)) );
}
/*--------------------------------------------------------------------------------*/
float TSelector_SusyNtuple_Truth::recalcMetRel(TLorentzVector metLV, TLorentzVector l1, TLorentzVector l2,
                                               const TruthJetVector jets, bool useForward)
{
    //copied from SusyNtTools.cxx and modified to work with TLorenzVector
    float dPhi = TMath::Pi()/2.;
    if( fabs(metLV.DeltaPhi(l1)) < dPhi ) dPhi = fabs(metLV.DeltaPhi(l1));
    if( fabs(metLV.DeltaPhi(l2)) < dPhi ) dPhi = fabs(metLV.DeltaPhi(l2));
    for(uint ij=0; ij<jets.size(); ++ij){
        const TruthJet* jet = jets.at(ij);
//     if( !useForward && SusyNtAna::isForwardJet(jet) ) continue; // Use only central jets
        if( fabs(metLV.DeltaPhi( *jet )) < dPhi ) dPhi = fabs(metLV.DeltaPhi( *jet ));
    }
    return metLV.Et() * sin(dPhi);
}
/*--------------------------------------------------------------------------------*/
float TSelector_SusyNtuple_Truth::calcMT2(TLorentzVector metlv, TLorentzVector l0, TLorentzVector l1)
{
  //copied from SusyNtTools.cxx and modified to work with TLorenzVector
    double pTMiss[3] = {0.0, metlv.Px(), metlv.Py()};
    double pA[3] = {0.0, l0.Px(), l0.Py()};
    double pB[3] = {0.0, l1.Px(), l1.Py()};
    // Create Mt2 object
    mt2_bisect::mt2 mt2_event;
    mt2_event.set_momenta(pA,pB,pTMiss);
    mt2_event.set_mn(0); // LSP mass = 0 is Generic
    return mt2_event.get_mt2();
}
/*--------------------------------------------------------------------------------*/
float TSelector_SusyNtuple_Truth::calcMT2J(TLorentzVector metlv, TLorentzVector l0, TLorentzVector l1,
                                           TLorentzVector j0, TLorentzVector j1)
{
//   l0 -> l0 + jet_i
//   l1 -> l1 + jet_j
//   minimize for jet
    //copied from SusyNtTools.cxx and modified to work with TLorentzVector
    TLorentzVector alpha_a, alpha_b;
    //case 1:
    alpha_a = l0 + j0;
    alpha_b = l1 + j1;
    double pTMiss1[3] = {0.0, metlv.Px(), metlv.Py()};
    double pA1[3] = {alpha_a.M(), alpha_a.Px(), alpha_a.Py()};
    double pB1[3] = {alpha_b.M(), alpha_b.Px(), alpha_b.Py()};
    // Create Mt2 object
    mt2_bisect::mt2 mt2_event1;
    mt2_event1.set_momenta(pA1,pB1,pTMiss1);
    mt2_event1.set_mn(0); // LSP mass = 0 is Generic

    //case 2:
    alpha_a = l0 + j1;
    alpha_b = l1 + j0;
    double pTMiss2[3] = {0.0, metlv.Px(), metlv.Py()};
    double pA2[3] = {alpha_a.M(), alpha_a.Px(), alpha_a.Py()};
    double pB2[3] = {alpha_b.M(), alpha_b.Px(), alpha_b.Py()};
    // Create Mt2 object
    mt2_bisect::mt2 mt2_event2;
    mt2_event2.set_momenta(pA2,pB2,pTMiss2);
    mt2_event2.set_mn(0); // LSP mass = 0 is Generic
    double min_mt2 = min(mt2_event1.get_mt2(), mt2_event2.get_mt2());
    double return_value = (min_mt2 > 0.) ? min_mt2 : -1.;
    return return_value;   
}
/*--------------------------------------------------------------------------------*/
float TSelector_SusyNtuple_Truth::calcMZTauTau_coll(const TLorentzVector &signal_lep_0,
                                                    const TLorentzVector &signal_lep_1,
                                                    const TLorentzVector &met)
{   
    // IMPLEMENTATION OF COLLINEAR APPROXIMATION - validate code!
    float px0(signal_lep_0.Px()), py0(signal_lep_0.Py());
    float px1(signal_lep_1.Px()), py1(signal_lep_1.Py());
    float pxm(met.Px()), pym(met.Py());
    float num( px0*py1 - py0*px1 );
    float den1( py1*pxm - px1*pym + px0*py1 - py0*px1 );
    float den2( px0*pym - py0*pxm + px0*py1 - py0*px1 );
    float x1 = ( den1 != 0.0 ? (num/den1) : 0.0);
    float x2 = ( den2 != 0.0 ? (num/den2) : 0.0);
    // not guaranteed that this configuration is kinematically possible
    float returnvalue = x1*x2 > 0.0 ? (signal_lep_0+signal_lep_1).M() / std::sqrt(x1*x2) : -1.0;
    return returnvalue;
}
/*--------------------------------------------------------------------------------*/
float TSelector_SusyNtuple_Truth::calcMZTauTau_mmc(TLorentzVector lep1, TLorentzVector lep2,
                                                   int tau0_decay_type, int tau1_decay_type)
{
    return 111.;
}

/*--------------------------------------------------------------------------------*/
float TSelector_SusyNtuple_Truth::calcSumMv1(const JetVector &signalJets)
{
    //with or without leptons?
    float sumMv1 = 0.;
    for(uint j=0; j<signalJets.size(); ++j){
        Jet* jet = signalJets.at(j);
        sumMv1 += jet->mv1;
    }
    return sumMv1;
}
/*--------------------------------------------------------------------------------*/
float TSelector_SusyNtuple_Truth::calc_D0(bool unbiased, const Lepton* lep)
{
  float d0_branch;
  if(unbiased){
    d0_branch = lep->d0Unbiased;
  } else {
      d0_branch = lep->d0;
  }
  return d0_branch;
}
/*--------------------------------------------------------------------------------*/
bool TSelector_SusyNtuple_Truth::compareElecMomentum (Electron* e0, Electron* e1)
{
    return (e0->Pt() > e1->Pt());
}
/*--------------------------------------------------------------------------------*/
bool TSelector_SusyNtuple_Truth::compareMuonMomentum (Muon* mu0, Muon* mu1)
{
    return (mu0->Pt() > mu1->Pt());
}
/*--------------------------------------------------------------------------------*/
bool TSelector_SusyNtuple_Truth::isCMSJet(const Susy::Jet* jet)
{
    if(jet->Pt() < 30.) return false;
    //if(fabs(jet->Eta()) > JET_ETA_CUT_2L) return false;
    if(fabs(jet->detEta) > JET_ETA_CUT_2L) return false;
    if(jet->Pt() < JET_JVF_PT &&
       fabs(jet->jvf) - 1e-3 < JET_JVF_CUT_2L) return false;
    if(jet->mv1 > MV1_80) return false;
    return true;
}
/*--------------------------------------------------------------------------------*/
int TSelector_SusyNtuple_Truth::numberOfCMSJets(const JetVector& jets)
{
    int ans = 0;
    for(uint ij=0; ij<jets.size(); ++ij){
        const Jet* j = jets.at(ij);
        if(isCMSJet(j)){
            ans++;
        }
    }
    return ans;
}
/*--------------------------------------------------------------------------------*/
vector<TLorentzVector> TSelector_SusyNtuple_Truth::overlapRemoval(vector<TLorentzVector> objects_type1,
                                                                  vector<TLorentzVector> objects_type2,
                                                                  double dr, bool sameType, bool removeSoft)
{
    //copied from MultiLep/CutflowTools.h and modified to work simply with vector<TLorentzVector>
    vector<TLorentzVector> survivors;
    for(unsigned int i=0; i<objects_type1.size(); i++) {
        bool is_overlap = false;
        for(unsigned int j=0; j<objects_type2.size(); j++) {
            // Don't remove an object against itself (for electron-electron overlap)
            if(sameType && i==j) continue;
            if (removeSoft) {
            } else {
                if(objects_type1.at(i).DeltaR(objects_type2.at(j)) <= dr) {
                    is_overlap = true;
                }
            }
        }
        if(is_overlap) {
            continue;
        }
        survivors.push_back(objects_type1.at(i));
    }
    return survivors;
}
/*--------------------------------------------------------------------------------*/
bool TSelector_SusyNtuple_Truth::doEventCleaning_andFillHistos(TruthJetVector baseJets,
                                                               const TruthMet* met,
                                                               TruthParticleVector baseMuons,
                                                               TruthParticleVector baseElectrons,
                                                               int flag,
                                                               float weight_ALL_EE,
                                                               float weight_ALL_MM,
                                                               float weight_ALL_EM,
                                                               SusyNtSys SysSetting,
                                                               bool n0150BugFix)
{
    float cutnumber = 0.;
    fillHistos_EE(cutnumber, weight_ALL_EE);
    fillHistos_MM(cutnumber, weight_ALL_MM);
    fillHistos_EM(cutnumber, weight_ALL_EM); // all events in the sample  
    if( !(flag & ECut_GRL) ) return false;
    cutnumber = 1.;
    fillHistos_EE(cutnumber, weight_ALL_EE);
    fillHistos_MM(cutnumber, weight_ALL_MM);
    fillHistos_EM(cutnumber, weight_ALL_EM); //grl Cut

    if( !(ECut_TileTrip & flag) ) return false;
    cutnumber = 2.;
    fillHistos_EE(cutnumber, weight_ALL_EE);
    fillHistos_MM(cutnumber, weight_ALL_MM);
    fillHistos_EM(cutnumber, weight_ALL_EM); //TileTripReader

    // if ( !(flag & ECut_TTC)) return false;
    cutnumber = 3.;
    fillHistos_EE(cutnumber, weight_ALL_EE);
    fillHistos_MM(cutnumber, weight_ALL_MM);
    fillHistos_EM(cutnumber, weight_ALL_EM);//IncompleteEvents Veto

    if(!(ECut_LarErr & flag) || !(ECut_TileErr & flag)) return false;
    cutnumber = 4.;
    fillHistos_EE(cutnumber, weight_ALL_EE);
    fillHistos_MM(cutnumber, weight_ALL_MM);
    fillHistos_EM(cutnumber, weight_ALL_EM); //LAr/TileError

    // if( !(flag & ECut_HotSpot)) return false; //remove event where a jet points into hot TileCal module
    cutnumber = 5.;
    fillHistos_EE(cutnumber, weight_ALL_EE);
    fillHistos_MM(cutnumber, weight_ALL_MM);
    fillHistos_EM(cutnumber, weight_ALL_EM); //TileCalHotSpot

    // remove event with VeryLooseBad jets with pT>20 GeV. No eta
    // cut. Only consider jets which are not overlapping with
    // electrons or taus:
    //if(SusyNtAna::hasBadJet(baseJets))return kFALSE;
    cutnumber = 6.;
    fillHistos_EE(cutnumber, weight_ALL_EE);
    fillHistos_MM(cutnumber, weight_ALL_MM);
    fillHistos_EM(cutnumber, weight_ALL_EM);//BadJets

    // on top of smart veto, veto event with >=1 jets before
    // electron-jet overlap removal with pT>40 GeV, BCH_CORR_JET >
    // 0.05, DeltaPhi(met,jet)<0.3 (Anyes)
    // TruthJetVector prejets = getPreTruthJets(&nt);
    // if(!passDeadRegions(prejets, met, nt.evt()->run, nt.evt()->isMC)/* || !(flag & ECut_SmartVeto)*/) return false; // SusyNtTools: passDeadRegions(...)
    cutnumber = 7.;
    fillHistos_EE(cutnumber, weight_ALL_EE);
    fillHistos_MM(cutnumber, weight_ALL_MM);
    fillHistos_EM(cutnumber, weight_ALL_EM); //CaloJets

    // if( !(flag & ECut_GoodVtx)) return false;
    cutnumber = 8.;
    fillHistos_EE(cutnumber, weight_ALL_EE);
    fillHistos_MM(cutnumber, weight_ALL_MM);
    fillHistos_EM(cutnumber, weight_ALL_EM);//PrimaryVertex

    // MuonVector preMuons = SusyNtAna::getPreMuons(&nt, SysSetting, n0150BugFix);
    // if( SusyNtAna::hasBadMuon(preMuons)) return false;
    cutnumber = 9.;
    fillHistos_EE(cutnumber, weight_ALL_EE);
    fillHistos_MM(cutnumber, weight_ALL_MM);
    fillHistos_EM(cutnumber, weight_ALL_EM);//BadMuons

    // if(hasCosmicMuon(baseMuons)) return false; // !(flag & ECut_Cosmic) - no longer guarantee the event flags that are stored :-(
    // preMuons.clear();
    cutnumber = 10.;
    fillHistos_EE(cutnumber, weight_ALL_EE);
    fillHistos_MM(cutnumber, weight_ALL_MM);
    fillHistos_EM(cutnumber, weight_ALL_EM);//Cosmic Muons


    // Sherpa WW fix, remove radiative b-quark processes that overlap
    // with single top: already done upstream in SusyCommon
    // SusyNtMaker

    if(nt.evt()->hfor == 4) return false; //remove events where same heavy flavor final states arise in multiple samples when combining ALPGEN samples
    cutnumber = 11.;
    fillHistos_EE(cutnumber, weight_ALL_EE);
    fillHistos_MM(cutnumber, weight_ALL_MM);
    fillHistos_EM(cutnumber, weight_ALL_EM);//hfor veto

    if(!(baseMuons.size() >= 2 ||
         baseElectrons.size() >= 2 ||
         (baseElectrons.size()+baseMuons.size()) >= 2)) return false;
    cutnumber = 12.;
    fillHistos_EE(cutnumber, weight_ALL_EE);
    fillHistos_MM(cutnumber, weight_ALL_MM);
    fillHistos_EM(cutnumber, weight_ALL_EM);//at least 2 base leptons

    if(!(baseMuons.size() == 2 ||
         baseElectrons.size() == 2 ||
         (baseElectrons.size()+baseMuons.size()) == 2)) return false;
    cutnumber = 13.;
    fillHistos_EE(cutnumber, weight_ALL_EE);
    fillHistos_MM(cutnumber, weight_ALL_MM);
    fillHistos_EM(cutnumber, weight_ALL_EM);//exactly 2 base leptons

    // if(!(m_baseElectrons.size()==2 || baseMuons.size()==2 || (m_baseElectrons.size()+baseMuons.size())==2)) return false; //only count leptons where no Mll < 20 GeV. Mll < 12 GeV veto ALREADY DONE FOR SELECTING BASELINE LEPTONS in performOverlapRemoval()


    // discard any SFOS baseline lepton pair with M_ll < 12 GeV unnecessary: already done when skimming/slimming ntuples
    // Any lepton pairs are required to have an invariant mass, m``, above 20 GeV such that to remove low-mass dilepton resonances
    if(baseMuons.size() == 2 && SusyNtAna::Mll(baseMuons[0], baseMuons[1]) < 20 )return false;
    if(baseElectrons.size() == 2 && SusyNtAna::Mll(baseElectrons[0], baseElectrons[1]) < 20 )return false;
    if((baseElectrons.size() ==1 && baseMuons.size() == 1) &&
       SusyNtAna::Mll(baseElectrons[0], baseMuons[0]) < 20 )return false;
    cutnumber = 14.;
    fillHistos_EE(cutnumber, weight_ALL_EE);
    fillHistos_MM(cutnumber, weight_ALL_MM);
    fillHistos_EM(cutnumber, weight_ALL_EM);//Mll
    return true;
}
//----------------------------------------------------------
bool TSelector_SusyNtuple_Truth::initTupleMaker(const std::string &outFilename, const std::string &treename)
{
    /*
      if(file_ && file_->IsOpen() && tree_) {
      cout<<"TupleMaker::init: already initialized"<<endl;
      return false;
      }
      else
    */
    return (initFile(outFilename) && initTree(treename));
}
//----------------------------------------------------------
bool TSelector_SusyNtuple_Truth::initFile(const std::string &outFilename)
{
    file_ = TFile::Open(outFilename.c_str(), "recreate");
    if(!file_) cout<<"TupleMaker::initFile('"<<outFilename<<"') : failed to create file"<<endl;
    cout << "init " << outFilename << endl;
    return (file_ && file_->IsOpen());
}
//----------------------------------------------------------
bool TSelector_SusyNtuple_Truth::initTree(const std::string &treename)
{
    bool initialized(false);
    TDirectory *startingDir = gDirectory;
    if(file_) {
        file_->cd();
        string title("TupleMaker tree");
        tree_ = new TTree(treename.c_str(), title.c_str());
        tree_->SetDirectory(file_);
        initTreeBranches();
        initialized = true;
	cout << "init " << treename << endl;
    } else {
        cout<<"TupleMaker::initTree: invalid file, failed to create tree"<<endl;
    }
    if(startingDir) startingDir->cd(); // root is easily confused by pwd; cd back to where we were
    return initialized;
}
//----------------------------------------------------------
bool TSelector_SusyNtuple_Truth::initTreeBranches()
{
    bool initialized(false);
    if(tree_) {
        tree_->Branch("l0",    &l0_); //l0_, l1_, met_: FourMom() : px(0), py(0), pz(0), E(0), isMu(false), isEl(false), isJet(false) {}
        tree_->Branch("l1",    &l1_);
        tree_->Branch("met",   &met_);
        tree_->Branch("jets",  &jets_); //vector<FourMom> jets_, lowptLepts_
        tree_->Branch("lepts", &lowptLepts_);
        tree_->Branch("pars",  &eventPars_); //EventParameters() : weight(0), eventNumber(0), runNumber(0) {}
    } else {
        cout<<"TupleMaker::initTreeBranches : invalid tree, failed to init branches"<<endl;
    }
    return initialized;
}
//----------------------------------------------------------
FourMom lepton2FourMom (const Lepton *l)
{
    return (  l && l->isMu()  ? FourMom().setMu(*l)
              : l && l->isEle() ? FourMom().setEl(*l)
              : FourMom());
}
//----------------------------------------------------------
FourMom jet2FourMom (const Jet *j) { return (j ? FourMom().setJet(*j) : FourMom()); }
//----------------------------------------------------------
bool TSelector_SusyNtuple_Truth::fillTupleMaker(const double weight,
                                                const unsigned int run,
                                                const unsigned int event,
                                                const bool isMc,
                                                const Susy::Lepton &l0,
                                                const Susy::Lepton &l1,
                                                const Susy::Met &met,
                                                const LeptonVector &otherLeptons,
                                                const JetVector &jets)
{
    bool someBytesWritten(false);
    if(tree_) {
//       cout << "fillTupleMaker" << endl;
        eventPars_.setWeight(weight).setRun(run).setEvent(event).setIsmc(isMc);
        l0.isMu() ? l0_.setMu(l0) : l0_.setEl(l0);
        l1.isMu() ? l1_.setMu(l1) : l1_.setEl(l1);
        met_.setMet(met);
        jets_.clear();
        lowptLepts_.clear();
        const LeptonVector &olps = otherLeptons;
        std::transform(jets.begin(), jets.end(), std::back_inserter(jets_),       jet2FourMom);
        std::transform(olps.begin(), olps.end(), std::back_inserter(lowptLepts_), lepton2FourMom);
        someBytesWritten = (tree_->Fill()>0);
    }
    return someBytesWritten;
}
/*--------------------------------------------------------------------------------*/
LeptonVector TSelector_SusyNtuple_Truth::getAnyElOrMu(SusyNtObject &susyNt,
                                                      const Lepton *l0,
                                                      const Lepton *l1)
{
    // DG 2013-12-02:
    // todo1 : re-implement with std algo
    // todo2 : re-implement with syst
    LeptonVector leptons;
    for(uint ie=0; ie<susyNt.ele()->size(); ++ie){
        if(Electron* e = & susyNt.ele()->at(ie)){ //e->setState(sys);
            if(fabs(e->d0Sig(true)) >= ELECTRON_D0SIG_CUT_WH) continue;
            if(fabs(e->z0SinTheta(true)) >= ELECTRON_Z0_SINTHETA_CUT) continue;
            if((e->q * l0->q)<0. || (e->q * l1->q)<0.){
                leptons.push_back(static_cast<Lepton*>(e));
            }
        }
    }
    for(uint im=0; im<susyNt.muo()->size(); ++im){
        if(Muon* m = & susyNt.muo()->at(im)){ //m->setState(sys);
            if(fabs(m->d0Sig(true)) >= MUON_D0SIG_CUT) continue;
            if(fabs(m->z0SinTheta(true)) >= MUON_Z0_SINTHETA_CUT) continue;
            if((m->q * l0->q)<0. || (m->q * l1->q)<0.){
                leptons.push_back(static_cast<Lepton*>(m));
            }
        }
    }
    return leptons;
}
/*--------------------------------------------------------------------------------*/
bool TSelector_SusyNtuple_Truth::closeTupleMaker()
{
    bool closed(false);
    if(file_) {
        file_->cd();
        file_->Write();
        file_->Close();
        file_->Delete();
        file_ = 0;
        closed = true;
        cout << "closeTupleMaker" << endl;
    } else {
        cout<<"TupleMaker::close : file not there"<<endl;
    }
    return closed;
}


/*--------------------------------------------------------------------------------*/
char *TSelector_SusyNtuple_Truth::GetID(Int_t type)
{
    return Form ("%d", type);
}

/*--------------------------------------------------------------------------------*/
// Debug event
/*--------------------------------------------------------------------------------*/
bool TSelector_SusyNtuple_Truth::debugEvent()
{
    uint run = nt.evt()->run;
    uint evt = nt.evt()->event;
    //if(run==191139 && evt==140644832) return true;
    if(run==180164&&evt==24769) return true;
    return false;
}
/*--------------------------------------------------------------------------------*/
float TSelector_SusyNtuple_Truth::getLeptonSF(const LeptonVector& leptons, int SFUncertType)
{
    float sf = 1.;
    for(uint i=0; i<leptons.size(); i++){
        const Lepton* lep = leptons[i];
        if(lep->isEle() && SFUncertType == 1)      sf *= (lep->effSF + lep->errEffSF);
        else if(lep->isEle() && SFUncertType == 2) sf *= (lep->effSF - lep->errEffSF);        
        else if(lep->isMu() && SFUncertType == 3)  sf *= (lep->effSF + lep->errEffSF);
        else if(lep->isMu() && SFUncertType == 4)  sf *= (lep->effSF - lep->errEffSF);
        else sf *= lep->effSF;
    }
    return sf;
}

/*--------------------------------------------------------------------------------*/
void TSelector_SusyNtuple_Truth::overlapProcedure(TruthParticleVector& elecs, TruthParticleVector& muons, TruthParticleVector& taus, TruthJetVector& jets)
{
    // Remove electrons from electrons
    e_e_overlap(elecs, E_E_DR);
    // Remove jets from electrons
    e_j_overlap(elecs, jets, J_E_DR, true);
    // Remove taus from electrons
    t_e_overlap(taus, elecs, T_E_DR);
    // Remove taus from muons
    t_m_overlap(taus, muons, T_M_DR);
    // Remove electrons from jets
    e_j_overlap(elecs, jets, E_J_DR, false);
    // Remove muons from jets
    m_j_overlap(muons, jets, M_J_DR);
    // Remove electrons and muons that overlap
    e_m_overlap(elecs, muons, E_M_DR);
    // Remove muons from muons
    m_m_overlap(muons, M_M_DR);
    // Remove jets from taus
    t_j_overlap(taus, jets, J_T_DR, true);
}
/*--------------------------------------------------------------------------------*/
void TSelector_SusyNtuple_Truth::e_e_overlap(TruthParticleVector& elecs, float minDr)
{
    uint nEl = elecs.size();
    if(nEl < 2) return;
    // Find all possible pairings
    static std::set<const TruthParticle*> removeElecs;
    removeElecs.clear();
    for(uint iEl=0; iEl<nEl; iEl++){
        const TruthParticle* e1 = elecs[iEl];
        for(uint jEl=iEl+1; jEl<nEl; jEl++){
            const TruthParticle* e2 = elecs[jEl];
            if(e1->DeltaR(*e2) < minDr){
                if(e1->Pt() < e2->Pt()){
                    removeElecs.insert(e1);
                    break;
                }else{
                    removeElecs.insert(e2);
                }
            } // dR
        } // e2 loop
    } // e1 loop

    // Remove electrons that overlap
    for(int iEl=nEl-1; iEl>=0; iEl--){
        if(removeElecs.find(elecs[iEl]) != removeElecs.end()){
            elecs.erase( elecs.begin() + iEl );
        }
    }
}
/*--------------------------------------------------------------------------------*/
void TSelector_SusyNtuple_Truth::e_j_overlap(TruthParticleVector& elecs, TruthJetVector& jets, float minDr, bool removeJets)
{
//   cout << nt.evt()->event << "  removeJets= " << removeJets << " minDr= " << minDr << endl;
  if(elecs.size()==0 || jets.size()==0) return;

  for(int ie=elecs.size()-1; ie>=0; ie--){
    const TruthParticle* e = elecs.at(ie);
    for(int ij=jets.size()-1; ij>=0; ij--){
        const TruthJet* j = jets.at(ij);
//       cout << "e->DeltaR(*j)= " << e->DeltaR(*j) << " e->Pt()= " << e->Pt() << " j->Pt()= " << j->Pt() << endl;
        if(e->DeltaR(*j) > minDr) continue;       
        if(removeJets){
            jets.erase( jets.begin() + ij );
        }else{
            elecs.erase( elecs.begin() + ie );
            break;
        }
    }// end loop over jets
  }// end loop over electrons
}
/*--------------------------------------------------------------------------------*/
void TSelector_SusyNtuple_Truth::m_j_overlap(TruthParticleVector& muons, TruthJetVector jets, float minDr)
{
  if(muons.size()==0 || jets.size()==0) return;
  for(int im=muons.size()-1; im>=0; im--){
      const TruthParticle* mu = muons.at(im);
      for(int ij=jets.size()-1; ij>=0; ij--){
          const TruthJet* j = jets.at(ij);
          if(mu->DeltaR(*j) > minDr) continue;
          muons.erase( muons.begin() + im );
          break;
      }// end loop over jets
  }// end loop over muons
}
/*--------------------------------------------------------------------------------*/
void TSelector_SusyNtuple_Truth::e_m_overlap(TruthParticleVector& elecs,
                                             TruthParticleVector& muons,
                                             float minDr)
{
    uint nEl = elecs.size();
    uint nMu = muons.size();
    if(nEl==0 || nMu==0) return;
    // Electron muon overlap should be pretty rare,
    // so we can take advantage of that and optimize
    static std::set<const TruthParticle*> removeElecs;
    static std::set<const TruthParticle*> removeMuons;
    removeElecs.clear();
    removeMuons.clear();
    // In this case we will want to remove both the electron and the muon
    for(uint iEl=0; iEl<nEl; iEl++){
        const TruthParticle* e = elecs[iEl];
        for(uint iMu=0; iMu<nMu; iMu++){
            const TruthParticle* mu = muons[iMu];
            if(e->DeltaR(*mu) < minDr){
                removeElecs.insert(e);
                removeMuons.insert(mu);
            }
        }
    }
    // Remove those electrons flagged for removal
    if(removeElecs.size()){
        for(int iEl=nEl-1; iEl>=0; iEl--){
            if(removeElecs.find(elecs[iEl])!=removeElecs.end()){
                elecs.erase( elecs.begin() + iEl );
            }
        }
    }
    // Remove those muons flagged for removal
    if(removeMuons.size()){
        for(int iMu=nMu-1; iMu>=0; iMu--){
            if(removeMuons.find(muons[iMu])!=removeMuons.end()){
                muons.erase( muons.begin() + iMu );
            }
        }
    }
}
/*--------------------------------------------------------------------------------*/
void TSelector_SusyNtuple_Truth::m_m_overlap(TruthParticleVector& muons, float minDr)
{
    uint nMu = muons.size();
    if(nMu < 2) return;
    // If 2 muons overlap, toss them both!
    static std::set<const TruthParticle*> removeMuons;
    removeMuons.clear();
    for(uint iMu=0; iMu<nMu; iMu++){
        const TruthParticle* mu1 = muons[iMu];
        for(uint jMu=iMu+1; jMu<nMu; jMu++){
            const TruthParticle* mu2 = muons[jMu];
            if(mu1->DeltaR(*mu2) < minDr){
                removeMuons.insert(mu1);
                removeMuons.insert(mu2);
            }
        }
    }
    for(int iMu=nMu-1; iMu>=0; iMu--){
        const TruthParticle* mu = muons[iMu];
        if(removeMuons.find(mu) != removeMuons.end()){
            muons.erase( muons.begin() + iMu );
        }
    }
}
/*--------------------------------------------------------------------------------*/
void TSelector_SusyNtuple_Truth::t_e_overlap(TruthParticleVector& taus,
                                             TruthParticleVector& elecs,
                                             float minDr)
{
    uint nTau = taus.size();
    uint nEle = elecs.size();
    if(nTau==0 || nEle==0) return;    
    for(int iTau=nTau-1; iTau>=0; iTau--){
        const TruthParticle* tau = taus[iTau];
        for(int iEl=nEle-1; iEl>=0; iEl--){
            const TruthParticle* e = elecs[iEl];
            if(tau->DeltaR(*e) < minDr){
                taus.erase( taus.begin() + iTau );
                break;
            }            
        }
    }
}
/*--------------------------------------------------------------------------------*/
void TSelector_SusyNtuple_Truth::t_m_overlap(TruthParticleVector& taus,
                                             TruthParticleVector& muons,
                                             float minDr)
{
    uint nTau = taus.size();
    uint nMuo = muons.size();
    if(nTau==0 || nMuo==0) return;
    for(int iTau=nTau-1; iTau>=0; iTau--){
        const TruthParticle* tau = taus[iTau];
        for(int iMu=nMuo-1; iMu>=0; iMu--){
            const TruthParticle* mu = muons[iMu];            
            if(tau->DeltaR(*mu) < minDr){
                taus.erase( taus.begin() + iTau );
                break;
            }            
        }
    }
}
/*--------------------------------------------------------------------------------*/
void TSelector_SusyNtuple_Truth::t_j_overlap(TruthParticleVector& taus,
                                             TruthJetVector& jets,
                                             float minDr,
                                             bool removeJets)
{
    uint nTau = taus.size();
    uint nJet = jets.size();
    if(nTau==0 || nJet==0) return;    
    for(int iTau=nTau-1; iTau>=0; iTau--){
        const TruthParticle* tau = taus.at(iTau);
        for(int iJet=jets.size()-1; iJet>=0; iJet--){
            const TruthJet* jet = jets.at(iJet);            
            if(tau->DeltaR(*jet) < minDr){
                if(removeJets)
                    jets.erase( jets.begin() + iJet );
                else{
                    taus.erase( taus.begin() + iTau );
                    break;
                }
            }
        } // end loop over jets
    } // end loop over electrons
}
/*--------------------------------------------------------------------------------*/
TruthJetVector TSelector_SusyNtuple_Truth::getTruthJetsL20   (TruthJetVector& truthBaseJets)
{
    TruthJetVector truthSignalJets;
    for(uint index=0; index<truthBaseJets.size(); ++index){
        TruthJet* jet = truthBaseJets.at(index);
        bool bmatch = false;
        TruthJet* truthjet = truthBaseJets.at(index);
        for(uint ij=0; ij<nt.jet()->size(); ++ij){
            Jet* recoj = & nt.jet()->at(ij);
            if(recoj->truthLabel==5){
                if(recoj->DeltaR(*truthjet) < 0.4) bmatch = true;
            }
        }
        if(!bmatch &&  jet->Pt() >= 20. && fabs(jet->Eta()) <= 2.4 ) truthSignalJets.push_back(jet);
    }
    return truthSignalJets;
}
/*--------------------------------------------------------------------------------*/
TruthJetVector TSelector_SusyNtuple_Truth::getTruthJetsF30   (TruthJetVector& truthBaseJets)
{
    TruthJetVector truthSignalJets;
    for(uint index=0; index<truthBaseJets.size(); ++index){
        TruthJet* jet = truthBaseJets.at(index);
        if(jet->Pt() >= 30. && fabs(jet->Eta()) > 2.4 ) truthSignalJets.push_back(jet);
    }
    return truthSignalJets;
}
/*--------------------------------------------------------------------------------*/
TruthJetVector TSelector_SusyNtuple_Truth::getTruthJetsB20   (TruthJetVector& truthBaseJets)
{
    TruthJetVector truthSignalJets;
    for(uint index=0; index<truthBaseJets.size(); ++index){
        bool bmatch = false;
        TruthJet* truthjet = truthBaseJets.at(index);
        for(uint ij=0; ij<nt.jet()->size(); ++ij){
            Jet* recoj = & nt.jet()->at(ij);
            if(recoj->truthLabel==5){
                if(recoj->DeltaR(*truthjet) < 0.4) bmatch = true;
            }
        }
        if( bmatch && truthjet->Pt() >= 20. &&
            fabs(truthjet->Eta()) <= 2.4)
            truthSignalJets.push_back(truthjet);
    }
    return truthSignalJets;
}
/*--------------------------------------------------------------------------------*/
TruthJetVector TSelector_SusyNtuple_Truth::getTruthJetsSignal   (TruthJetVector& truthBaseJets)
{
    TruthJetVector truthSignalJets;
    for(uint index=0; index<truthBaseJets.size(); ++index){
        TruthJet* jet = truthBaseJets.at(index);
        if(jet->Pt() >= 20. && fabs(jet->Eta()) <= 2.4 ) truthSignalJets.push_back(jet);
        if(jet->Pt() >= 30. && fabs(jet->Eta()) > 2.4  ) truthSignalJets.push_back(jet);
    }
    return truthSignalJets;
}
/*--------------------------------------------------------------------------------*/
void TSelector_SusyNtuple_Truth::SlaveTerminate()
{
    // SusyNtAna::Terminate();
    if(makeNTuple) closeTupleMaker();
    if(m_dbg) cout << "TSelector_SusyNtuple_Truth::Terminate" << endl;    
    if(m_writeOut) { out.close(); }    

    char buffer[10];
    sprintf(buffer, "%d", isys);
    printf("%s\n", buffer);
    string str;
    // string str1 = "histos_ZN_tauveto_signal_";
    // str.append(outputfile);
    // str.append(buffer);
    string str2 = ".root";
    str.append(str2);
    cout <<"str= " << str << endl;

    TString outputfilename = str;
    TFile* output_file = new TFile(outputfilename, "recreate") ;//update or recreate?
    cout << "TFile created" << endl;
    output_file->cd();

    /*if(!calcSysUncert) */writeHistos();
    if(calcSysUncert) writeHistos_sysUncert();
    cout << "histos written" << endl;
    output_file->Write() ;
    cout << "outputfile written" << endl;
    output_file->Close();
   cout << "closed" << endl;
}
/*--------------------------------------------------------------------------------*/
// The Terminate() function is the last function to be called
/*--------------------------------------------------------------------------------*/
void TSelector_SusyNtuple_Truth::Terminate()
{
//   closeTupleMaker();
}
void TSelector_SusyNtuple_Truth::defineHistos()
{

}
/*--------------------------------------------------------------------------------*/
void TSelector_SusyNtuple_Truth::defineHistos_sysUncert()
{

}
void TSelector_SusyNtuple_Truth::writeHistos()
{

}
void TSelector_SusyNtuple_Truth::writeHistos_sysUncert()
{

}

void TSelector_SusyNtuple_Truth::calcJet_variables(TLorentzVector met_TLV, TruthJetVector signalJets)
{    
    TLorentzVector signalJet0_TLV;
    TLorentzVector signalJet1_TLV;
    nSignalJets = signalJets.size();

    if(nSignalJets>0){
        signalJet0_TLV.SetPtEtaPhiE(signalJets.at(0)->Pt(), // DG why twice?
                                    signalJets.at(0)->eta,
                                    signalJets.at(0)->phi,
                                    signalJets.at(0)->Pt()*cosh(signalJets.at(0)->eta));
        signalJet0_TLV.SetPtEtaPhiM(signalJets.at(0)->Pt(),
                                    signalJets.at(0)->eta,
                                    signalJets.at(0)->phi,
                                    signalJets.at(0)->m);
    }
    if(nSignalJets>1){
        signalJet1_TLV.SetPtEtaPhiE(signalJets.at(1)->Pt(),
                                    signalJets.at(1)->eta,
                                    signalJets.at(1)->phi,
                                    signalJets.at(1)->Pt()*cosh(signalJets.at(1)->eta));
        signalJet1_TLV.SetPtEtaPhiM(signalJets.at(1)->Pt(),
                                    signalJets.at(1)->eta,
                                    signalJets.at(1)->phi,
                                    signalJets.at(1)->m);
    }
    if(nSignalJets>0){
        pTj0 = signalJet0_TLV.Pt();
        pTj1 = (nSignalJets>1) ? signalJet1_TLV.Pt() : 0.;
        eta_j0 = signalJet0_TLV.Eta();
        eta_j1 = (nSignalJets>1) ? signalJet1_TLV.Eta() : 0.;
        mjj = (nSignalJets>1) ? (signalJet0_TLV + signalJet1_TLV).M() : signalJet0_TLV.M();
        DeltaPhijj = (nSignalJets>1) ? fabs(signalJet0_TLV.DeltaPhi(signalJet1_TLV)) : 0.;
        pTjj = (nSignalJets>1) ? (signalJet0_TLV + signalJet1_TLV).Pt() : signalJet0_TLV.Pt();
        DeltaPhiMETj0 = fabs(signalJet0_TLV.DeltaPhi(met_TLV));
        DeltaPhiMETj1 = (nSignalJets>1) ? fabs(signalJet1_TLV.DeltaPhi(met_TLV)) : 0.;
        DeltaRjj = (nSignalJets>1) ?  fabs(signalJet0_TLV.DeltaR(signalJet1_TLV)) : 0.;
        DeltaEtajj = (nSignalJets>1) ? fabs(signalJet0_TLV.Eta() - signalJet1_TLV.Eta()) : 0.;
        DeltaYjj = (nSignalJets>1) ? fabs(signalJet0_TLV.Rapidity() - signalJet1_TLV.Rapidity()) : 0.;
        DeltaPhiMETjj = (nSignalJets>1) ?
            fabs((signalJet0_TLV+ signalJet1_TLV).DeltaPhi(met_TLV)) :
            fabs(signalJet0_TLV.DeltaPhi(met_TLV));
//     NBJets = SusyNtAna::numberOfCBJets(signalJets);
//     NCJets = SusyNtAna::numberOfCLJets(signalJets);
//     NFJets = SusyNtAna::numberOfFJets(signalJets);
//     meff = calcMeff(met_TLV, signalJets);
    }
}
void TSelector_SusyNtuple_Truth::calc_EE_variables(TruthParticle* el0, TruthParticle* el1,
                                                   TLorentzVector el0_TLV, TLorentzVector el1_TLV,
                                                   TruthJetVector signalJets,
                                                   TLorentzVector signalJet0_TLV,
                                                   TLorentzVector signalJet1_TLV,
                                                   TLorentzVector met_TLV){

//   pTl0_EE = el0_TLV.Pt();
//   pTl1_EE = el1_TLV.Pt();
//   etal0_EE = el0_TLV.Eta();
//   etal1_EE = el1_TLV.Eta();
//   DeltaR_EE = fabs(el0_TLV.DeltaR(el1_TLV));
//   pTll_EE = (el0_TLV + el1_TLV).Pt();
//   Mll_EE = Mll(el0, el1);
//   METrel_EE = recalcMetRel(met_TLV, el0_TLV, el1_TLV, m_signalJets2Lep, useForwardJets);
//   MET_EE = met_TLV.Pt();
    HT_EE = calcHT(el0_TLV, el1_TLV, met_TLV, signalJets);
//   mTWW_EE = calcMt((el0_TLV + el1_TLV), met_TLV);
//   mT_EE = calcMt((el0_TLV+el1_TLV), met_TLV);
//   mTmin_EE = (calcMt(el0_TLV, met_TLV) > calcMt(el1_TLV, met_TLV)) ? calcMt(el1_TLV, met_TLV) : calcMt(el0_TLV, met_TLV);
    mTmax_EE = (calcMt(el0_TLV, met_TLV) < calcMt(el1_TLV, met_TLV)) ?
        calcMt(el1_TLV, met_TLV) :
        calcMt(el0_TLV, met_TLV);
//   mTl0MET_EE = calcMt(el0_TLV, met_TLV);
//   mTl1MET_EE = calcMt(el1_TLV, met_TLV);
//   mMET_EE = (el0_TLV + el1_TLV + met_TLV).M();
//   DeltaPhi_EE = fabs(el0_TLV.DeltaPhi(el1_TLV));
//   DeltaPhiMETl0_EE = fabs(el0_TLV.DeltaPhi(met_TLV));
//   DeltaPhiMETl1_EE = fabs(el1_TLV.DeltaPhi(met_TLV));
//   DeltaPhiMETll_EE = fabs((el0_TLV + el1_TLV).DeltaPhi(met_TLV));

    Mljj_EE = -1.;
    Mlj_EE = -1.;
    if(nSignalJets>1){
        //find dijet axis:
        double DeltaRDijetEl0 = el0_TLV.DeltaR(signalJet0_TLV + signalJet1_TLV);
        double DeltaRDijetEl1 = el1_TLV.DeltaR(signalJet0_TLV + signalJet1_TLV);
        TLorentzVector closestElecDijetAxis_TLV = (DeltaRDijetEl0 > DeltaRDijetEl1) ? el1_TLV : el0_TLV;
        Mljj_EE = (signalJet0_TLV + signalJet1_TLV + closestElecDijetAxis_TLV).M();
    }
    if(signalJets.size()>0){
        //find dijet axis:
        double DeltaRDijetEl0 = el0_TLV.DeltaR(signalJet0_TLV);
        double DeltaRDijetEl1 = el1_TLV.DeltaR(signalJet0_TLV);
        TLorentzVector closestElecDijetAxis_TLV = (DeltaRDijetEl0 > DeltaRDijetEl1) ? el1_TLV : el0_TLV;
        Mlj_EE = (signalJet0_TLV + closestElecDijetAxis_TLV).M();
    }
    //get all electrons in SusyNtuple which have pT > 6 GeV

    TruthParticleVector Electrons_all_vec;
    for(uint index=0; index<nt.tpr()->size(); ++index){
        TruthParticle* particle = & nt.tpr()->at(index);
        if(fabs(particle->pdgId) == 11) {
            if( particle->Pt() < 6.0 || fabs(particle->Eta()) > 2.47 ) continue;
            Electrons_all_vec.push_back(particle);
        }
    }

    mllZcandImpact_EE = -1.;
    mTllZcandImpact_EE = -1.;
    IClZcandImpact_EE = -5;
    pTlZcandImpact_EE = -1.;
    etalZcandImpact_EE = -1.;
    ptcone30lZcandImpact_EE = -1.;
    etcone30lZcandImpact_EE = -1.;
    d0SiglZcandImpact_EE = -1.;
    z0SinThetalZcandImpact_EE = -1.;

    ZcandLep_exists_EE = false;
    ZcandLep_passesPT_EE = true;
    ZcandLep_passesEta_EE = true;
    ZcandLep_passesPTcone_EE = true;
    ZcandLep_passesETcone_EE = true;
    ZcandLep_passesD0_EE = true;
    ZcandLep_passesZ0_EE = true;
    ZcandLep_PassesMedium_EE = true;
    ZcandLep_PassesTight_EE = true;
    ZcandLep_PassesORAndMllCut_EE = false;
    ZcandLep_PassesPR_EE = true;

    double DeltaMZ_lZcandImpact = 99999.;

    TruthParticleVector Electron_ZcandImpact_vec;
    TruthParticle* el_ZcandImpact_lost;

    TruthParticle* closest_signal_el;
    TLorentzVector closest_signal_el_TLV;
    bool foundCandidate = false;
    for(uint ie=0; ie<Electrons_all_vec.size(); ie++){
        TruthParticle* el_ZcandImpact = Electrons_all_vec.at(ie);        
        if((el_ZcandImpact->DeltaR(*el0) < 0.05) ||
           (el_ZcandImpact->DeltaR(*el1) < 0.05)) continue; //no overlap w/ signal lepton

//       if(fabs(el_ZcandImpact->d0Sig(true)) >= ELECTRON_D0SIG_CUT_WH) continue;
//       if(fabs(el_ZcandImpact->z0SinTheta(true)) >= ELECTRON_Z0_SINTHETA_CUT) continue;

        TLorentzVector ZcandImpactElec_TLV;
        ZcandImpactElec_TLV.SetPtEtaPhiE(el_ZcandImpact->Pt(),
                                         el_ZcandImpact->eta,
                                         el_ZcandImpact->phi,
                                         el_ZcandImpact->Pt()*cosh(el_ZcandImpact->eta));
        ZcandImpactElec_TLV.SetPtEtaPhiM(el_ZcandImpact->Pt(),
                                         el_ZcandImpact->eta,
                                         el_ZcandImpact->phi,
                                         el_ZcandImpact->m);
//       SFOS pair with leading or subleading signal electron closer to Zmass?
        float DeltaMZ_l0_ZcandImpact = 999.;
        if((fabs(MZ - SusyNtAna::Mll(el0, el_ZcandImpact)) < DeltaMZ_lZcandImpact) &&
           ((el_ZcandImpact->charge * el0->charge)<0.)){
            DeltaMZ_l0_ZcandImpact = fabs(MZ - SusyNtAna::Mll(el0, el_ZcandImpact));
            foundCandidate = true;
        }
        float DeltaMZ_l1_ZcandImpact = 999.;
        if((fabs(MZ - SusyNtAna::Mll(el1, el_ZcandImpact)) < DeltaMZ_lZcandImpact) &&
           ((el_ZcandImpact->charge * el1->charge)<0.)){
            DeltaMZ_l1_ZcandImpact = fabs(MZ - SusyNtAna::Mll(el1, el_ZcandImpact));
            foundCandidate = true;
        }
        if(foundCandidate &&
           ((fabs(MZ - SusyNtAna::Mll(el1, el_ZcandImpact)) < DeltaMZ_lZcandImpact) ||
            (fabs(MZ - SusyNtAna::Mll(el0, el_ZcandImpact)) < DeltaMZ_lZcandImpact) )){
            el_ZcandImpact_lost = el_ZcandImpact;
            if(DeltaMZ_l0_ZcandImpact < DeltaMZ_l1_ZcandImpact){
                closest_signal_el = el0;
                closest_signal_el_TLV = el0_TLV;
            }else{
                closest_signal_el = el1;
                closest_signal_el_TLV = el1_TLV;
            }
            mllZcandImpact_EE = SusyNtAna::Mll(closest_signal_el, el_ZcandImpact);
            DeltaMZ_lZcandImpact = fabs(MZ - SusyNtAna::Mll(closest_signal_el, el_ZcandImpact));            
            Electron_ZcandImpact_vec.push_back(el_ZcandImpact);
        }
    }
}
//----------------------------------------------------------
void TSelector_SusyNtuple_Truth::calc_MM_variables(TruthParticle* mu0, TruthParticle* mu1,
                                                   TLorentzVector mu0_TLV, TLorentzVector mu1_TLV,
                                                   TruthJetVector signalJets,
                                                   TLorentzVector signalJet0_TLV,
                                                   TLorentzVector signalJet1_TLV,
                                                   TLorentzVector met_TLV){
    HT_MM = calcHT(mu0_TLV, mu1_TLV, met_TLV, signalJets);
    mTmax_MM = (calcMt(mu0_TLV, met_TLV) < calcMt(mu1_TLV, met_TLV)) ?
        calcMt(mu1_TLV, met_TLV) : calcMt(mu0_TLV, met_TLV);
    Mljj_MM = -1.;
    Mlj_MM = -1.;
    if(nSignalJets>1){
        double DeltaRDijetMu0 = mu0_TLV.DeltaR(signalJet0_TLV + signalJet1_TLV);
        double DeltaRDijetMu1 = mu1_TLV.DeltaR(signalJet0_TLV + signalJet1_TLV);
        TLorentzVector closestMuonDijetAxis_TLV = (DeltaRDijetMu0 > DeltaRDijetMu1) ? mu1_TLV : mu0_TLV;
        Mljj_MM = (signalJet0_TLV + signalJet1_TLV + closestMuonDijetAxis_TLV).M();
    }
    if(nSignalJets>0){
        double DeltaRDijetMu0 = mu0_TLV.DeltaR(signalJet0_TLV);
        double DeltaRDijetMu1 = mu1_TLV.DeltaR(signalJet0_TLV);
        TLorentzVector closestMuonDijetAxis_TLV = (DeltaRDijetMu0 > DeltaRDijetMu1) ? mu1_TLV : mu0_TLV;
        Mlj_MM = (signalJet0_TLV + closestMuonDijetAxis_TLV).M();
    }
    DeltaEtall_MM = fabs(mu0_TLV.Eta() - mu1_TLV.Eta());
    TruthParticleVector Muons_all_vec;
    for(uint index=0; index<nt.tpr()->size(); ++index){
        TruthParticle* particle = & nt.tpr()->at(index);
        if(fabs(particle->pdgId) == 13) {
            if( particle->Pt() < 6.0 || fabs(particle->Eta()) > 2.4 ) continue;
            Muons_all_vec.push_back(particle);
        }
    }

    //ZcandImpact muons: all muons, only check for distance to signal muons
    mllZcandImpact_MM = -1.;
    mTllZcandImpact_MM = -1.;
    IClZcandImpact_MM = -5;
    pTlZcandImpact_MM = -1.;
    etalZcandImpact_MM = -1.;
    ptcone30lZcandImpact_MM = -1.;
    etcone30lZcandImpact_MM = -1.;
    d0SiglZcandImpact_MM = -1.;
    z0SinThetalZcandImpact_MM = -1.;

    ZcandLep_exists_MM = false;
    ZcandLep_passesPT_MM = true;
    ZcandLep_passesEta_MM = true;
    ZcandLep_passesPTcone_MM = true;
    ZcandLep_passesETcone_MM = true;
    ZcandLep_passesD0_MM = true;
    ZcandLep_passesZ0_MM = true;
    ZcandLep_PassesMedium_MM = true;
    ZcandLep_PassesTight_MM = true;
    ZcandLep_PassesORAndMllCut_MM = false;
    ZcandLep_PassesPR_MM = true;

    double DeltaMZ_lZcandImpact = 99999.;
    TruthParticleVector Muon_ZcandImpact_vec;
    TruthParticle* mu_ZcandImpact_lost;
    TruthParticle* closest_signal_mu;
    TLorentzVector closest_signal_mu_TLV;
    bool foundCandidate = false;
    for(uint im=0; im<Muons_all_vec.size(); im++){
        TruthParticle* mu_ZcandImpact = Muons_all_vec.at(im);
        if((mu_ZcandImpact->DeltaR(*mu0) < 0.05) ||
           (mu_ZcandImpact->DeltaR(*mu1) < 0.05)) continue; //only check for separation of signal leptons

        TLorentzVector ZcandImpact_TLV;
        ZcandImpact_TLV.SetPtEtaPhiE(mu_ZcandImpact->Pt(),
                                     mu_ZcandImpact->eta,
                                     mu_ZcandImpact->phi,
                                     mu_ZcandImpact->Pt()*cosh(mu_ZcandImpact->eta));
        ZcandImpact_TLV.SetPtEtaPhiM(mu_ZcandImpact->Pt(),
                                     mu_ZcandImpact->eta,
                                     mu_ZcandImpact->phi,
                                     mu_ZcandImpact->m);
        float DeltaMZ_l0_ZcandImpact = 999.;
        if(((fabs(MZ - SusyNtAna::Mll(mu0, mu_ZcandImpact)) < DeltaMZ_lZcandImpact) &&
            ((mu_ZcandImpact->charge * mu0->charge)<0.))){
            DeltaMZ_l0_ZcandImpact = fabs(MZ - SusyNtAna::Mll(mu0, mu_ZcandImpact));
            foundCandidate = true;
        }
        float DeltaMZ_l1_ZcandImpact = 999.;
        if(((fabs(MZ - SusyNtAna::Mll(mu1, mu_ZcandImpact)) < DeltaMZ_lZcandImpact) &&
            ((mu_ZcandImpact->charge * mu1->charge)<0.))){
            DeltaMZ_l1_ZcandImpact = fabs(MZ - SusyNtAna::Mll(mu1, mu_ZcandImpact));
            foundCandidate = true;
        }
        if(foundCandidate &&
           ((fabs(MZ - SusyNtAna::Mll(mu1, mu_ZcandImpact)) < DeltaMZ_lZcandImpact) ||
            (fabs(MZ - SusyNtAna::Mll(mu0, mu_ZcandImpact)) < DeltaMZ_lZcandImpact))){
            mu_ZcandImpact_lost = mu_ZcandImpact;
            if(DeltaMZ_l0_ZcandImpact < DeltaMZ_l1_ZcandImpact){
                closest_signal_mu = mu0;
                closest_signal_mu_TLV = mu0_TLV;
            }else{
                closest_signal_mu = mu1;
                closest_signal_mu_TLV = mu1_TLV;
            }
            mllZcandImpact_MM = SusyNtAna::Mll(closest_signal_mu, mu_ZcandImpact);
            DeltaMZ_lZcandImpact = fabs(MZ - SusyNtAna::Mll(closest_signal_mu, mu_ZcandImpact));
            Muon_ZcandImpact_vec.push_back(mu_ZcandImpact);
        }        
    }
}
//----------------------------------------------------------
void TSelector_SusyNtuple_Truth::calc_EM_variables(TruthParticle* mu, TruthParticle* el,
                                                   TLorentzVector mu_TLV, TLorentzVector el_TLV,
                                                   TruthJetVector signalJets,
                                                   TLorentzVector signalJet0_TLV,
                                                   TLorentzVector signalJet1_TLV,
                                                   TLorentzVector met_TLV){
    HT_EM = calcHT(mu_TLV, el_TLV, met_TLV, signalJets);
    mTmax_EM = (calcMt(mu_TLV, met_TLV) < calcMt(el_TLV, met_TLV)) ?
        calcMt(el_TLV, met_TLV) : calcMt(mu_TLV, met_TLV);
    Mljj_EM = -1.;
    Mlj_EM = -1.;
    if(nSignalJets>1){
        //find dijet axis:
        double DeltaRDijetMu = mu_TLV.DeltaR(signalJet0_TLV + signalJet1_TLV);
        double DeltaRDijetEl = el_TLV.DeltaR(signalJet0_TLV + signalJet1_TLV);
        TLorentzVector closestLepDijetAxis_TLV = (DeltaRDijetMu > DeltaRDijetEl) ? el_TLV : mu_TLV;        
        Mljj_EM = (signalJet0_TLV + signalJet1_TLV + closestLepDijetAxis_TLV).M();
    }
    if(nSignalJets>0){
        //find dijet axis:
        double DeltaRDijetMu = mu_TLV.DeltaR(signalJet0_TLV);
        double DeltaRDijetEl = el_TLV.DeltaR(signalJet0_TLV);
        TLorentzVector closestLepDijetAxis_TLV = (DeltaRDijetMu > DeltaRDijetEl) ? el_TLV : mu_TLV;
        Mlj_EM = (signalJet0_TLV + closestLepDijetAxis_TLV).M();
    }
    DeltaEtall_EM = fabs(mu_TLV.Eta() - el_TLV.Eta());    
    TruthParticleVector Muons_all_vec;
    for(uint index=0; index<nt.tpr()->size(); ++index){
        TruthParticle* particle = & nt.tpr()->at(index);
        if(fabs(particle->pdgId) == 13) {
            if( particle->Pt() < 6.0 || fabs(particle->Eta()) > 2.4 ) continue;
            Muons_all_vec.push_back(particle);
        }
    }
    TruthParticleVector Electrons_all_vec;
    for(uint index=0; index<nt.tpr()->size(); ++index){
        TruthParticle* particle = & nt.tpr()->at(index);
        if(fabs(particle->pdgId) == 11) {
            if( particle->Pt() < 6.0 || fabs(particle->Eta()) > 2.47 ) continue;
            Electrons_all_vec.push_back(particle);
        }
    }

    mllZcandImpact_mu_EM = -1.;
    mTllZcandImpact_mu_EM = -1.;
    IClZcandImpact_mu_EM = -5;
    pTlZcandImpact_mu_EM = -1.;
    etalZcandImpact_mu_EM = -1.;
    ptcone30lZcandImpact_mu_EM = -1.;
    etcone30lZcandImpact_mu_EM = -1.;
    d0SiglZcandImpact_mu_EM = -1.;
    z0SinThetalZcandImpact_mu_EM = -1.;

    mllZcandImpact_el_EM = -1.;
    mTllZcandImpact_el_EM = -1.;
    IClZcandImpact_el_EM = -5;
    pTlZcandImpact_el_EM = -1.;
    etalZcandImpact_el_EM = -1.;
    ptcone30lZcandImpact_el_EM = -1.;
    etcone30lZcandImpact_el_EM = -1.;
    d0SiglZcandImpact_el_EM = -1.;
    z0SinThetalZcandImpact_el_EM = -1.;

    double DeltaMZ_lZcandImpact_mu = 99999.;
    double DeltaMZ_lZcandImpact_el = 99999.;
    TruthParticleVector Muon_ZcandImpact_vec;
    TruthParticleVector Electron_ZcandImpact_vec;
    TruthParticle* mu_ZcandImpact_lost=NULL;
    TruthParticle* el_ZcandImpact_lost=NULL;

    bool isMu = false;
    bool isEl = false;

    TruthParticle* mu_ZcandImpact;
    //loop over all muons and electrons, check for separation from
    //signal leptons (in the meaning of DeltaR). If one SFOS pair is
    //closer to Zmass, use this third lepton. Mark if it is electron
    //or muon with 'isEl' and 'isMu'
    for(uint im=0; im<Muons_all_vec.size(); im++){
        mu_ZcandImpact = Muons_all_vec.at(im);
        if((mu_ZcandImpact->DeltaR(*mu) < 0.05) ||
           mu_ZcandImpact->DeltaR(*el) < 0.05) continue;  //only check for separation of signal leptons        
        if((mu_ZcandImpact->charge * mu->charge)<0.){
            if((fabs(MZ - SusyNtAna::Mll(mu, mu_ZcandImpact)) < DeltaMZ_lZcandImpact_mu)){
                DeltaMZ_lZcandImpact_mu = fabs(MZ - SusyNtAna::Mll(mu, mu_ZcandImpact));
                mu_ZcandImpact_lost = mu_ZcandImpact;
                Muon_ZcandImpact_vec.push_back(mu_ZcandImpact);
            }
        }        
    }

    TruthParticle* el_ZcandImpact;
    for(uint ie=0; ie<Electrons_all_vec.size(); ie++){
        el_ZcandImpact = Electrons_all_vec.at(ie);
        if((el_ZcandImpact->DeltaR(*mu) < 0.05) ||
           (el_ZcandImpact->DeltaR(*el) < 0.05)) continue; //only check for separation of signal leptons
        if((el_ZcandImpact->charge * el->charge)<0.){
            if(fabs(MZ - SusyNtAna::Mll(el, el_ZcandImpact)) < DeltaMZ_lZcandImpact_el){
                DeltaMZ_lZcandImpact_el = fabs(MZ - SusyNtAna::Mll(el, el_ZcandImpact));
                el_ZcandImpact_lost = el_ZcandImpact;
                Electron_ZcandImpact_vec.push_back(el_ZcandImpact);
            }
        }
    }

    if(DeltaMZ_lZcandImpact_mu < DeltaMZ_lZcandImpact_el) isMu = true;
    else if (DeltaMZ_lZcandImpact_mu > DeltaMZ_lZcandImpact_el) isEl = true;

    if(isMu){ mllZcandImpact_mu_EM = SusyNtAna::Mll(mu, mu_ZcandImpact_lost); }
    if(isEl){ mllZcandImpact_el_EM = SusyNtAna::Mll(el, el_ZcandImpact_lost); }
}

void TSelector_SusyNtuple_Truth::fillHistos_EE(int cutnumber, float weight) { }
void TSelector_SusyNtuple_Truth::fillHistos_MM(int cutnumber, float weight) { }
void TSelector_SusyNtuple_Truth::fillHistos_EM(int cutnumber, float weight) { }
void TSelector_SusyNtuple_Truth::fillHistos_EE_SRSS1(float cut_EE, float weight_ALL_EE)
{

}

void TSelector_SusyNtuple_Truth::fillHistos_MM_SRSS1(float cut_MM, float weight_ALL_MM)
{

}

void TSelector_SusyNtuple_Truth::fillHistos_EM_SRSS1(float cut_EM, float weight_ALL_EM)
{

}

void TSelector_SusyNtuple_Truth::print_counters() const
{
    cout<<endl
        <<"Counts summary:"<<endl
        <<"counter_input          : "<<counter_input          <<endl
        <<"counter_event_cleaning : "<<counter_event_cleaning <<endl
        <<"counter_EE             : "<<counter_EE             <<endl
        <<"counter_EM             : "<<counter_EM             <<endl
        <<"counter_MM             : "<<counter_MM             <<endl  
        <<"counter_EE_SRSS1 : "<<counter_EE_SRSS1<<" counter_EE_SRSS2 : "<<counter_EE_SRSS2<<endl
        <<"counter_EM_SRSS1 : "<<counter_EM_SRSS1<<" counter_EM_SRSS2 : "<<counter_EM_SRSS2<<endl
        <<"counter_MM_SRSS1 : "<<counter_MM_SRSS1<<" counter_MM_SRSS2 : "<<counter_MM_SRSS2<<endl
        <<endl;
}
