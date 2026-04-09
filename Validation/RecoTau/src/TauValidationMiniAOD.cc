// Class used to do the Validation of the Tau in miniAOD and reco/hlt
// Created by Aniello Spiezia on August 13, 2019
// Updated By Ece Asilar, Gage DeZoort on April 6th, 2020
// Updated By Gourab Saha on July 4th, 2023
// Updated By Elena Vernazza on April 9th, 2026

#include "Validation/RecoTau/interface/TauValidationMiniAOD.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace edm;
using namespace std;
using namespace reco;

TauValidationMiniAOD::TauValidationMiniAOD(const edm::ParameterSet &iConfig) {
  // Flag for the definition of the tau type (mini or reco):
  TauType = iConfig.getParameter<string>("TauType");
  isMini = (std::string("mini") == TauType);  // <pat::TauCollection>
  isReco = (std::string("reco") == TauType);  // <reco::PFTauCollection>
  // Input collection of legitimate taus:
  if (isMini) {
    patTauToken_ = consumes<pat::TauCollection>(iConfig.getParameter<InputTag>("TauCollection"));
  } else if (isReco) {
    pfTauToken_ = consumes<reco::PFTauCollection>(iConfig.getParameter<InputTag>("TauCollection"));
  } else {
    throw cms::Exception("Configuration") << "Unknown tau type: " << TauType << "\nPlease use 'mini', or 'reco'.";
  }
  // Input collection to compare to taus:
  genRefToken_ = consumes<edm::View<reco::Candidate>>(iConfig.getParameter<InputTag>("RefCollection"));
  // Input collection of gen particles:
  prunedGenToken_ = consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<InputTag>("GenCollection"));
  // Input collection of primary vertices:
  pvToken_ = consumes<std::vector<reco::Vertex>>(iConfig.getParameter<InputTag>("PVCollection"));
  // Information about reference collection:
  extensionName_ = iConfig.getParameter<string>("ExtensionName");
  // List of discriminators and their cuts:
  discriminators_ = iConfig.getParameter<std::vector<edm::ParameterSet>>("discriminators");
  isHLT = iConfig.getParameter<bool>("isHLT");
}

//------------------------------------------------------------------------------
// ~TauValidationMiniAOD
//------------------------------------------------------------------------------
TauValidationMiniAOD::~TauValidationMiniAOD() {}

//------------------------------------------------------------------------------
// bookHistograms
//------------------------------------------------------------------------------

void TauValidationMiniAOD::bookHistograms(DQMStore::IBooker &ibooker, edm::Run const &iRun, edm::EventSetup const &iSetup) {
 
  std::string main_folder;
  if (isHLT) {
    main_folder = "HLT/Tau/TauValidation";
  } else {
    main_folder = "RecoTauV/miniAODValidation/" + extensionName_;
  }

  // ---------------------------- Book Summary Histograms -------------------------------

  ibooker.setCurrentFolder(main_folder + "/Summary");

  // book histograms for taus matched to reference gen particle
  for (auto& hVar : histoVars) {
    auto [nBins, hMin, hMax] = hVar.second;
    h_tausMatchedToRef_[hVar.first] = ibooker.book1D("tau_" + hVar.first, "#tau;" + hVar.first, nBins, hMin, hMax);
  }

  // book deepTau and decay mode histograms
  decayModeTauReco = ibooker.book1D("tau_decayMode_reco", "Reco #tau;Decay Mode", 15, -0.5, 14.5);
  decayModeTauGen = ibooker.book1D("tau_decayMode_gen", "Gen #tau;Decay Mode", 15, -0.5, 14.5);
  dmMigration = ibooker.book2D("dmMigration", ";Gen #tau DM;Reco #tau DM", 15, -0.5, 14.5, 15, -0.5, 14.5);
  nTau_vs_dm = ibooker.book2D("ntau_vs_dm", ";nTau;#tau DM", 15, 0, 15, 15, 0, 15);

  if (isMini) {
    // book histograms for deepTau discriminators and decay mode finding
    DeepTau2018v2p5VSele = ibooker.book1D("tau_byDeepTau2018v2p5VSeraw", ";DeepTau vs Ele", 200, 0., 1.);
    DeepTau2018v2p5VSjet = ibooker.book1D("tau_byDeepTau2018v2p5VSjetraw", ";DeepTau vs Jet", 200, 0., 1.);
    DeepTau2018v2p5VSmuo = ibooker.book1D("tau_byDeepTau2018v2p5VSmuraw", ";DeepTau vs Muo", 200, 0., 1.);
    decayModeFinding = ibooker.book1D("tau_decayModeFinding", ";Decay Mode Finding", 2, -0.5, 1.5);
    
    // book summary histograms for discriminators
    summary_num = ibooker.book1D("summaryPlotNum", "Summary Num", 15, -0.5, 14.5);
    summary_den = ibooker.book1D("summaryPlotDen", "Summary Den", 15, -0.5, 14.5);
    // add discriminator labels to summary plots
    for (uint j = 0; j < discriminators_.size(); j++) {
      string DiscriminatorLabel = discriminators_[j].getParameter<string>("discriminator");
      summary_num->setBinLabel(j + 1, DiscriminatorLabel);
      summary_den->setBinLabel(j + 1, DiscriminatorLabel);
    }
  }
  
  // book histograms for tau pt and prong pt for different decay modes
  h_TauMass_dm_.resize(dm_list.size());
  h_pTOverProng_dm_.resize(dm_list.size());
  std::string pT_Prong = "p_{T}^{#tau};p_{T} lead ch had";
  for (uint i = 0; i < dm_list.size(); i++) {
    const auto& [dm, label] = dm_list[i];;
    h_TauMass_dm_[i] = ibooker.book1D("mtau_" + label, "DM = " + label + ";m_{#tau}", 20, 0.0, 2.0);
    if (isMini)
      h_pTOverProng_dm_[i] = ibooker.book2D("pTOverProng_" + label, "DM = " + label + ";" + pT_Prong, 50, 0, 1000, 50, 0, 1000);
  }

  if (isMini) {
    // ---------------------------- /vsJet/ ---------------------------------------------
    if (extensionName_.compare(qcd) == 0 || extensionName_.compare(real_data) == 0 || extensionName_.compare(ztt) == 0) {

      ibooker.setCurrentFolder(main_folder + "/vsJet/tight");
      for (auto& hVar : histoVars) {
        auto [nBins, hMin, hMax] = hVar.second;
        h_tausMatchedToRef_TightvsJet[hVar.first] = ibooker.book1D("tau_tightvsJet_" + hVar.first, "#tau;" + hVar.first, nBins, hMin, hMax);
      } 

      ibooker.setCurrentFolder(main_folder + "/vsJet/medium");
      for (auto& hVar : histoVars) {
        auto [nBins, hMin, hMax] = hVar.second;
        h_tausMatchedToRef_MediumvsJet[hVar.first] = ibooker.book1D("tau_mediumvsJet_" + hVar.first, "#tau;" + hVar.first, nBins, hMin, hMax);
      } 

      ibooker.setCurrentFolder(main_folder + "/vsJet/loose");
      for (auto& hVar : histoVars) {
        auto [nBins, hMin, hMax] = hVar.second;
        h_tausMatchedToRef_LoosevsJet[hVar.first] = ibooker.book1D("tau_loosevsJet_" + hVar.first, "#tau;" + hVar.first, nBins, hMin, hMax);
      } 
    }

    // ---------------------------- /vsEle/ ---------------------------------------------
    if (extensionName_.compare(real_eledata) == 0 || extensionName_.compare(zee) == 0 || extensionName_.compare(ztt) == 0) {

      ibooker.setCurrentFolder(main_folder + "/vsEle/tight");
      for (auto& hVar : histoVars) {
        auto [nBins, hMin, hMax] = hVar.second;
        h_tausMatchedToRef_TightvsEle[hVar.first] = ibooker.book1D("tau_tightvsEle_" + hVar.first, "#tau;" + hVar.first, nBins, hMin, hMax);
      }

      ibooker.setCurrentFolder(main_folder + "/vsEle/medium");
      for (auto& hVar : histoVars) {
        auto [nBins, hMin, hMax] = hVar.second;
        h_tausMatchedToRef_MediumvsEle[hVar.first] = ibooker.book1D("tau_mediumvsEle_" + hVar.first, "#tau;" + hVar.first, nBins, hMin, hMax);
      }

      ibooker.setCurrentFolder(main_folder + "/vsEle/loose");
      for (auto& hVar : histoVars) {
        auto [nBins, hMin, hMax] = hVar.second;
        h_tausMatchedToRef_LoosevsEle[hVar.first] = ibooker.book1D("tau_loosevsEle_" + hVar.first, "#tau;" + hVar.first, nBins, hMin, hMax);
      }
    }

    // ---------------------------- /vsMuo/ ---------------------------------------------
    if (extensionName_.compare(real_mudata) == 0 || extensionName_.compare(zmm) == 0 || extensionName_.compare(ztt) == 0) {

      ibooker.setCurrentFolder(main_folder + "/vsMuo/tight");
      for (auto& hVar : histoVars) {
        auto [nBins, hMin, hMax] = hVar.second;
        h_tausMatchedToRef_TightvsMuo[hVar.first] = ibooker.book1D("tau_tightvsMuo_" + hVar.first, "#tau;" + hVar.first, nBins, hMin, hMax);
      }

      ibooker.setCurrentFolder(main_folder + "/vsMuo/medium");
      for (auto& hVar : histoVars) {
        auto [nBins, hMin, hMax] = hVar.second;
        h_tausMatchedToRef_MediumvsMuo[hVar.first] = ibooker.book1D("tau_mediumvsMuo_" + hVar.first, "#tau;" + hVar.first, nBins, hMin, hMax);
      }

      ibooker.setCurrentFolder(main_folder + "/vsMuo/loose");
      for (auto& hVar : histoVars) {
        auto [nBins, hMin, hMax] = hVar.second;
        h_tausMatchedToRef_LoosevsMuo[hVar.first] = ibooker.book1D("tau_loosevsMuo_" + hVar.first, "#tau;" + hVar.first, nBins, hMin, hMax);
      }
    }
  }
}

//------------------------------------------------------------------------------
// analyze
//------------------------------------------------------------------------------

void TauValidationMiniAOD::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup) {

  // create a handle to the tau collection
  edm::Handle<pat::TauCollection> patTaus;
  edm::Handle<reco::PFTauCollection> pfTaus;
  if (isMini) {
    iEvent.getByToken(patTauToken_, patTaus);
    if (!patTaus.isValid()) {
      edm::LogWarning("TauValidationMiniAOD") << " Pat Tau collection not found while running TauValidationMiniAOD.cc ";
      return;
    }
  } else if (isReco) {
    iEvent.getByToken(pfTauToken_, pfTaus);
    if (!pfTaus.isValid()) {
      edm::LogWarning("TauValidationMiniAOD") << " Pf Tau collection not found while running TauValidationMiniAOD.cc ";
      return;
    }
  }

  // create a handle to the reference collection
  typedef edm::View<reco::Candidate> refCandidateCollection;
  edm::Handle<refCandidateCollection> genRefParticles;
  bool isRef = iEvent.getByToken(genRefToken_, genRefParticles);
  if (!isRef) {
    edm::LogWarning("TauValidationMiniAOD") << " Reference collection not found while running TauValidationMiniAOD.cc ";
    return;
  }

  // create a handle to the gen Part collection
  edm::Handle<std::vector<reco::GenParticle>> genParticles;
  iEvent.getByToken(prunedGenToken_, genParticles);
  if (!genParticles.isValid()) {
    edm::LogWarning("TauValidationMiniAOD") << " GenParticle collection not found while running TauValidationMiniAOD.cc ";
  }

  // create a handle to the primary vertex collection
  Handle<std::vector<reco::Vertex>> primaryVertices;
  bool isPV = iEvent.getByToken(pvToken_, primaryVertices);
  if (!isPV) {
    edm::LogWarning("TauValidationMiniAOD") << " PV collection not found while running TauValidationMiniAOD.cc ";
  }

  float dR2min = 0.15;
  int matchedTauIndex = -99;
  float gendR2min = 0.15;
  int genmatchedTauIndex = -99;

  // dR match between gen reference object (taus, electrons, jets, muons, ...) and tau
  for (refCandidateCollection::const_iterator genRef = genRefParticles->begin(); genRef != genRefParticles->end(); genRef++) {

    float tau_pt = -1;
    float tau_eta = -99;
    float tau_phi = -99;
    float tau_mass = -1;
    int tau_decayMode = -99;
    uint nTaus = -1;

    // find best matched tau to reference gen particle
    if (isMini) {
      nTaus = patTaus->size();
      for (unsigned iTau = 0; iTau < nTaus; iTau++) {
        pat::TauRef tau(patTaus, iTau);
        float dR2 = deltaR2(tau->eta(), tau->phi(), genRef->eta(), genRef->phi());
        if (dR2 < dR2min) {
          dR2min = dR2;
          matchedTauIndex = iTau;
          tau_pt = tau->pt();
          tau_eta = tau->eta();
          tau_phi = tau->phi();
          tau_mass = tau->mass();
          tau_decayMode = tau->decayMode();
        }
      }
    } else if (isReco) {
      nTaus = pfTaus->size();
      for (unsigned iTau = 0; iTau < nTaus; iTau++) {
        reco::PFTauRef tau(pfTaus, iTau);
        float dR2 = deltaR2(tau->eta(), tau->phi(), genRef->eta(), genRef->phi());
        if (dR2 < dR2min) {
          dR2min = dR2;
          matchedTauIndex = iTau;
          tau_pt = tau->pt();
          tau_eta = tau->eta();
          tau_phi = tau->phi();
          tau_mass = tau->mass();
          tau_decayMode = tau->decayMode();
        }
      }
    }

    if (dR2min < 0.15) {

      // fill kinematics with matched Tau quantities
      h_tausMatchedToRef_["pt"]->Fill(tau_pt);
      h_tausMatchedToRef_["eta"]->Fill(tau_eta);
      h_tausMatchedToRef_["phi"]->Fill(tau_phi);
      h_tausMatchedToRef_["mass"]->Fill(tau_mass);
      h_tausMatchedToRef_["pu"]->Fill(primaryVertices->size());
      decayModeTauReco->Fill(tau_decayMode);

      // fill decay mode population plot
      nTau_vs_dm->Fill(nTaus, tau_decayMode);

      // fill tau mass for decay modes 0,1+2,5,6,7,10,11
      for (size_t i = 0; i < dm_list.size(); ++i) {
          const auto& dms = dm_list[i].first;
          if (std::find(dms.begin(), dms.end(), tau_decayMode) != dms.end()) {
              h_TauMass_dm_[i]->Fill(tau_mass);
          }
      }

      // find corresponding tau gen particle
      for (unsigned iGen = 0; iGen < genParticles->size(); iGen++) {
        const reco::GenParticle &gentau = genParticles->at(iGen);
        if (abs(gentau.pdgId()) == 15) {
          float gendR2 = deltaR2(tau_eta, tau_phi, gentau.eta(), gentau.phi());
          if (gendR2 < gendR2min) {
            gendR2min = gendR2;
            genmatchedTauIndex = iGen;
          }
        }
      }

      // fill decay mode migration 2D histogragms
      if (gendR2min < 0.15) {
        int nPi0s = 0;
        int nPis = 0;
        const reco::GenParticle &gentau = genParticles->at(genmatchedTauIndex);
        for (unsigned idtr = 0; idtr < gentau.numberOfDaughters(); idtr++) {
          const reco::GenParticle *dtr = dynamic_cast<const reco::GenParticle *>(gentau.daughter(idtr));
          int dtrpdgID = std::abs(dtr->pdgId());
          int dtrstatus = dtr->status();
          if (dtrpdgID == 12 || dtrpdgID == 14 || dtrpdgID == 16) // neutrinos
            continue;
          if (dtrpdgID == 111 || dtrpdgID == 311) // neutral pions and kaons
            nPi0s++;
          else if (dtrpdgID == 211 || dtrpdgID == 321) // charged pions and kaons
            nPis++;
          else if (dtrpdgID == 15 && dtrstatus == 2) {
            for (unsigned idtr2 = 0; idtr2 < dtr->numberOfDaughters(); idtr2++) {
              const reco::GenParticle *dtr2 = dynamic_cast<const reco::GenParticle *>(dtr->daughter(idtr2));
              int dtr2pdgID = std::abs(dtr2->pdgId());
              if (dtr2pdgID == 12 || dtr2pdgID == 14 || dtr2pdgID == 16) // neutrinos
                continue;
              if (dtr2pdgID == 111 || dtr2pdgID == 311) // neutral pions and kaons
                nPi0s++;
              else if (dtr2pdgID == 211 || dtr2pdgID == 321) // charged pions and kaons
                nPis++;
            }
          }
        }
        int genTau_decayMode = findDecayMode(nPis, nPi0s);
        decayModeTauGen->Fill(genTau_decayMode);
        dmMigration->Fill(genTau_decayMode, tau_decayMode);
      }

      // fill histograms for pat tau collection variables
      if (isMini) {
        pat::TauRef matchedTau(patTaus, matchedTauIndex);

        // fill select discriminators with matched Tau quantities
        if (matchedTau->isTauIDAvailable("decayModeFindingNewDMs"))
          decayModeFinding->Fill(matchedTau->tauID("decayModeFindingNewDMs"));
        if (matchedTau->isTauIDAvailable("byDeepTau2018v2p5VSeraw"))
          DeepTau2018v2p5VSele->Fill(matchedTau->tauID("byDeepTau2018v2p5VSeraw"));
        if (matchedTau->isTauIDAvailable("byDeepTau2018v2p5VSjetraw"))
          DeepTau2018v2p5VSjet->Fill(matchedTau->tauID("byDeepTau2018v2p5VSjetraw"));
        if (matchedTau->isTauIDAvailable("byDeepTau2018v2p5VSmuraw"))
          DeepTau2018v2p5VSmuo->Fill(matchedTau->tauID("byDeepTau2018v2p5VSmuraw"));

        // fill tau mass for decay modes 0,1+2,5,6,7,10,11
        float tau_ptLeadChargedCand = matchedTau->ptLeadChargedCand();
        for (size_t i = 0; i < dm_list.size(); ++i) {
            const auto& dms = dm_list[i].first;
            if (std::find(dms.begin(), dms.end(), tau_decayMode) != dms.end()) {
                h_pTOverProng_dm_[i]->Fill(tau_pt, tau_ptLeadChargedCand);
            }
        }

        // count number of taus passing each discriminator's selection cut
        for (uint j = 0; j < discriminators_.size(); j++) {
          string currentDiscriminator = discriminators_[j].getParameter<string>("discriminator");
          double selectionCut = discriminators_[j].getParameter<double>("selectionCut");
          summary_den->Fill(j);
          if (matchedTau->isTauIDAvailable(currentDiscriminator) &&
              matchedTau->tauID(currentDiscriminator) >= selectionCut)
            summary_num->Fill(j);
        }

        // fill the discriminator histograms
        if (extensionName_.compare(qcd) == 0 || extensionName_.compare(real_data) == 0 ||
            extensionName_.compare(ztt) == 0) { // VS JET
          // vsJet/tight
          if (matchedTau->isTauIDAvailable("byTightDeepTau2018v2p5VSjet") &&
              matchedTau->tauID("byTightDeepTau2018v2p5VSjet") >= 0.5) {
            h_tausMatchedToRef_TightvsJet["pt"]->Fill(tau_pt);
            h_tausMatchedToRef_TightvsJet["eta"]->Fill(tau_eta);
            h_tausMatchedToRef_TightvsJet["phi"]->Fill(tau_phi);
            h_tausMatchedToRef_TightvsJet["mass"]->Fill(tau_mass);
            h_tausMatchedToRef_TightvsJet["pu"]->Fill(primaryVertices->size());
          }
          // vsJet/medium
          if (matchedTau->isTauIDAvailable("byMediumDeepTau2018v2p5VSjet") &&
              matchedTau->tauID("byMediumDeepTau2018v2p5VSjet") >= 0.5) {
            h_tausMatchedToRef_MediumvsJet["pt"]->Fill(tau_pt);
            h_tausMatchedToRef_MediumvsJet["eta"]->Fill(tau_eta);
            h_tausMatchedToRef_MediumvsJet["phi"]->Fill(tau_phi);
            h_tausMatchedToRef_MediumvsJet["mass"]->Fill(tau_mass);
            h_tausMatchedToRef_MediumvsJet["pu"]->Fill(primaryVertices->size());
          }
          // vsJet/loose
          if (matchedTau->isTauIDAvailable("byLooseDeepTau2018v2p5VSjet") &&
              matchedTau->tauID("byLooseDeepTau2018v2p5VSjet") >= 0.5) {
            h_tausMatchedToRef_LoosevsJet["pt"]->Fill(tau_pt);
            h_tausMatchedToRef_LoosevsJet["eta"]->Fill(tau_eta);
            h_tausMatchedToRef_LoosevsJet["phi"]->Fill(tau_phi);
            h_tausMatchedToRef_LoosevsJet["mass"]->Fill(tau_mass);
            h_tausMatchedToRef_LoosevsJet["pu"]->Fill(primaryVertices->size());
          }
        }

        if (extensionName_.compare(real_eledata) == 0 || extensionName_.compare(zee) == 0 ||
            extensionName_.compare(ztt) == 0) {  // VS ELE
          // vsEle/tight
          if (matchedTau->isTauIDAvailable("byTightDeepTau2018v2p5VSe") &&
              matchedTau->tauID("byTightDeepTau2018v2p5VSe") >= 0.5) {
            h_tausMatchedToRef_TightvsEle["pt"]->Fill(tau_pt);
            h_tausMatchedToRef_TightvsEle["eta"]->Fill(tau_eta);
            h_tausMatchedToRef_TightvsEle["phi"]->Fill(tau_phi);
            h_tausMatchedToRef_TightvsEle["mass"]->Fill(tau_mass);
            h_tausMatchedToRef_TightvsEle["pu"]->Fill(primaryVertices->size());
          }
          // vsEle/medium
          if (matchedTau->isTauIDAvailable("byMediumDeepTau2018v2p5VSe") &&
              matchedTau->tauID("byMediumDeepTau2018v2p5VSe") >= 0.5) {
            h_tausMatchedToRef_MediumvsEle["pt"]->Fill(tau_pt);
            h_tausMatchedToRef_MediumvsEle["eta"]->Fill(tau_eta);
            h_tausMatchedToRef_MediumvsEle["phi"]->Fill(tau_phi);
            h_tausMatchedToRef_MediumvsEle["mass"]->Fill(tau_mass);
            h_tausMatchedToRef_MediumvsEle["pu"]->Fill(primaryVertices->size());
          }
          // vsEle/loose
          if (matchedTau->isTauIDAvailable("byLooseDeepTau2018v2p5VSe") &&
              matchedTau->tauID("byLooseDeepTau2018v2p5VSe") >= 0.5) {
            h_tausMatchedToRef_LoosevsEle["pt"]->Fill(tau_pt);
            h_tausMatchedToRef_LoosevsEle["eta"]->Fill(tau_eta);
            h_tausMatchedToRef_LoosevsEle["phi"]->Fill(tau_phi);
            h_tausMatchedToRef_LoosevsEle["mass"]->Fill(tau_mass);
            h_tausMatchedToRef_LoosevsEle["pu"]->Fill(primaryVertices->size());
          }
        }

        if (extensionName_.compare(real_mudata) == 0 || extensionName_.compare(zmm) == 0 ||
            extensionName_.compare(ztt) == 0) {  // VS MU
          // vsMuo/tight
          if (matchedTau->isTauIDAvailable("byTightDeepTau2018v2p5VSmu") &&
              matchedTau->tauID("byTightDeepTau2018v2p5VSmu") >= 0.5) {
            h_tausMatchedToRef_TightvsMuo["pt"]->Fill(tau_pt);
            h_tausMatchedToRef_TightvsMuo["eta"]->Fill(tau_eta);
            h_tausMatchedToRef_TightvsMuo["phi"]->Fill(tau_phi);
            h_tausMatchedToRef_TightvsMuo["mass"]->Fill(tau_mass);
            h_tausMatchedToRef_TightvsMuo["pu"]->Fill(primaryVertices->size());
          }
          // vsMuo/medium
          if (matchedTau->isTauIDAvailable("byMediumDeepTau2018v2p5VSmu") &&
              matchedTau->tauID("byMediumDeepTau2018v2p5VSmu") >= 0.5) {
            h_tausMatchedToRef_MediumvsMuo["pt"]->Fill(tau_pt);
            h_tausMatchedToRef_MediumvsMuo["eta"]->Fill(tau_eta);
            h_tausMatchedToRef_MediumvsMuo["phi"]->Fill(tau_phi);
            h_tausMatchedToRef_MediumvsMuo["mass"]->Fill(tau_mass);
            h_tausMatchedToRef_MediumvsMuo["pu"]->Fill(primaryVertices->size());
          }
          // vsMuo/loose
          if (matchedTau->isTauIDAvailable("byLooseDeepTau2018v2p5VSmu") &&
              matchedTau->tauID("byLooseDeepTau2018v2p5VSmu") >= 0.5) {
            h_tausMatchedToRef_LoosevsMuo["pt"]->Fill(tau_pt);
            h_tausMatchedToRef_LoosevsMuo["eta"]->Fill(tau_eta);
            h_tausMatchedToRef_LoosevsMuo["phi"]->Fill(tau_phi);
            h_tausMatchedToRef_LoosevsMuo["mass"]->Fill(tau_mass);
            h_tausMatchedToRef_LoosevsMuo["pu"]->Fill(primaryVertices->size());
          }
        }
      }

    }
  }
}

//------------------------------------------------------------------------------
// fill description
//------------------------------------------------------------------------------

void TauValidationMiniAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  // Default tau validation offline
  desc.add<std::string>("TauType", "mini");
  desc.add<edm::InputTag>("TauCollection", edm::InputTag("slimmedTaus"));
  desc.add<edm::InputTag>("RefCollection", edm::InputTag("kinematicSelectedTauValDenominatorZTT"));
  desc.add<edm::InputTag>("PVCollection", edm::InputTag("offlineSlimmedPrimaryVertices"));
  desc.add<edm::InputTag>("GenCollection", edm::InputTag("prunedGenParticles"));
  desc.add<std::string>("ExtensionName", "ZTT");
  desc.add<bool>("isHLT", false);

  edm::ParameterSetDescription discrDesc;
  discrDesc.add<std::string>("discriminator");
  discrDesc.add<double>("selectionCut");

  desc.addVPSet("discriminators", discrDesc);
  
  descriptions.addWithDefaultLabel(desc);
}
