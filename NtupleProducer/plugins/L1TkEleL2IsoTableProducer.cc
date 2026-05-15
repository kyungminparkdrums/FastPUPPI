// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"

#include "L1Trigger/Phase2L1ParticleFlow/interface/L1TPFUtils.h"

#include "DataFormats/L1TCorrelator/interface/TkElectron.h"

#include <algorithm>

class L1TkEleL2IsoTableProducer : public edm::global::EDProducer<>  {
public:
    explicit L1TkEleL2IsoTableProducer(const edm::ParameterSet&);
    ~L1TkEleL2IsoTableProducer();

private:
    virtual void produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const override;
    
    StringCutObjectSelector<reco::Candidate> sel_;

    struct ExtraVar {
        std::string name, expr;
        StringObjectFunction<reco::Candidate> func;
        ExtraVar(const std::string & n, const std::string & expr) : name(n), expr(expr), func(expr, true) {}
    };
    std::vector<ExtraVar> extraVars_;

    struct CandRecord {
    public:
        std::string coll;
        edm::EDGetTokenT<reco::CandidateView> src;
        StringCutObjectSelector<reco::Candidate> sel;

        CandRecord(const std::string & name, const edm::EDGetTokenT<reco::CandidateView> & tag, const edm::ParameterSet & pset) :
            coll(name), src(tag),
            sel(pset.existsAs<std::string>(name+"_sel") ? pset.getParameter<std::string>(name+"_sel") : "", true) {}
    };

    // I should put this under private but boh
    std::vector<CandRecord> pf_cands_;
    std::vector<CandRecord> tkele_cands_;
};

L1TkEleL2IsoTableProducer::L1TkEleL2IsoTableProducer(const edm::ParameterSet& iConfig) :
    sel_(iConfig.getParameter<std::string>("commonSel"), true)
{

    produces<nanoaod::FlatTable>(); 
    // I should take these from config but boh
    pf_cands_.emplace_back("L1PFCands", consumes<reco::CandidateView>(edm::InputTag("l1tLayer1:PF")), iConfig);
    tkele_cands_.emplace_back("TkEleL2", consumes<reco::CandidateView>(edm::InputTag("l1tLayer2EG:L1CtTkElectron")), iConfig);
}

L1TkEleL2IsoTableProducer::~L1TkEleL2IsoTableProducer() { }

// ------------ method called for each event  ------------
void
L1TkEleL2IsoTableProducer::produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const
{
    edm::Handle<reco::CandidateView> src;
    std::vector<const reco::Candidate *> pf_selected;
    std::vector<const reco::Candidate *> tkele_selected;
    std::vector<float> vals_isoRaw, vals_isoRel, vals_isoRelUncorrPt; // using regressed pt vs not
    
    // for debugging purposes, add everything after each type of veto to the sum
    std::vector<float> vals_isoRawSumAll, vals_isoRawSelfVetoOnly;
    std::vector<float> vals_isoRelSumAll, vals_isoRelSelfVetoOnly;
    std::vector<float> vals_isoRelSumAllUncorrPt, vals_isoRelSelfVetoOnlyUncorrPt;
    std::vector<float> vals_nPfAll, vals_nPfDr0p3, vals_nPfSelfVetoOnly, vals_nPfDz;

    // Get PF candidates and TkElectron candidates
    for (auto & pf_cands : pf_cands_) {
        // get and select
        iEvent.getByToken(pf_cands.src, src);
        for (const auto & k : *src) {
            pf_selected.push_back(&k);
        }
    }

    for (auto & tkele_cands : tkele_cands_) {
        // get and select
        iEvent.getByToken(tkele_cands.src, src);
        for (const auto & j : *src) {
            tkele_selected.push_back(&j);
        }
    }

    // create the table
    unsigned int ncands = tkele_selected.size();
    unsigned int ncands_pf = pf_selected.size();
    auto out = std::make_unique<nanoaod::FlatTable>(ncands, "TkEleL2", false, true);

    // resize the vectors per electron candidate size
    vals_isoRaw.resize(ncands);
    vals_isoRel.resize(ncands);
    vals_isoRelUncorrPt.resize(ncands);
 
    vals_isoRawSumAll.resize(ncands);
    vals_isoRawSelfVetoOnly.resize(ncands);
    vals_isoRelSumAll.resize(ncands);
    vals_isoRelSelfVetoOnly.resize(ncands);
    vals_isoRelSumAllUncorrPt.resize(ncands);
    vals_isoRelSelfVetoOnlyUncorrPt.resize(ncands);

    vals_nPfAll.resize(ncands);
    vals_nPfDr0p3.resize(ncands);
    vals_nPfSelfVetoOnly.resize(ncands);
    vals_nPfDz.resize(ncands);

    const float bz = 3.8112; // for caloeta/phi calculation
    
    // loop over electrons
    for (unsigned int iEle = 0; iEle < ncands; ++iEle) {
        float isoRawSumAll = 0.;
        float isoRawSelfVetoOnly = 0.;
        float isoRaw = 0.; // after dz veto

        float nPfAll = ncands_pf;
        float nPfDr0p3 = 0;
        float nPfSelfVetoOnly = 0;
        float nPfDz = 0;

	// recast to TkElectron object to access some variables
        const auto * tkEle = dynamic_cast<const l1t::TkElectron*>(tkele_selected[iEle]);
         
        math::XYZTLorentzVector vertex(tkEle->vx(),tkEle->vy(),tkEle->vz(),0.);
        auto caloetaphi = l1tpf::propagateToCalo(tkEle->p4(),vertex,tkEle->charge(),bz);
        float caloeta = caloetaphi.first;
        float calophi = caloetaphi.second;

        //std::cout << "Electron eta = " << tkEle->eta() << ", caloeta = " << caloeta << ", phi = " << tkEle->phi() << ", calophi = " << calophi << ", default caloeta = " << tkEle->egCaloPtr()->eta() <<  std::endl; // FIXME: the default caloeta vs. manually propagated caloeta are different
       
	// for each electron, get the nearby PFs by checking deltaR(ele, PF) < 0.3
        for (unsigned int iPF = 0; iPF < ncands_pf; ++iPF) {
            // use electron caloeta/phi for NEUTRAL pf candidates (use the default caloeta/phi stored in tkEle->egCaloPtr()
            float eta = (pf_selected[iPF]->charge() != 0) ? tkEle->eta() : tkEle->egCaloPtr()->eta();
            float phi = (pf_selected[iPF]->charge() != 0) ? tkEle->phi() : tkEle->egCaloPtr()->phi();

            float dR_ele_pf = reco::deltaR(pf_selected[iPF]->eta(), pf_selected[iPF]->phi(), eta, phi); 
        
            if (dR_ele_pf > 0.3) continue;

            nPfDr0p3 += 1;
            isoRawSumAll += pf_selected[iPF]->pt();

            // self-veto; if deltaR(ele, PF) < 0.05, then do not add to the sum
            if (dR_ele_pf < 0.05) continue;

            nPfSelfVetoOnly += 1;
            isoRawSelfVetoOnly += pf_selected[iPF]->pt();

            // same vertex requirement; for charged PF, add to the sum only if delta vz (ele, charged PF) < 0.5
            if (pf_selected[iPF]->charge() != 0) {
                float dz = std::abs(tkEle->trkzVtx() - pf_selected[iPF]->vz()); // tkEle->vz() gives zero
                if (dz > 0.5) continue;
            }

            nPfDz += 1;
            isoRaw += pf_selected[iPF]->pt();
        }
            
        vals_nPfAll[iEle] = nPfAll; // should be the same for all ele in the event
        vals_nPfDr0p3[iEle] = nPfDr0p3;
        vals_nPfSelfVetoOnly[iEle] = nPfSelfVetoOnly;
        vals_nPfDz[iEle] = nPfDz;

        vals_isoRawSumAll[iEle] = isoRawSumAll; 
        vals_isoRawSelfVetoOnly[iEle] = isoRawSelfVetoOnly; 
        vals_isoRaw[iEle] = isoRaw; 

        float ptCorr = tkEle->userFloat("ptCorr");
	vals_isoRelSumAll[iEle] = isoRawSumAll / ptCorr;
        vals_isoRelSelfVetoOnly[iEle] = isoRawSelfVetoOnly / ptCorr;
        vals_isoRel[iEle] = isoRaw / ptCorr;
   
        float pt = tkEle->pt();	
        vals_isoRelSumAllUncorrPt[iEle] = isoRawSumAll / pt;
        vals_isoRelSelfVetoOnlyUncorrPt[iEle] = isoRawSelfVetoOnly / pt;
        vals_isoRelUncorrPt[iEle] = isoRaw / pt;
    
        //std::cout << "nPFDz = " << nPfDz << ", raw isolation after dz = " << isoRaw << ", relative isolation with uncorrected pt = " << isoRaw/pt << ", w/ regressed pt = " << isoRaw/ptCorr << std::endl;
    }
  
    out->addColumn<float>("nPfAll", vals_nPfAll, "number of PF candidates in the event");
    out->addColumn<float>("nPfDr0p3", vals_nPfDr0p3, "number of PF candidates within dR < 0.3");
    out->addColumn<float>("nPfSelfVetoOnly", vals_nPfSelfVetoOnly, "number of PF candidates within dR < 0.3 & self veto");
    out->addColumn<float>("nPfDz", vals_nPfDz, "number of PF candidates within dR < 0.3 & self veto & dz");
    
    out->addColumn<float>("customPfIsoRawSumAll", vals_isoRawSumAll, "custom PF iso (no veto)");
    out->addColumn<float>("customPfIsoRawSelfVetoOnly", vals_isoRawSelfVetoOnly, "custom PF iso (self veto only)");
    out->addColumn<float>("customPfIsoRaw", vals_isoRaw, "custom PF iso");
    
    out->addColumn<float>("customPfIsoRelSumAll", vals_isoRelSumAll, "custom PF iso relative (no veto)");
    out->addColumn<float>("customPfIsoRelSelfVetoOnly", vals_isoRelSelfVetoOnly, "custom PF iso relative (self veto only)");
    out->addColumn<float>("customPfIsoRel", vals_isoRel, "custom PF iso relative");

    out->addColumn<float>("customPfIsoRelSumAllUncorrPt", vals_isoRelSumAllUncorrPt, "custom PF iso relative (no veto); using uncorrected pt");
    out->addColumn<float>("customPfIsoRelSelfVetoOnlyUncorrPt", vals_isoRelSelfVetoOnlyUncorrPt, "custom PF iso relative (self veto only); using uncorrected pt");
    out->addColumn<float>("customPfIsoRelUncorrPt", vals_isoRelUncorrPt, "custom PF iso relative; using uncorrected pt");

    // save to the event branches
    iEvent.put(std::move(out));

    // clear
    tkele_selected.clear();
    pf_selected.clear();
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1TkEleL2IsoTableProducer);
