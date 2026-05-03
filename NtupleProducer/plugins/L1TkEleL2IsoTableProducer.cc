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
    std::vector<float> vals_isoRaw, vals_isoRawOtherEleVeto, vals_isoRel, vals_isoRelOtherEleVeto;
    std::vector<float> vals_isoRawNoCaloEtaPhi, vals_isoRawOtherEleVetoNoCaloEtaPhi, vals_isoRelNoCaloEtaPhi, vals_isoRelOtherEleVetoNoCaloEtaPhi; // using eta-phi only (without calo eta/phi for neutrals) // variable names are getting longer and longer

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

	// create the table
	unsigned int ncands = tkele_selected.size();
	unsigned int ncands_pf = pf_selected.size();
        auto out = std::make_unique<nanoaod::FlatTable>(ncands, "TkEleL2", false, true);

	// fill in the table
	vals_isoRaw.resize(ncands);
	vals_isoRel.resize(ncands);
	vals_isoRawOtherEleVeto.resize(ncands);
	vals_isoRelOtherEleVeto.resize(ncands);

        const float bz = 3.8112; // for caloeta/phi calculation
	
	// loop over electrons
	for (unsigned int iEle = 0; iEle < ncands; ++iEle) {
	    float isoRaw = 0.;
	    float isoRaw_otherEleVeto = 0.;

	    float isoRawNoCaloEtaPhi = 0.;
	    float isoRaw_otherEleVetoNoCaloEtaPhi = 0.;

	    // for each electron, get the nearby PFs by checking deltaR(ele, PF) < 0.3
	    for (unsigned int iPF = 0; iPF < ncands_pf; ++iPF) {
		// use caloeta/phi for NEUTRAL pf candidates
                math::XYZTLorentzVector vertex(pf_selected[iPF]->vx(),pf_selected[iPF]->vy(),pf_selected[iPF]->vz(),0.);
		auto caloetaphi = l1tpf::propagateToCalo(pf_selected[iPF]->p4(),vertex,pf_selected[iPF]->charge(),bz);
                float caloeta = caloetaphi.first;
                float calophi = caloetaphi.second;
   
                float eta = (pf_selected[iPF]->charge() != 0) ? pf_selected[iPF]->eta() : caloeta;
                float phi = (pf_selected[iPF]->charge() != 0) ? pf_selected[iPF]->phi() : calophi;

		float dR_ele_pf = reco::deltaR(tkele_selected[iEle]->eta(), tkele_selected[iEle]->phi(), eta, phi);	
                
		if (dR_ele_pf > 0.3) continue;

		// self-veto; if deltaR(ele, PF) < 0.05, then do not add to the sum
		if (dR_ele_pf < 0.05) continue;

		// same vertex requirement; for charged PF, add to the sum only if delta vz (ele, charged PF) < 0.5
                if (pf_selected[iPF]->charge() != 0) {
                    float dz = std::abs(tkele_selected[iEle]->vz() - pf_selected[iPF]->vz());
		    if (dz > 0.5) continue;
		}

		isoRaw += pf_selected[iPF]->pt();

		bool veto_by_other_ele = false;
		// additionally, other electron veto; if deltaR(another ele, PF) < 0.02, then do not add to the sum
	        for (unsigned int jEle = 0; jEle < ncands; ++jEle) {
                    if (jEle == iEle) continue;

		    float dR_other_pf = reco::deltaR(tkele_selected[jEle]->eta(), tkele_selected[jEle]->phi(), eta, phi);
		    if (dR_other_pf < 0.02) {
                        veto_by_other_ele = true;
			break;
	            }
	        }
		if (veto_by_other_ele) continue;

	        isoRaw_otherEleVeto += pf_selected[iPF]->pt();
	    }

	    vals_isoRaw[iEle] = isoRaw; 
	    vals_isoRawOtherEleVeto[iEle] = isoRaw_otherEleVeto; 
	
	    const auto * tkEle = dynamic_cast<const l1t::TkElectron*>(tkele_selected[iEle]);
	    float ptCorr = tkEle->userFloat("ptCorr");

	    vals_isoRel[iEle] = isoRaw / ptCorr;
	    vals_isoRelOtherEleVeto[iEle] = isoRaw_otherEleVeto / ptCorr;
	
	}
        out->addColumn<float>("customPfIsoRaw", vals_isoRaw, "custom PF iso");
        out->addColumn<float>("customPfIsoRawOtherEleVeto", vals_isoRawOtherEleVeto, "custom PF iso (w/ other electron veto)");
        out->addColumn<float>("customPfIsoRel", vals_isoRel, "custom PF iso relative");
        out->addColumn<float>("customPfIsoRelOtherEleVeto", vals_isoRelOtherEleVeto, "custom PF iso relative (w/ other electron veto)");

	// save to the event branches
        iEvent.put(std::move(out));

        // clear
        tkele_selected.clear();
        pf_selected.clear();
    } 
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1TkEleL2IsoTableProducer);
