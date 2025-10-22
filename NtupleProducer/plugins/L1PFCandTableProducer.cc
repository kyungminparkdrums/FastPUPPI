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

#include <algorithm>

#include "L1Trigger/Phase2L1ParticleFlow/interface/L1TPFUtils.h"

class L1PFCandTableProducer : public edm::global::EDProducer<>  {
    public:
        explicit L1PFCandTableProducer(const edm::ParameterSet&);
        ~L1PFCandTableProducer();

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
        std::vector<CandRecord> cands_;

	std::vector<CandRecord> gencands_;
};

L1PFCandTableProducer::L1PFCandTableProducer(const edm::ParameterSet& iConfig) :
    sel_(iConfig.getParameter<std::string>("commonSel"), true)
{
    edm::ParameterSet cands = iConfig.getParameter<edm::ParameterSet>("cands");
    auto candnames = cands.getParameterNamesForType<edm::InputTag>();
    for (const std::string & name : candnames) {
        cands_.emplace_back(name, consumes<reco::CandidateView>(cands.getParameter<edm::InputTag>(name)), cands);
        produces<nanoaod::FlatTable>(name+"Cands");
    }

    gencands_.emplace_back("GenCands", consumes<reco::CandidateView>(edm::InputTag("genParticlesForMETAllVisible")), iConfig);

    if (iConfig.existsAs<edm::ParameterSet>("moreVariables")) {
        edm::ParameterSet vars = iConfig.getParameter<edm::ParameterSet>("moreVariables");
        auto morenames = vars.getParameterNamesForType<std::string>();
        for (const std::string & name : morenames) {
            extraVars_.emplace_back(name, vars.getParameter<std::string>(name));
        }
    }
 }

double calculate_deltaR(double eta1, double phi1, double eta2, double phi2) {
    // Build pseudo-vectors (PtEtaPhiM) with dummy pT/M
    reco::Candidate::LorentzVector v1(1.0, eta1, phi1, 0.0);
    reco::Candidate::LorentzVector v2(1.0, eta2, phi2, 0.0);

    return reco::deltaR(v1, v2);
}

L1PFCandTableProducer::~L1PFCandTableProducer() { }

// ------------ method called for each event  ------------
    void
L1PFCandTableProducer::produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const
{
    edm::Handle<reco::CandidateView> src;
    std::vector<const reco::Candidate *> selected;
    std::vector<const reco::Candidate *> gen_selected;
    std::vector<float> vals_pt, vals_eta, vals_phi, vals_mass;

    for (auto & gencands : gencands_) {
        // get and select
        iEvent.getByToken(gencands.src, src);
        for (const auto & j : *src) {
            if (sel_(j) && gencands.sel(j)) {
                gen_selected.push_back(&j);
            }
        }
    } 

    for (auto & cands : cands_) {
        // get and select
        iEvent.getByToken(cands.src, src);
        for (const auto & j : *src) {
            if (sel_(j) && cands.sel(j)) {
                selected.push_back(&j);
            }
        }
        
        // create the table
        unsigned int nGenCands = gen_selected.size();
        unsigned int ncands = selected.size();
        auto out = std::make_unique<nanoaod::FlatTable>(ncands, cands.coll+"Cands", false);

        // fill basic info
        vals_pt.resize(ncands); 
        vals_eta.resize(ncands); 
        vals_phi.resize(ncands); 
        vals_mass.resize(ncands); 
        for (unsigned int i = 0; i < ncands; ++i) {
            vals_pt[i] = selected[i]->pt();
            vals_eta[i] = selected[i]->eta();
            vals_phi[i] = selected[i]->phi();
            vals_mass[i] = selected[i]->mass();
        }
        out->addColumn<float>("pt", vals_pt, "pt of cand");
        out->addColumn<float>("eta", vals_eta, "eta of cand");
        out->addColumn<float>("phi", vals_phi, "phi of cand");
        out->addColumn<float>("mass", vals_mass, "mass of cand");

        // fill extra vars
        for (const auto & evar : extraVars_) {
            for (unsigned int i = 0; i < ncands; ++i) {
                vals_pt[i] = evar.func(*selected[i]);
            }
            out->addColumn<float>(evar.name, vals_pt, evar.expr);
        }

        // Add caloeta, calophi
	const float bz = 3.8112;

        std::vector<float> vals_caloeta, vals_calophi;
        std::vector<float> vals_genPtSum, vals_genNeutralPtSum, vals_recoPtSum, vals_genRecoPtRatio;
        vals_caloeta.resize(ncands);
        vals_calophi.resize(ncands);
        
	vals_genNeutralPtSum.resize(ncands);
	vals_genPtSum.resize(ncands);
	vals_recoPtSum.resize(ncands);
	vals_genRecoPtRatio.resize(ncands);

        for (unsigned int i = 0; i < ncands; ++i) {
            math::XYZTLorentzVector vertex(selected[i]->vx(),selected[i]->vy(),selected[i]->vz(),0.);
            auto caloetaphi = l1tpf::propagateToCalo(selected[i]->p4(),vertex,selected[i]->charge(),bz);
            vals_caloeta[i] = caloetaphi.first;
            vals_calophi[i] = caloetaphi.second;

            double eta1 = (selected[i]->charge() == 0) ? caloetaphi.first : selected[i]->eta();
            double phi1 = (selected[i]->charge() == 0) ? caloetaphi.second : selected[i]->phi();

            double sum_pT_gen = 0;
            double sum_pT_genNeutral = 0;
            double sum_pT_reco = 0;

	    // pT sum around a cone
            for (unsigned int j = 0; j < ncands; ++j) {
                math::XYZTLorentzVector vertex2(selected[j]->vx(),selected[j]->vy(),selected[j]->vz(),0.);
                auto caloetaphi2 = l1tpf::propagateToCalo(selected[j]->p4(),vertex2,selected[j]->charge(),bz);
                double eta2 = (selected[j]->charge() == 0) ? caloetaphi2.first : selected[j]->eta();
                double phi2 = (selected[j]->charge() == 0) ? caloetaphi2.first : selected[j]->phi();

                double deltaR = calculate_deltaR(eta1, phi1, eta2, phi2);

                if (deltaR < 0.2) sum_pT_reco += selected[j]->pt();
            }

	    // pT sum around a gen cone
            for (unsigned int k = 0; k < nGenCands; ++k) {
                math::XYZTLorentzVector vertex3(gen_selected[k]->vx(),gen_selected[k]->vy(),gen_selected[k]->vz(),0.);
                auto caloetaphi3 = l1tpf::propagateToCalo(gen_selected[k]->p4(),vertex3,gen_selected[k]->charge(),bz);
                double eta3 = (gen_selected[k]->charge() == 0) ? caloetaphi3.first : gen_selected[k]->eta();
                double phi3 = (gen_selected[k]->charge() == 0) ? caloetaphi3.first : gen_selected[k]->phi();

                double deltaR = calculate_deltaR(eta1, phi1, eta3, phi3);

                if (deltaR < 0.2) sum_pT_gen += gen_selected[k]->pt();
                if ((deltaR < 0.2) & (gen_selected[k]->charge() == 0)) sum_pT_genNeutral += gen_selected[k]->pt();
            }

            vals_genNeutralPtSum[i] = sum_pT_genNeutral;
            vals_genPtSum[i] = sum_pT_gen;
            vals_recoPtSum[i] = sum_pT_reco;
            vals_genRecoPtRatio[i] = sum_pT_gen / sum_pT_reco;
        }

        out->addColumn<float>("genPtSum0p2", vals_genPtSum, "");
        out->addColumn<float>("genNeutralPtSum0p2", vals_genNeutralPtSum, "");
        out->addColumn<float>("recoPtSum0p2", vals_recoPtSum, "");
        out->addColumn<float>("genRecoRatio0p2", vals_genRecoPtRatio, "");

        out->addColumn<float>("caloeta", vals_caloeta, "");
        out->addColumn<float>("calophi", vals_calophi, "");

        // save to the event branches
        iEvent.put(std::move(out), cands.coll+"Cands");

        // clear
        selected.clear();
    }
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1PFCandTableProducer);
