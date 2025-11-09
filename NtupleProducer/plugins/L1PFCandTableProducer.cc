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

struct PtSums {
    double genPtSum = 0, genNeutralPtSum = 0, genChargedPtSum = 0, genChargedPtHadSum = 0, genNeutralPtHadSum = 0;
    double recoPtSum = 0, recoNeutralPtSum = 0, recoChargedPtSum = 0, recoChargedPtHadSum = 0, recoNeutralPtHadSum = 0;
    double genRecoRatio = 0;
    int isGenMatched = 0;
};

// Helper function to compute pT sums in a cone
PtSums computePtSumsForCone(
    const reco::Candidate* cand,
    const std::vector<const reco::Candidate*>& selected,
    const std::vector<const reco::Candidate*>& gen_selected,
    double coneSize,
    double bz
) {
    PtSums sums;

    math::XYZTLorentzVector vertex(cand->vx(), cand->vy(), cand->vz(), 0.);
    auto caloetaphi = l1tpf::propagateToCalo(cand->p4(), vertex, cand->charge(), bz);
    double eta1 = (cand->charge() == 0) ? caloetaphi.first : cand->eta();
    double phi1 = (cand->charge() == 0) ? caloetaphi.second : cand->phi();

    // --- RECO loop ---
    for (unsigned int j = 0; j < selected.size(); ++j) {
        const auto* other = selected[j];
        if (cand == other) continue;

        math::XYZTLorentzVector vertex2(other->vx(), other->vy(), other->vz(), 0.);
        auto caloetaphi2 = l1tpf::propagateToCalo(other->p4(), vertex2, other->charge(), bz);
        double eta2 = (other->charge() == 0) ? caloetaphi2.first : other->eta();
        double phi2 = (other->charge() == 0) ? caloetaphi2.second : other->phi();

        double deltaR = calculate_deltaR(eta1, phi1, eta2, phi2);
        if (deltaR >= coneSize) continue;

        sums.recoPtSum += other->pt();
        if (other->charge() == 0) sums.recoNeutralPtSum += other->pt();
        else sums.recoChargedPtSum += other->pt();
        if (abs(other->pdgId()) == 211) sums.recoChargedPtHadSum += other->pt();
        if (abs(other->pdgId()) == 130) sums.recoNeutralPtHadSum += other->pt();
    }

    // include self
    sums.recoPtSum += cand->pt();
    if (cand->charge() == 0) sums.recoNeutralPtSum += cand->pt();
    else sums.recoChargedPtSum += cand->pt();
    if (abs(cand->pdgId()) == 211) sums.recoChargedPtHadSum += cand->pt();
    if (abs(cand->pdgId()) == 130) sums.recoNeutralPtHadSum += cand->pt();

    // --- GEN loop ---
    double min_dR = 999.;
    int idx_min_dR = -1;

    for (unsigned int k = 0; k < gen_selected.size(); ++k) {
        const auto* gen = gen_selected[k];

        math::XYZTLorentzVector vertex3(gen->vx(), gen->vy(), gen->vz(), 0.);
        auto caloetaphi3 = l1tpf::propagateToCalo(gen->p4(), vertex3, gen->charge(), bz);
        double eta3 = (gen->charge() == 0) ? caloetaphi3.first : gen->eta();
        double phi3 = (gen->charge() == 0) ? caloetaphi3.second : gen->phi();

        double deltaR = calculate_deltaR(eta1, phi1, eta3, phi3);
        if (deltaR < min_dR) {
            min_dR = deltaR;
            idx_min_dR = k;
        }

        if (deltaR < coneSize) {
            sums.genPtSum += gen->pt();
            if (gen->charge() == 0) sums.genNeutralPtSum += gen->pt();
            else sums.genChargedPtSum += gen->pt();
            if (abs(gen->pdgId()) == 211) sums.genChargedPtHadSum += gen->pt();
            if (abs(gen->pdgId()) == 130) sums.genNeutralPtHadSum += gen->pt();
        }
    }

    // --- Explicitly add closest gen match for charged sum
    if (idx_min_dR >= 0 && min_dR < coneSize) {
        const auto* bestMatch = gen_selected[idx_min_dR];
        sums.genChargedPtSum     += bestMatch->pt();
        sums.genChargedPtHadSum  += bestMatch->pt();
    }

    // --- Gen match flag
    if (min_dR < 0.1) sums.isGenMatched = 1;

    // --- Ratio
    sums.genRecoRatio = (sums.recoPtSum > 0) ? sums.genPtSum / sums.recoPtSum : 0.0;

    return sums;
}


void
L1PFCandTableProducer::produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const
{
    edm::Handle<reco::CandidateView> src;
    std::vector<const reco::Candidate *> selected;
    std::vector<const reco::Candidate *> gen_selected;
    std::vector<float> vals_pt, vals_eta, vals_phi, vals_mass;

    for (auto & gencands : gencands_) {
        iEvent.getByToken(gencands.src, src);
        for (const auto & j : *src) {
            if (sel_(j) && gencands.sel(j)) {
                gen_selected.push_back(&j);
            }
        }
    }

    for (auto & cands : cands_) {
        iEvent.getByToken(cands.src, src);
        for (const auto & j : *src) {
            if (sel_(j) && cands.sel(j)) {
                selected.push_back(&j);
            }
        }

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

        const float bz = 3.8112;

        // Declare vectors for both cones
        std::vector<float> vals_caloeta(ncands), vals_calophi(ncands);
        std::vector<int> vals_isGenMatched(ncands);

        std::vector<float> vals_genPtSum0p2(ncands), vals_genNeutralPtSum0p2(ncands),
                           vals_genChargedPtSum0p2(ncands), vals_genChargedPtHadSum0p2(ncands),
                           vals_genNeutralPtHadSum0p2(ncands);

        std::vector<float> vals_recoPtSum0p2(ncands), vals_recoNeutralPtSum0p2(ncands),
                           vals_recoChargedPtSum0p2(ncands), vals_recoChargedPtHadSum0p2(ncands),
                           vals_recoNeutralPtHadSum0p2(ncands);

        std::vector<float> vals_genRecoPtRatio0p2(ncands);

        std::vector<float> vals_genPtSum0p3(ncands), vals_genNeutralPtSum0p3(ncands),
                           vals_genChargedPtSum0p3(ncands), vals_genChargedPtHadSum0p3(ncands),
                           vals_genNeutralPtHadSum0p3(ncands);

        std::vector<float> vals_recoPtSum0p3(ncands), vals_recoNeutralPtSum0p3(ncands),
                           vals_recoChargedPtSum0p3(ncands), vals_recoChargedPtHadSum0p3(ncands),
                           vals_recoNeutralPtHadSum0p3(ncands);

        std::vector<float> vals_genRecoPtRatio0p3(ncands);

        // --- main loop per candidate
        for (unsigned int i = 0; i < ncands; ++i) {
            const auto* cand = selected[i];

            math::XYZTLorentzVector vertex(cand->vx(), cand->vy(), cand->vz(), 0.);
            auto caloetaphi = l1tpf::propagateToCalo(cand->p4(), vertex, cand->charge(), bz);
            vals_caloeta[i] = caloetaphi.first;
            vals_calophi[i] = caloetaphi.second;

            // compute sums for both cone sizes
            auto sums0p2 = computePtSumsForCone(cand, selected, gen_selected, 0.2, bz);
            auto sums0p3 = computePtSumsForCone(cand, selected, gen_selected, 0.3, bz);

            vals_genPtSum0p2[i] = sums0p2.genPtSum;
            vals_genNeutralPtSum0p2[i] = sums0p2.genNeutralPtSum;
            vals_genChargedPtSum0p2[i] = sums0p2.genChargedPtSum;
            vals_genChargedPtHadSum0p2[i] = sums0p2.genChargedPtHadSum;
            vals_genNeutralPtHadSum0p2[i] = sums0p2.genNeutralPtHadSum;

            vals_recoPtSum0p2[i] = sums0p2.recoPtSum;
            vals_recoNeutralPtSum0p2[i] = sums0p2.recoNeutralPtSum;
            vals_recoChargedPtSum0p2[i] = sums0p2.recoChargedPtSum;
            vals_recoChargedPtHadSum0p2[i] = sums0p2.recoChargedPtHadSum;
            vals_recoNeutralPtHadSum0p2[i] = sums0p2.recoNeutralPtHadSum;

            vals_genRecoPtRatio0p2[i] = sums0p2.genRecoRatio;
            vals_isGenMatched[i] = sums0p2.isGenMatched;

            vals_genPtSum0p3[i] = sums0p3.genPtSum;
            vals_genNeutralPtSum0p3[i] = sums0p3.genNeutralPtSum;
            vals_genChargedPtSum0p3[i] = sums0p3.genChargedPtSum;
            vals_genChargedPtHadSum0p3[i] = sums0p3.genChargedPtHadSum;
            vals_genNeutralPtHadSum0p3[i] = sums0p3.genNeutralPtHadSum;

            vals_recoPtSum0p3[i] = sums0p3.recoPtSum;
            vals_recoNeutralPtSum0p3[i] = sums0p3.recoNeutralPtSum;
            vals_recoChargedPtSum0p3[i] = sums0p3.recoChargedPtSum;
            vals_recoChargedPtHadSum0p3[i] = sums0p3.recoChargedPtHadSum;
            vals_recoNeutralPtHadSum0p3[i] = sums0p3.recoNeutralPtHadSum;

            vals_genRecoPtRatio0p3[i] = sums0p3.genRecoRatio;
        }

        // --- Add columns for 0p2 cone
        out->addColumn<float>("genPtSum0p2", vals_genPtSum0p2, "");
        out->addColumn<float>("genNeutralPtSum0p2", vals_genNeutralPtSum0p2, "");
        out->addColumn<float>("genChargedPtSum0p2", vals_genChargedPtSum0p2, "");
        out->addColumn<float>("genChargedHadPtSum0p2", vals_genChargedPtHadSum0p2, "");
        out->addColumn<float>("genNeutralHadPtSum0p2", vals_genNeutralPtHadSum0p2, "");

        out->addColumn<float>("recoPtSum0p2", vals_recoPtSum0p2, "");
        out->addColumn<float>("recoNeutralPtSum0p2", vals_recoNeutralPtSum0p2, "");
        out->addColumn<float>("recoChargedPtSum0p2", vals_recoChargedPtSum0p2, "");
        out->addColumn<float>("recoChargedHadPtSum0p2", vals_recoChargedPtHadSum0p2, "");
        out->addColumn<float>("recoNeutralHadPtSum0p2", vals_recoNeutralPtHadSum0p2, "");

        out->addColumn<float>("genRecoRatio0p2", vals_genRecoPtRatio0p2, "");

        // --- Add columns for 0p3 cone
        out->addColumn<float>("genPtSum0p3", vals_genPtSum0p3, "");
        out->addColumn<float>("genNeutralPtSum0p3", vals_genNeutralPtSum0p3, "");
        out->addColumn<float>("genChargedPtSum0p3", vals_genChargedPtSum0p3, "");
        out->addColumn<float>("genChargedHadPtSum0p3", vals_genChargedPtHadSum0p3, "");
        out->addColumn<float>("genNeutralHadPtSum0p3", vals_genNeutralPtHadSum0p3, "");

        out->addColumn<float>("recoPtSum0p3", vals_recoPtSum0p3, "");
        out->addColumn<float>("recoNeutralPtSum0p3", vals_recoNeutralPtSum0p3, "");
        out->addColumn<float>("recoChargedPtSum0p3", vals_recoChargedPtSum0p3, "");
        out->addColumn<float>("recoChargedHadPtSum0p3", vals_recoChargedPtHadSum0p3, "");
        out->addColumn<float>("recoNeutralHadPtSum0p3", vals_recoNeutralPtHadSum0p3, "");

        out->addColumn<float>("genRecoRatio0p3", vals_genRecoPtRatio0p3, "");

        out->addColumn<int>("isGenMatched", vals_isGenMatched, "");

        out->addColumn<float>("caloeta", vals_caloeta, "");
        out->addColumn<float>("calophi", vals_calophi, "");

        // save to the event branches
        iEvent.put(std::move(out), cands.coll+"Cands");

        selected.clear();
    }
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1PFCandTableProducer);
