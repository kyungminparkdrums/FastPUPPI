// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"

#include "L1Trigger/Phase2L1ParticleFlow/interface/L1TPFUtils.h"

#include <algorithm>

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
    return reco::deltaR(eta1, phi1, eta2, phi2);
}

L1PFCandTableProducer::~L1PFCandTableProducer() { }

struct GenCounts {
    int nGenInCone = 0;
    int nGenStatus1InCone = 0;
    int nGenPtThrInCone = 0;
    int nGenStatus1PtThrInCone = 0;

    int nChargedGenInCone = 0;
    int nChargedGenStatus1InCone = 0;
    int nChargedGenPt2InCone = 0;
    int nChargedGenStatus1Pt2InCone = 0;

    int nNeutralGenInCone = 0;
    int nNeutralGenStatus1InCone = 0;
    int nNeutralGenPt1InCone = 0;
    int nNeutralGenStatus1Pt1InCone = 0;

    int nChargedHadGenInCone = 0;
    int nChargedHadGenStatus1InCone = 0;
    int nChargedHadGenPt2InCone = 0;
    int nChargedHadGenStatus1Pt2InCone = 0;

    int nNeutralHadGenInCone = 0;
    int nNeutralHadGenStatus1InCone = 0;
    int nNeutralHadGenPt1InCone = 0;
    int nNeutralHadGenStatus1Pt1InCone = 0;
};

struct PtSums {
    double genPtSum = 0, genNeutralPtSum = 0, genChargedPtSum = 0, genChargedPtHadSum = 0, genNeutralPtHadSum = 0;
    double recoPtSum = 0, recoNeutralPtSum = 0, recoChargedPtSum = 0, recoChargedPtHadSum = 0, recoNeutralPtHadSum = 0;
    double genRecoRatio = 0;
    int isGenMatched = 0;
    GenCounts counts;
};

// ----------------------------------------------------------------------
// Mask that defines which GEN-count columns to create for each category
// ----------------------------------------------------------------------
struct GenColumnMask {
    bool makePtThr = false;  // GEN-only
    bool makePt1   = false;  // neutrals
    bool makePt2   = false;  // charged
};


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

    // ---- RECO loop ----
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

    // ---- GEN loop (status==1 only) ----
    double min_dR = 999.;
    int idx_min_dR = -1;

    for (unsigned int k = 0; k < gen_selected.size(); ++k) {
        const auto* gen = gen_selected[k];
        const reco::GenParticle* gp = dynamic_cast<const reco::GenParticle*>(gen);
        if (!gp || gp->status() != 1) continue; // only stable

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

	  // --- classification ---
	  bool isNeutral = (gen->charge() == 0);
	  bool passNeutral = (isNeutral && gen->pt() > 1);
	  bool passCharged = (!isNeutral && gen->pt() > 2);
	  bool isChargedHad = (abs(gp->pdgId()) == 211);
	  bool isNeutralHad = (abs(gp->pdgId()) == 130);
	  bool passNeutralHad = (isNeutralHad && gen->pt() > 1);

	  // --- generic counts ---
	  sums.counts.nGenInCone++;
	  sums.counts.nGenStatus1InCone++;

	  if (passNeutral || passCharged) {
	    sums.counts.nGenPtThrInCone++;
	    sums.counts.nGenStatus1PtThrInCone++; 
	  }
    	  
	  // --- charged vs neutral ---
	  if (!isNeutral) {
	    sums.counts.nChargedGenInCone++;
	    sums.counts.nChargedGenStatus1InCone++;
	    
	    if (passCharged) {
	      sums.counts.nChargedGenPt2InCone++;
	      sums.counts.nChargedGenStatus1Pt2InCone++;
	    }
	  } else {
	    sums.counts.nNeutralGenInCone++;
	    sums.counts.nNeutralGenStatus1InCone++;

	    if (passNeutral) {
	      sums.counts.nNeutralGenPt1InCone++;
	      sums.counts.nNeutralGenStatus1Pt1InCone++;
	    }
	  }
	  
	  // --- hadron IDs ---
	  if (isChargedHad) {
	    sums.counts.nChargedHadGenInCone++;
	    sums.counts.nChargedHadGenStatus1InCone++;
	    
	    if (passCharged) {
	      sums.counts.nChargedHadGenPt2InCone++;
	      sums.counts.nChargedHadGenStatus1Pt2InCone++;
	    }
	  }

	  if (isNeutralHad) {
	    sums.counts.nNeutralHadGenInCone++;
	    sums.counts.nNeutralHadGenStatus1InCone++;
	    
	    if (passNeutralHad) {
	      sums.counts.nNeutralHadGenPt1InCone++;
	      sums.counts.nNeutralHadGenStatus1Pt1InCone++;
	    }
	  }

	    
	  // --- pT sums (status==1 only) ---
	  if (passNeutral || passCharged)
	    sums.genPtSum += gen->pt();
	  
	  if (passNeutral)
	    sums.genNeutralPtSum += gen->pt();
	  
	  if (passCharged)
	    sums.genChargedPtSum += gen->pt();
	  
	  if (abs(gen->pdgId()) == 211 && passCharged)
	    sums.genChargedPtHadSum += gen->pt();
	  
	  if (abs(gen->pdgId()) == 130 && passNeutral)
	    sums.genNeutralPtHadSum += gen->pt();
	  
        }
    }

    // closest match among stable particles
    if (idx_min_dR >= 0 && min_dR < coneSize) {
        const auto* bestMatch = gen_selected[idx_min_dR];
        sums.genChargedPtSum     += bestMatch->pt();
        sums.genChargedPtHadSum  += bestMatch->pt();
    }

    if (min_dR < 0.1) sums.isGenMatched = 1;
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

    // ---- collect GEN candidates ----
    for (auto & gencands : gencands_) {
        iEvent.getByToken(gencands.src, src);
        for (const auto & j : *src) {
            if (sel_(j) && gencands.sel(j)) gen_selected.push_back(&j);
        }
    }

    // ---- loop over candidate collections ----
    for (auto & cands : cands_) {
        iEvent.getByToken(cands.src, src);
        for (const auto & j : *src) {
            if (sel_(j) && cands.sel(j)) selected.push_back(&j);
        }

        unsigned int ncands = selected.size();
        auto out = std::make_unique<nanoaod::FlatTable>(ncands, cands.coll+"Cands", false);

        // ---- fill basic info ----
        vals_pt.resize(ncands);
        vals_eta.resize(ncands);
        vals_phi.resize(ncands);
        vals_mass.resize(ncands);
        for (unsigned int i = 0; i < ncands; ++i) {
            vals_pt[i]   = selected[i]->pt();
            vals_eta[i]  = selected[i]->eta();
            vals_phi[i]  = selected[i]->phi();
            vals_mass[i] = selected[i]->mass();
        }
        out->addColumn<float>("pt",   vals_pt,   "pt of cand");
        out->addColumn<float>("eta",  vals_eta,  "eta of cand");
        out->addColumn<float>("phi",  vals_phi,  "phi of cand");
        out->addColumn<float>("mass", vals_mass, "mass of cand");

        // ---- extra user-defined variables ----
        for (const auto & evar : extraVars_) {
            for (unsigned int i = 0; i < ncands; ++i) vals_pt[i] = evar.func(*selected[i]);
            out->addColumn<float>(evar.name, vals_pt, evar.expr);
        }

        const float bz = 3.8112;

        // ---- allocate output vectors ----
        std::vector<int> vals_isGenMatched(ncands);
        std::vector<float> vals_genRecoPtRatio0p2(ncands), vals_genRecoPtRatio0p3(ncands);

        std::vector<float> vals_genPtSum0p2(ncands), vals_genNeutralPtSum0p2(ncands),
                           vals_genChargedPtSum0p2(ncands), vals_genChargedPtHadSum0p2(ncands),
                           vals_genNeutralPtHadSum0p2(ncands);
        std::vector<float> vals_recoPtSum0p2(ncands), vals_recoNeutralPtSum0p2(ncands),
                           vals_recoChargedPtSum0p2(ncands), vals_recoChargedPtHadSum0p2(ncands),
                           vals_recoNeutralPtHadSum0p2(ncands);

        std::vector<float> vals_genPtSum0p3(ncands), vals_genNeutralPtSum0p3(ncands),
                           vals_genChargedPtSum0p3(ncands), vals_genChargedPtHadSum0p3(ncands),
                           vals_genNeutralPtHadSum0p3(ncands);
        std::vector<float> vals_recoPtSum0p3(ncands), vals_recoNeutralPtSum0p3(ncands),
                           vals_recoChargedPtSum0p3(ncands), vals_recoChargedPtHadSum0p3(ncands),
                           vals_recoNeutralPtHadSum0p3(ncands);

        std::vector<float> vals_caloeta(ncands), vals_calophi(ncands);

        // ---- prepare vectors for gen counters ----
        std::vector<int> nGenInCone0p1(ncands,0), nGenStatus1InCone0p1(ncands,0),
                         nGenPtThrInCone0p1(ncands,0), nGenStatus1PtThrInCone0p1(ncands,0);
        std::vector<int> nChargedGenInCone0p1(ncands,0), nChargedGenStatus1InCone0p1(ncands,0),
                         nChargedGenPt2InCone0p1(ncands,0), nChargedGenStatus1Pt2InCone0p1(ncands,0);
        std::vector<int> nNeutralGenInCone0p1(ncands,0), nNeutralGenStatus1InCone0p1(ncands,0),
                         nNeutralGenPt1InCone0p1(ncands,0), nNeutralGenStatus1Pt1InCone0p1(ncands,0);
        std::vector<int> nChargedHadGenInCone0p1(ncands,0), nChargedHadGenStatus1InCone0p1(ncands,0),
                         nChargedHadGenPt2InCone0p1(ncands,0), nChargedHadGenStatus1Pt2InCone0p1(ncands,0);
        std::vector<int> nNeutralHadGenInCone0p1(ncands,0), nNeutralHadGenStatus1InCone0p1(ncands,0),
                         nNeutralHadGenPt1InCone0p1(ncands,0), nNeutralHadGenStatus1Pt1InCone0p1(ncands,0);

        // replicate for 0.2 and 0.3
        std::vector<int> nGenInCone0p2(ncands,0), nGenStatus1InCone0p2(ncands,0),
                         nGenPtThrInCone0p2(ncands,0), nGenStatus1PtThrInCone0p2(ncands,0);
        std::vector<int> nChargedGenInCone0p2(ncands,0), nChargedGenStatus1InCone0p2(ncands,0),
                         nChargedGenPt2InCone0p2(ncands,0), nChargedGenStatus1Pt2InCone0p2(ncands,0);
        std::vector<int> nNeutralGenInCone0p2(ncands,0), nNeutralGenStatus1InCone0p2(ncands,0),
                         nNeutralGenPt1InCone0p2(ncands,0), nNeutralGenStatus1Pt1InCone0p2(ncands,0);
        std::vector<int> nChargedHadGenInCone0p2(ncands,0), nChargedHadGenStatus1InCone0p2(ncands,0),
                         nChargedHadGenPt2InCone0p2(ncands,0), nChargedHadGenStatus1Pt2InCone0p2(ncands,0);
        std::vector<int> nNeutralHadGenInCone0p2(ncands,0), nNeutralHadGenStatus1InCone0p2(ncands,0),
                         nNeutralHadGenPt1InCone0p2(ncands,0), nNeutralHadGenStatus1Pt1InCone0p2(ncands,0);

        std::vector<int> nGenInCone0p3(ncands,0), nGenStatus1InCone0p3(ncands,0),
                         nGenPtThrInCone0p3(ncands,0), nGenStatus1PtThrInCone0p3(ncands,0);
        std::vector<int> nChargedGenInCone0p3(ncands,0), nChargedGenStatus1InCone0p3(ncands,0),
                         nChargedGenPt2InCone0p3(ncands,0), nChargedGenStatus1Pt2InCone0p3(ncands,0);
        std::vector<int> nNeutralGenInCone0p3(ncands,0), nNeutralGenStatus1InCone0p3(ncands,0),
                         nNeutralGenPt1InCone0p3(ncands,0), nNeutralGenStatus1Pt1InCone0p3(ncands,0);
        std::vector<int> nChargedHadGenInCone0p3(ncands,0), nChargedHadGenStatus1InCone0p3(ncands,0),
                         nChargedHadGenPt2InCone0p3(ncands,0), nChargedHadGenStatus1Pt2InCone0p3(ncands,0);
        std::vector<int> nNeutralHadGenInCone0p3(ncands,0), nNeutralHadGenStatus1InCone0p3(ncands,0),
                         nNeutralHadGenPt1InCone0p3(ncands,0), nNeutralHadGenStatus1Pt1InCone0p3(ncands,0);

        // ---- main candidate loop ----
        for (unsigned int i = 0; i < ncands; ++i) {
            const auto* cand = selected[i];

            math::XYZTLorentzVector vertex(cand->vx(), cand->vy(), cand->vz(), 0.);
            auto caloetaphi = l1tpf::propagateToCalo(cand->p4(), vertex, cand->charge(), bz);
            vals_caloeta[i] = caloetaphi.first;
            vals_calophi[i] = caloetaphi.second;

            auto sums0p1 = computePtSumsForCone(cand, selected, gen_selected, 0.1, bz);
            auto sums0p2 = computePtSumsForCone(cand, selected, gen_selected, 0.2, bz);
            auto sums0p3 = computePtSumsForCone(cand, selected, gen_selected, 0.3, bz);

            // --- store main numeric results (0.2/0.3) ---
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

            // --- store gen counters (0.1/0.2/0.3) ---
            nGenInCone0p1[i]                = sums0p1.counts.nGenInCone;
            nGenStatus1InCone0p1[i]         = sums0p1.counts.nGenStatus1InCone;
            nGenPtThrInCone0p1[i]             = sums0p1.counts.nGenPtThrInCone;
            nGenStatus1PtThrInCone0p1[i]      = sums0p1.counts.nGenStatus1PtThrInCone;
            nChargedGenInCone0p1[i]         = sums0p1.counts.nChargedGenInCone;
            nChargedGenStatus1InCone0p1[i]  = sums0p1.counts.nChargedGenStatus1InCone;
            nChargedGenPt2InCone0p1[i]      = sums0p1.counts.nChargedGenPt2InCone;
            nChargedGenStatus1Pt2InCone0p1[i]=sums0p1.counts.nChargedGenStatus1Pt2InCone;
            nNeutralGenInCone0p1[i]         = sums0p1.counts.nNeutralGenInCone;
            nNeutralGenStatus1InCone0p1[i]  = sums0p1.counts.nNeutralGenStatus1InCone;
            nNeutralGenPt1InCone0p1[i]      = sums0p1.counts.nNeutralGenPt1InCone;
            nNeutralGenStatus1Pt1InCone0p1[i]=sums0p1.counts.nNeutralGenStatus1Pt1InCone;
            nChargedHadGenInCone0p1[i]      = sums0p1.counts.nChargedHadGenInCone;
            nChargedHadGenStatus1InCone0p1[i]=sums0p1.counts.nChargedHadGenStatus1InCone;
            nChargedHadGenPt2InCone0p1[i]   = sums0p1.counts.nChargedHadGenPt2InCone;
            nChargedHadGenStatus1Pt2InCone0p1[i]=sums0p1.counts.nChargedHadGenStatus1Pt2InCone;
            nNeutralHadGenInCone0p1[i]      = sums0p1.counts.nNeutralHadGenInCone;
            nNeutralHadGenStatus1InCone0p1[i]=sums0p1.counts.nNeutralHadGenStatus1InCone;
            nNeutralHadGenPt1InCone0p1[i]   = sums0p1.counts.nNeutralHadGenPt1InCone;
            nNeutralHadGenStatus1Pt1InCone0p1[i]=sums0p1.counts.nNeutralHadGenStatus1Pt1InCone;

            // repeat for 0.2 / 0.3
            nGenInCone0p2[i]                = sums0p2.counts.nGenInCone;
            nGenStatus1InCone0p2[i]         = sums0p2.counts.nGenStatus1InCone;
            nGenPtThrInCone0p2[i]             = sums0p2.counts.nGenPtThrInCone;
            nGenStatus1PtThrInCone0p2[i]      = sums0p2.counts.nGenStatus1PtThrInCone;
            nChargedGenInCone0p2[i]         = sums0p2.counts.nChargedGenInCone;
            nChargedGenStatus1InCone0p2[i]  = sums0p2.counts.nChargedGenStatus1InCone;
            nChargedGenPt2InCone0p2[i]      = sums0p2.counts.nChargedGenPt2InCone;
            nChargedGenStatus1Pt2InCone0p2[i]=sums0p2.counts.nChargedGenStatus1Pt2InCone;
            nNeutralGenInCone0p2[i]         = sums0p2.counts.nNeutralGenInCone;
            nNeutralGenStatus1InCone0p2[i]  = sums0p2.counts.nNeutralGenStatus1InCone;
            nNeutralGenPt1InCone0p2[i]      = sums0p2.counts.nNeutralGenPt1InCone;
            nNeutralGenStatus1Pt1InCone0p2[i]=sums0p2.counts.nNeutralGenStatus1Pt1InCone;
            nChargedHadGenInCone0p2[i]      = sums0p2.counts.nChargedHadGenInCone;
            nChargedHadGenStatus1InCone0p2[i]=sums0p2.counts.nChargedHadGenStatus1InCone;
            nChargedHadGenPt2InCone0p2[i]   = sums0p2.counts.nChargedHadGenPt2InCone;
            nChargedHadGenStatus1Pt2InCone0p2[i]=sums0p2.counts.nChargedHadGenStatus1Pt2InCone;
            nNeutralHadGenInCone0p2[i]      = sums0p2.counts.nNeutralHadGenInCone;
            nNeutralHadGenStatus1InCone0p2[i]=sums0p2.counts.nNeutralHadGenStatus1InCone;
            nNeutralHadGenPt1InCone0p2[i]   = sums0p2.counts.nNeutralHadGenPt1InCone;
            nNeutralHadGenStatus1Pt1InCone0p2[i]=sums0p2.counts.nNeutralHadGenStatus1Pt1InCone;

            nGenInCone0p3[i]                = sums0p3.counts.nGenInCone;
            nGenStatus1InCone0p3[i]         = sums0p3.counts.nGenStatus1InCone;
            nGenPtThrInCone0p3[i]             = sums0p3.counts.nGenPtThrInCone;
            nGenStatus1PtThrInCone0p3[i]      = sums0p3.counts.nGenStatus1PtThrInCone;
            nChargedGenInCone0p3[i]         = sums0p3.counts.nChargedGenInCone;
            nChargedGenStatus1InCone0p3[i]  = sums0p3.counts.nChargedGenStatus1InCone;
            nChargedGenPt2InCone0p3[i]      = sums0p3.counts.nChargedGenPt2InCone;
            nChargedGenStatus1Pt2InCone0p3[i]=sums0p3.counts.nChargedGenStatus1Pt2InCone;
            nNeutralGenInCone0p3[i]         = sums0p3.counts.nNeutralGenInCone;
            nNeutralGenStatus1InCone0p3[i]  = sums0p3.counts.nNeutralGenStatus1InCone;
            nNeutralGenPt1InCone0p3[i]      = sums0p3.counts.nNeutralGenPt1InCone;
            nNeutralGenStatus1Pt1InCone0p3[i]=sums0p3.counts.nNeutralGenStatus1Pt1InCone;
            nChargedHadGenInCone0p3[i]      = sums0p3.counts.nChargedHadGenInCone;
            nChargedHadGenStatus1InCone0p3[i]=sums0p3.counts.nChargedHadGenStatus1InCone;
            nChargedHadGenPt2InCone0p3[i]   = sums0p3.counts.nChargedHadGenPt2InCone;
            nChargedHadGenStatus1Pt2InCone0p3[i]=sums0p3.counts.nChargedHadGenStatus1Pt2InCone;
            nNeutralHadGenInCone0p3[i]      = sums0p3.counts.nNeutralHadGenInCone;
            nNeutralHadGenStatus1InCone0p3[i]=sums0p3.counts.nNeutralHadGenStatus1InCone;
            nNeutralHadGenPt1InCone0p3[i]   = sums0p3.counts.nNeutralHadGenPt1InCone;
            nNeutralHadGenStatus1Pt1InCone0p3[i]=sums0p3.counts.nNeutralHadGenStatus1Pt1InCone;
        }

        // ---- add columns for 0p2 and 0p3 cones ----
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

	// ----------------------------------------------------------------------
	// Adds GEN count columns conditionally, based on the column mask.
	// ----------------------------------------------------------------------
	auto addCounts = [&](const std::string &prefix,
			     const std::string &cone,
			     const std::vector<int> &InCone,
			     const std::vector<int> &Status1InCone,
			     const std::vector<int> &PtThrInCone,
			     const std::vector<int> &Pt1InCone,
			     const std::vector<int> &Pt2InCone,
			     const std::vector<int> &Status1Pt1InCone,
			     const GenColumnMask &mask)
	{
	  out->addColumn<int>(prefix+"InCone"+cone,        InCone,        "");
	  out->addColumn<int>(prefix+"Status1InCone"+cone, Status1InCone, "");
	  
	  if (mask.makePtThr)
	    out->addColumn<int>(prefix+"PtThrInCone"+cone, PtThrInCone, "");
	  
	  if (mask.makePt1) {
	    out->addColumn<int>(prefix+"Pt1InCone"+cone, Pt1InCone, "");
	    out->addColumn<int>(prefix+"Status1Pt1InCone"+cone, Status1Pt1InCone, "");
	  }
	  
	  if (mask.makePt2)
	    out->addColumn<int>(prefix+"Pt2InCone"+cone, Pt2InCone, "");
	};

	// -------------- GEN total (PtThr only) -----------------
	{
	  GenColumnMask mask; 
	  mask.makePtThr = true;
	  
	  addCounts("nGen", "0p1",
		    nGenInCone0p1, nGenStatus1InCone0p1,
		    nGenPtThrInCone0p1,
		    {}, {}, 
		    nGenStatus1PtThrInCone0p1, mask);
	  
	  addCounts("nGen", "0p2",
		    nGenInCone0p2, nGenStatus1InCone0p2,
		    nGenPtThrInCone0p2,
		    {}, {}, 
		    nGenStatus1PtThrInCone0p2, mask);
	  
	  addCounts("nGen", "0p3",
		    nGenInCone0p3, nGenStatus1InCone0p3,
		    nGenPtThrInCone0p3,
		    {}, {}, 
		    nGenStatus1PtThrInCone0p3, mask);
	}
	
	// -------------- Charged (Pt2 only) -----------------
	{
	  GenColumnMask mask;
	  mask.makePt2 = true;
	  
	  addCounts("nChargedGen", "0p1",
		    nChargedGenInCone0p1, nChargedGenStatus1InCone0p1,
		    {}, {}, 
		    nChargedGenPt2InCone0p1,
		    {}, mask);
	  
	  addCounts("nChargedGen", "0p2",
		    nChargedGenInCone0p2, nChargedGenStatus1InCone0p2,
		    {}, {}, 
		    nChargedGenPt2InCone0p2,
		    {}, mask);
	  
	  addCounts("nChargedGen", "0p3",
		    nChargedGenInCone0p3, nChargedGenStatus1InCone0p3,
		    {}, {}, 
		    nChargedGenPt2InCone0p3,
		    {}, mask);
	}
	
	// -------------- Neutral (Pt1 only) -----------------
	{
	  GenColumnMask mask;
	  mask.makePt1 = true;
	  
	  addCounts("nNeutralGen", "0p1",
		    nNeutralGenInCone0p1, nNeutralGenStatus1InCone0p1,
		    {}, 
		    nNeutralGenPt1InCone0p1,
		    {}, 
		    nNeutralGenStatus1Pt1InCone0p1, mask);
	  
	  addCounts("nNeutralGen", "0p2",
		    nNeutralGenInCone0p2, nNeutralGenStatus1InCone0p2,
		    {}, 
		    nNeutralGenPt1InCone0p2,
		    {}, 
		    nNeutralGenStatus1Pt1InCone0p2, mask);
	  
	  addCounts("nNeutralGen", "0p3",
		    nNeutralGenInCone0p3, nNeutralGenStatus1InCone0p3,
		    {}, 
		    nNeutralGenPt1InCone0p3,
		    {}, 
		    nNeutralGenStatus1Pt1InCone0p3, mask);
	}
	
	// -------------- Charged Had (Pt2 only) -----------------
	{
	  GenColumnMask mask;
	  mask.makePt2 = true;
	  
	  addCounts("nChargedHadGen", "0p1",
		    nChargedHadGenInCone0p1, nChargedHadGenStatus1InCone0p1,
		    {}, {}, 
		    nChargedHadGenPt2InCone0p1,
		    {}, mask);
	  
	  addCounts("nChargedHadGen", "0p2",
		    nChargedHadGenInCone0p2, nChargedHadGenStatus1InCone0p2,
		    {}, {}, 
		    nChargedHadGenPt2InCone0p2,
		    {}, mask);
	  
	  addCounts("nChargedHadGen", "0p3",
		    nChargedHadGenInCone0p3, nChargedHadGenStatus1InCone0p3,
		    {}, {}, 
		    nChargedHadGenPt2InCone0p3,
		    {}, mask);
	}
	
	// -------------- Neutral Had (Pt1 only) -----------------
	{
	  GenColumnMask mask;
	  mask.makePt1 = true;
	  
	  addCounts("nNeutralHadGen", "0p1",
		    nNeutralHadGenInCone0p1, nNeutralHadGenStatus1InCone0p1,
		    {}, 
		    nNeutralHadGenPt1InCone0p1,
		    {}, 
		    nNeutralHadGenStatus1Pt1InCone0p1, mask);
	  
	  addCounts("nNeutralHadGen", "0p2",
		    nNeutralHadGenInCone0p2, nNeutralHadGenStatus1InCone0p2,
		    {}, 
		    nNeutralHadGenPt1InCone0p2,
		    {}, 
		    nNeutralHadGenStatus1Pt1InCone0p2, mask);
	  
	  addCounts("nNeutralHadGen", "0p3",
		    nNeutralHadGenInCone0p3, nNeutralHadGenStatus1InCone0p3,
		    {}, 
		    nNeutralHadGenPt1InCone0p3,
		    {}, 
		    nNeutralHadGenStatus1Pt1InCone0p3, mask);
	}
	
        // ---- add new gen-count columns ----
        /*auto addCounts = [&](const std::string &prefix, const std::string &cone,
	  const std::vector<int> &a, const std::vector<int> &b,
                             const std::vector<int> &c, const std::vector<int> &d) {
            out->addColumn<int>(prefix+"InCone"+cone, a, "");
            out->addColumn<int>(prefix+"Status1InCone"+cone, b, "");
            out->addColumn<int>(prefix+"PtThrInCone"+cone, c, "");
            out->addColumn<int>(prefix+"Pt1InCone"+cone, c, "");
            out->addColumn<int>(prefix+"Pt2InCone"+cone, c, "");
            out->addColumn<int>(prefix+"Status1Pt1InCone"+cone, d, "");
	    };
	
        addCounts("nGen", "0p1", nGenInCone0p1, nGenStatus1InCone0p1, nGenPtThrInCone0p1, nGenStatus1PtThrInCone0p1);
        addCounts("nGen", "0p2", nGenInCone0p2, nGenStatus1InCone0p2, nGenPtThrInCone0p2, nGenStatus1PtThrInCone0p2);
        addCounts("nGen", "0p3", nGenInCone0p3, nGenStatus1InCone0p3, nGenPtThrInCone0p3, nGenStatus1PtThrInCone0p3);

        addCounts("nChargedGen", "0p1", nChargedGenInCone0p1, nChargedGenStatus1InCone0p1, nChargedGenPt2InCone0p1, nChargedGenStatus1Pt2InCone0p1);
        addCounts("nChargedGen", "0p2", nChargedGenInCone0p2, nChargedGenStatus1InCone0p2, nChargedGenPt2InCone0p2, nChargedGenStatus1Pt2InCone0p2);
        addCounts("nChargedGen", "0p3", nChargedGenInCone0p3, nChargedGenStatus1InCone0p3, nChargedGenPt2InCone0p3, nChargedGenStatus1Pt2InCone0p3);

        addCounts("nNeutralGen", "0p1", nNeutralGenInCone0p1, nNeutralGenStatus1InCone0p1, nNeutralGenPt1InCone0p1, nNeutralGenStatus1Pt1InCone0p1);
        addCounts("nNeutralGen", "0p2", nNeutralGenInCone0p2, nNeutralGenStatus1InCone0p2, nNeutralGenPt1InCone0p2, nNeutralGenStatus1Pt1InCone0p2);
        addCounts("nNeutralGen", "0p3", nNeutralGenInCone0p3, nNeutralGenStatus1InCone0p3, nNeutralGenPt1InCone0p3, nNeutralGenStatus1Pt1InCone0p3);

        addCounts("nChargedHadGen", "0p1", nChargedHadGenInCone0p1, nChargedHadGenStatus1InCone0p1, nChargedHadGenPt2InCone0p1, nChargedHadGenStatus1Pt2InCone0p1);
        addCounts("nChargedHadGen", "0p2", nChargedHadGenInCone0p2, nChargedHadGenStatus1InCone0p2, nChargedHadGenPt2InCone0p2, nChargedHadGenStatus1Pt2InCone0p2);
        addCounts("nChargedHadGen", "0p3", nChargedHadGenInCone0p3, nChargedHadGenStatus1InCone0p3, nChargedHadGenPt2InCone0p3, nChargedHadGenStatus1Pt2InCone0p3);

        addCounts("nNeutralHadGen", "0p1", nNeutralHadGenInCone0p1, nNeutralHadGenStatus1InCone0p1, nNeutralHadGenPt1InCone0p1, nNeutralHadGenStatus1Pt1InCone0p1);
        addCounts("nNeutralHadGen", "0p2", nNeutralHadGenInCone0p2, nNeutralHadGenStatus1InCone0p2, nNeutralHadGenPt1InCone0p2, nNeutralHadGenStatus1Pt1InCone0p2);
        addCounts("nNeutralHadGen", "0p3", nNeutralHadGenInCone0p3, nNeutralHadGenStatus1InCone0p3, nNeutralHadGenPt1InCone0p3, nNeutralHadGenStatus1Pt1InCone0p3);
	*/
       
        // ---- save to event ----
        iEvent.put(std::move(out), cands.coll+"Cands");
        selected.clear();
    }
}

// define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1PFCandTableProducer);

