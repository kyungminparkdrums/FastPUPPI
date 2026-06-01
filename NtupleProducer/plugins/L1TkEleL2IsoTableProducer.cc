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
#include "DataFormats/L1TParticleFlow/interface/PFTrack.h"

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

    std::vector<edm::EDGetTokenT<std::vector<l1t::PFTrack>>> tk_tokens_; // tkIso check
};

L1TkEleL2IsoTableProducer::L1TkEleL2IsoTableProducer(const edm::ParameterSet& iConfig) :
    sel_(iConfig.getParameter<std::string>("commonSel"), true)
{

    produces<nanoaod::FlatTable>(); 
    // I should take these from config but boh
    pf_cands_.emplace_back("L1PFCands", consumes<reco::CandidateView>(edm::InputTag("l1tLayer1:PF")), iConfig);
    tkele_cands_.emplace_back("TkEleL2", consumes<reco::CandidateView>(edm::InputTag("l1tLayer2EG:L1CtTkElectron")), iConfig);

    tk_tokens_.push_back(consumes<std::vector<l1t::PFTrack>>(edm::InputTag("l1tLayer1Barrel", "DecodedTK")));
    tk_tokens_.push_back(consumes<std::vector<l1t::PFTrack>>(edm::InputTag("l1tLayer1HGCal", "DecodedTK")));
}

L1TkEleL2IsoTableProducer::~L1TkEleL2IsoTableProducer() { }

// ------------ method called for each event  ------------
void
L1TkEleL2IsoTableProducer::produce(edm::StreamID id, edm::Event& iEvent, const edm::EventSetup& iSetup) const
{
    edm::Handle<reco::CandidateView> src;

    std::vector<const reco::Candidate *> pf_selected;
    std::vector<const reco::Candidate *> tkele_selected;

    std::vector<const l1t::PFTrack *> tk_selected;

    // for debugging purposes, add everything after each type of veto to the sum
    std::vector<float> vals_isoRawSumAll, vals_isoRawAllSelfVeto, vals_isoRawAllDzVeto, vals_isoRawAllPFeleVeto, vals_isoRawAllPFphoVeto, vals_isoRawAllPFegmVeto;
    std::vector<float> vals_isoRelSumAll, vals_isoRelAllSelfVeto, vals_isoRelAllDzVeto, vals_isoRelAllPFeleVeto, vals_isoRelAllPFphoVeto, vals_isoRelAllPFegmVeto;
    std::vector<float> vals_isoRelSumAllUncorrPt, vals_isoRelAllSelfVetoUncorrPt, vals_isoRelAllDzVetoUncorrPt, vals_isoRelAllPFeleVetoUncorrPt, vals_isoRelAllPFphoVetoUncorrPt, vals_isoRelAllPFegmVetoUncorrPt;

    // consider charged only
    std::vector<float> vals_isoRawSumChg, vals_isoRawChgSelfVeto, vals_isoRawChgDzVeto, vals_isoRawChgPFeleVeto;
    std::vector<float> vals_isoRelSumChg, vals_isoRelChgSelfVeto, vals_isoRelChgDzVeto, vals_isoRelChgPFeleVeto;
    std::vector<float> vals_isoRelSumChgUncorrPt, vals_isoRelChgSelfVetoUncorrPt, vals_isoRelChgDzVetoUncorrPt, vals_isoRelChgPFeleVetoUncorrPt;

    // charged PF with tkIso-like settings
    std::vector<float> vals_isoRawChgTkLike;
    std::vector<float> vals_isoRelChgTkLike;
    std::vector<float> vals_isoRelChgTkLikeUncorrPt;

    // consider neutral only
    std::vector<float> vals_isoRawSumNeu, vals_isoRawNeuSelfVeto, vals_isoRawNeuPFphoVeto;
    std::vector<float> vals_isoRelSumNeu, vals_isoRelNeuSelfVeto, vals_isoRelNeuPFphoVeto;
    std::vector<float> vals_isoRelSumNeuUncorrPt, vals_isoRelNeuSelfVetoUncorrPt, vals_isoRelNeuPFphoVetoUncorrPt;

    // neutrals with (my attempt for) brem-treatment
    std::vector<float> vals_isoRawHybrid;
    std::vector<float> vals_isoRelHybrid;
    std::vector<float> vals_isoRelHybridUncorrPt; 

    // multiplicity
    std::vector<float> vals_nPFall, vals_nPFallDr0p3, vals_nPFallSelfVeto, vals_nPFallDz, vals_nPFallEleVeto, vals_nPFallPhoVeto, vals_nPFallEgmVeto;
    std::vector<float> vals_nPFchg, vals_nPFchgDr0p3, vals_nPFchgSelfVeto, vals_nPFchgDz, vals_nPFchgEleVeto;
    std::vector<float> vals_nPFneu, vals_nPFneuDr0p3, vals_nPFneuSelfVeto, vals_nPFneuPhoVeto;

    // Get PF candidates
    for (const auto & pf_cands : pf_cands_) {
        iEvent.getByToken(pf_cands.src, src);
        pf_selected.reserve(pf_selected.size() + src->size());

        for (const auto & k : *src) {
            if (pf_cands.sel(k)) {
                pf_selected.push_back(&k);
            }
        }
    }

    // Get TkElectron candidates
    for (const auto & tkele_cands : tkele_cands_) {
        iEvent.getByToken(tkele_cands.src, src);
        tkele_selected.reserve(tkele_selected.size() + src->size());

        for (const auto & j : *src) {
            if (tkele_cands.sel(j) && sel_(j)) {
                tkele_selected.push_back(&j);
            }
        }
    }

    // Get decoded L1 tracks used by tk isolation debugging
    for (const auto & tk_token : tk_tokens_) {
        edm::Handle<std::vector<l1t::PFTrack>> tracks;
        iEvent.getByToken(tk_token, tracks);

        if (!tracks.isValid()) continue;

        tk_selected.reserve(tk_selected.size() + tracks->size());

        for (const auto & tk : *tracks) {
            tk_selected.push_back(&tk);
        }
    }

    // create the table
    unsigned int ncands = tkele_selected.size();
    unsigned int ncands_pf = pf_selected.size();

    auto out = std::make_unique<nanoaod::FlatTable>(ncands, "TkEleL2", false, true);

    // resize the vectors per electron candidate size
    vals_isoRawSumAll.resize(ncands);
    vals_isoRawAllSelfVeto.resize(ncands);
    vals_isoRawAllDzVeto.resize(ncands);
    vals_isoRawAllPFeleVeto.resize(ncands);
    vals_isoRawAllPFphoVeto.resize(ncands);
    vals_isoRawAllPFegmVeto.resize(ncands);

    vals_isoRelSumAll.resize(ncands);
    vals_isoRelAllSelfVeto.resize(ncands);
    vals_isoRelAllDzVeto.resize(ncands);
    vals_isoRelAllPFeleVeto.resize(ncands);
    vals_isoRelAllPFphoVeto.resize(ncands);
    vals_isoRelAllPFegmVeto.resize(ncands);

    vals_isoRelSumAllUncorrPt.resize(ncands);
    vals_isoRelAllSelfVetoUncorrPt.resize(ncands);
    vals_isoRelAllDzVetoUncorrPt.resize(ncands);
    vals_isoRelAllPFeleVetoUncorrPt.resize(ncands);
    vals_isoRelAllPFphoVetoUncorrPt.resize(ncands);
    vals_isoRelAllPFegmVetoUncorrPt.resize(ncands);

    vals_isoRawSumChg.resize(ncands);
    vals_isoRawChgSelfVeto.resize(ncands);
    vals_isoRawChgDzVeto.resize(ncands);
    vals_isoRawChgPFeleVeto.resize(ncands);

    vals_isoRelSumChg.resize(ncands);
    vals_isoRelChgSelfVeto.resize(ncands);
    vals_isoRelChgDzVeto.resize(ncands);
    vals_isoRelChgPFeleVeto.resize(ncands);

    vals_isoRelSumChgUncorrPt.resize(ncands);
    vals_isoRelChgSelfVetoUncorrPt.resize(ncands);
    vals_isoRelChgDzVetoUncorrPt.resize(ncands);
    vals_isoRelChgPFeleVetoUncorrPt.resize(ncands);

    vals_isoRawChgTkLike.resize(ncands);
    vals_isoRelChgTkLike.resize(ncands);
    vals_isoRelChgTkLikeUncorrPt.resize(ncands);

    vals_isoRawSumNeu.resize(ncands);
    vals_isoRawNeuSelfVeto.resize(ncands);
    vals_isoRawNeuPFphoVeto.resize(ncands);

    vals_isoRelSumNeu.resize(ncands);
    vals_isoRelNeuSelfVeto.resize(ncands);
    vals_isoRelNeuPFphoVeto.resize(ncands);

    vals_isoRelSumNeuUncorrPt.resize(ncands);
    vals_isoRelNeuSelfVetoUncorrPt.resize(ncands);
    vals_isoRelNeuPFphoVetoUncorrPt.resize(ncands);

    vals_isoRawHybrid.resize(ncands);
    vals_isoRelHybrid.resize(ncands);
    vals_isoRelHybridUncorrPt.resize(ncands);

    vals_nPFall.resize(ncands);
    vals_nPFallDr0p3.resize(ncands);
    vals_nPFallSelfVeto.resize(ncands);
    vals_nPFallDz.resize(ncands);
    vals_nPFallEleVeto.resize(ncands);
    vals_nPFallPhoVeto.resize(ncands);
    vals_nPFallEgmVeto.resize(ncands);

    vals_nPFchg.resize(ncands);
    vals_nPFchgDr0p3.resize(ncands);
    vals_nPFchgSelfVeto.resize(ncands);
    vals_nPFchgDz.resize(ncands);
    vals_nPFchgEleVeto.resize(ncands);

    vals_nPFneu.resize(ncands);
    vals_nPFneuDr0p3.resize(ncands);
    vals_nPFneuSelfVeto.resize(ncands);
    vals_nPFneuPhoVeto.resize(ncands);

    const float bz = 3.8112; // for caloeta/phi calculation

    //bool doPrint = false;
    // loop over electrons
    for (unsigned int iEle = 0; iEle < ncands; iEle++) {
        //if (iEle == 0) doPrint = true;
        

        const auto * tkEle = dynamic_cast<const l1t::TkElectron*>(tkele_selected[iEle]);

        const float eleEta = tkEle->eta();
        const float elePhi = tkEle->phi();
        const float eleZ = tkEle->trkzVtx(); // simple vz() method returns 0 for tkEle object somehow

        const float egEta = tkEle->egCaloPtr()->eta();
        const float egPhi = tkEle->egCaloPtr()->phi();

        // Keep this propagated calo calculation for future reference
        // FIXME: neutral PF matching below uses tkEle->egCaloPtr()->eta()/phi(), but this seems to be different from manually propagated values
        math::XYZTLorentzVector vertex(tkEle->vx(), tkEle->vy(), tkEle->vz(), 0.);
        auto caloetaphi = l1tpf::propagateToCalo(tkEle->p4(), vertex, tkEle->charge(), bz);
        float caloeta = caloetaphi.first;
        float calophi = caloetaphi.second;

        // std::cout << "Electron eta = " << tkEle->eta() << ", caloeta = " << caloeta << ", phi = " << tkEle->phi() << ", calophi = " << calophi << ", default caloeta = " << tkEle->egCaloPtr()->eta() << std::endl;

        const float ptCorr = tkEle->userFloat("ptCorr"); // regressed pt
        const float pt = tkEle->pt();

        //if (doPrint) std::cout << "\nElectron regressed pT = " << ptCorr << ", default pT = " << pt << ", eta = " << eleEta << ", phi = " << elePhi << ", vz = " << eleZ << std::endl;
        
        // raw isolation sums: all PF
        float isoRawSumAll = 0.;
        float isoRawAllSelfVeto = 0.;
        float isoRawAllDzVeto = 0.;
        float isoRawAllPFeleVeto = 0.;
        float isoRawAllPFphoVeto = 0.;
        float isoRawAllPFegmVeto = 0.;

        // raw isolation sums: charged PF
        float isoRawSumChg = 0.;
        float isoRawChgSelfVeto = 0.;
        float isoRawChgDzVeto = 0.;
        float isoRawChgPFeleVeto = 0.;

        float isoRawChgTkLike = 0.;

        // raw isolation sums: neutral PF
        float isoRawSumNeu = 0.;
        float isoRawNeuSelfVeto = 0.;
        float isoRawNeuPFphoVeto = 0.;

        float isoRawHybrid = 0.;

        // multiplicities: all PF
        float nPFall = ncands_pf;
        float nPFallDr0p3 = 0.;
        float nPFallSelfVeto = 0.;
        float nPFallDz = 0.;
        float nPFallEleVeto = 0.;
        float nPFallPhoVeto = 0.;
        float nPFallEgmVeto = 0.;

        // multiplicities: charged PF
        float nPFchg = 0.;
        float nPFchgDr0p3 = 0.;
        float nPFchgSelfVeto = 0.;
        float nPFchgDz = 0.;
        float nPFchgEleVeto = 0.;

        // multiplicities: neutral PF
        float nPFneu = 0.;
        float nPFneuDr0p3 = 0.;
        float nPFneuSelfVeto = 0.;
        float nPFneuPhoVeto = 0.;

        for (const auto * tk : tk_selected) {
            const float tkPt = tk->pt();

            const float tkEta = tk->eta();
            const float tkPhi = tk->phi();
            const float tkVz  = tk->vz();

            const float dR = reco::deltaR(tkEta, tkPhi, eleEta, elePhi);

            if (dR > 0.3f) continue;

            //if (doPrint) std::cout << "L1 Track within dR (e,trk) < 0.3: pT = " << tkPt << ", eta = " << tkEta << ", phi = " << tkPhi << ", vz = " << tkVz << ", dR(e, PF) = " << dR << std::endl;
        }

        // loop over PF candidates
        for (const auto * pf : pf_selected) {
            const float pfPt = pf->pt();
            const int pfCharge = pf->charge();

            const bool isChg = pfCharge != 0;
            const bool isNeu = pfCharge == 0;

            if (isChg) nPFchg++;
            if (isNeu) nPFneu++;

            // charged PF: use TkElectron track eta/phi
            // neutral PF: use EG calo eta/phi
            const float refEta = isChg ? eleEta : egEta;
            const float refPhi = isChg ? elePhi : egPhi;

            const float dR_ele_pf = reco::deltaR(pf->eta(), pf->phi(), refEta, refPhi);

            // tkIso-like charged PF iso 
            if (isChg) {
                bool passTkLike = true;

                if (pfPt < 2.) passTkLike = false;
                if (dR_ele_pf < 0.03f) passTkLike = false;
                if (dR_ele_pf > 0.20f) passTkLike = false;

                const float dzTkLike = std::abs(eleZ - pf->vz());

                if (dzTkLike > 0.6f) passTkLike = false;

                if (passTkLike) isoRawChgTkLike += pfPt;
            }

            // special treatment for neutrals for brem treatment
            bool addToHybrid = false;

            if (isChg) {
                const float dz = std::abs(eleZ - pf->vz());

                addToHybrid = (dR_ele_pf > 0.05f) && (dR_ele_pf < 0.30f) && (dz < 0.5f);
            } 
            else {
                const bool bremsLike = (std::abs(pf->eta() - egEta) < 0.03f) && (std::abs(reco::deltaPhi(pf->phi(), egPhi)) < 0.30f); // some loose definition
                //const bool bremsLike = (std::abs(pf->eta() - egEta) < 0.03f) && (std::abs(reco::deltaPhi(pf->phi(), egPhi)) < 0.30f) && ((pfPt / std::max(ptCorr, 0.1f)) < 1.0f); // some loose definition

                addToHybrid = (dR_ele_pf > 0.05f) && (dR_ele_pf < 0.30f) && (!bremsLike);
            }
            if (addToHybrid) isoRawHybrid += pfPt; 

            // Usual isolations
            if (dR_ele_pf > 0.3f) continue;

            nPFallDr0p3++;
            isoRawSumAll += pfPt;

            if (isChg) {
                //if (doPrint) std::cout << "Charged PF within dR (e,PF) < 0.3: pT = " << pfPt << ", eta = " << pf->eta() << ", phi = " << pf->phi() << ", vz = " << pf->vz() << ", dR(e, PF) = " << dR_ele_pf << std::endl;
                nPFchgDr0p3++;
                isoRawSumChg += pfPt;
            } else {
                nPFneuDr0p3++;
                isoRawSumNeu += pfPt;
            }

            // self-veto
            if (dR_ele_pf < 0.05) continue;

            nPFallSelfVeto++;
            isoRawAllSelfVeto += pfPt;

            if (isChg) {
                nPFchgSelfVeto++;
                isoRawChgSelfVeto += pfPt;
            } else {
                nPFneuSelfVeto++;
                isoRawNeuSelfVeto += pfPt;
            }

            // same vertex requirement for charged PFs only
            bool passDz = true;
            if (isChg) {
                const float dz = std::abs(eleZ - pf->vz());
                passDz = dz <= 0.5;
            }

            if (!passDz) continue;

            nPFallDz++;
            isoRawAllDzVeto += pfPt;

            if (isChg) {
                nPFchgDz++;
                isoRawChgDzVeto += pfPt;
            }

            // PF electron / photon / EGM veto
            const int absPdgId = std::abs(pf->pdgId());
            const bool isPFEle = absPdgId == 11;
            const bool isPFPho = absPdgId == 22;
            const bool isPFEgm = isPFEle || isPFPho;

            // PF electron veto
            if (!isPFEle) {
                nPFallEleVeto++;
                isoRawAllPFeleVeto += pfPt;

                if (isChg) {
                    nPFchgEleVeto++;
                    isoRawChgPFeleVeto += pfPt;
                }
            }

            // PF photon veto
            if (!isPFPho) {
                nPFallPhoVeto++;
                isoRawAllPFphoVeto += pfPt;

                if (isNeu) {
                    nPFneuPhoVeto++;
                    isoRawNeuPFphoVeto += pfPt;
                }
            }

            // PF EGM veto
            if (!isPFEgm) {
                nPFallEgmVeto++;
                isoRawAllPFegmVeto += pfPt;
            }
        }

        // save multiplicities
        vals_nPFall[iEle] = nPFall;
        vals_nPFallDr0p3[iEle] = nPFallDr0p3;
        vals_nPFallSelfVeto[iEle] = nPFallSelfVeto;
        vals_nPFallDz[iEle] = nPFallDz;
        vals_nPFallEleVeto[iEle] = nPFallEleVeto;
        vals_nPFallPhoVeto[iEle] = nPFallPhoVeto;
        vals_nPFallEgmVeto[iEle] = nPFallEgmVeto;

        vals_nPFchg[iEle] = nPFchg;
        vals_nPFchgDr0p3[iEle] = nPFchgDr0p3;
        vals_nPFchgSelfVeto[iEle] = nPFchgSelfVeto;
        vals_nPFchgDz[iEle] = nPFchgDz;
        vals_nPFchgEleVeto[iEle] = nPFchgEleVeto;

        vals_nPFneu[iEle] = nPFneu;
        vals_nPFneuDr0p3[iEle] = nPFneuDr0p3;
        vals_nPFneuSelfVeto[iEle] = nPFneuSelfVeto;
        vals_nPFneuPhoVeto[iEle] = nPFneuPhoVeto;

        // save raw isolation
        vals_isoRawSumAll[iEle] = isoRawSumAll;
        vals_isoRawAllSelfVeto[iEle] = isoRawAllSelfVeto;
        vals_isoRawAllDzVeto[iEle] = isoRawAllDzVeto;
        vals_isoRawAllPFeleVeto[iEle] = isoRawAllPFeleVeto;
        vals_isoRawAllPFphoVeto[iEle] = isoRawAllPFphoVeto;
        vals_isoRawAllPFegmVeto[iEle] = isoRawAllPFegmVeto;

        vals_isoRawSumChg[iEle] = isoRawSumChg;
        vals_isoRawChgSelfVeto[iEle] = isoRawChgSelfVeto;
        vals_isoRawChgDzVeto[iEle] = isoRawChgDzVeto;
        vals_isoRawChgPFeleVeto[iEle] = isoRawChgPFeleVeto;

        vals_isoRawChgTkLike[iEle] = isoRawChgTkLike;

        vals_isoRawSumNeu[iEle] = isoRawSumNeu;
        vals_isoRawNeuSelfVeto[iEle] = isoRawNeuSelfVeto;
        vals_isoRawNeuPFphoVeto[iEle] = isoRawNeuPFphoVeto;

        vals_isoRawHybrid[iEle] = isoRawHybrid;
        
        // relative isolation using corrected pt
        vals_isoRelSumAll[iEle] = isoRawSumAll / ptCorr;
        vals_isoRelAllSelfVeto[iEle] = isoRawAllSelfVeto / ptCorr;
        vals_isoRelAllDzVeto[iEle] = isoRawAllDzVeto / ptCorr;
        vals_isoRelAllPFeleVeto[iEle] = isoRawAllPFeleVeto / ptCorr;
        vals_isoRelAllPFphoVeto[iEle] = isoRawAllPFphoVeto / ptCorr;
        vals_isoRelAllPFegmVeto[iEle] = isoRawAllPFegmVeto / ptCorr;

        vals_isoRelSumChg[iEle] = isoRawSumChg / ptCorr;
        vals_isoRelChgSelfVeto[iEle] = isoRawChgSelfVeto / ptCorr;
        vals_isoRelChgDzVeto[iEle] = isoRawChgDzVeto / ptCorr;
        vals_isoRelChgPFeleVeto[iEle] = isoRawChgPFeleVeto / ptCorr;

        vals_isoRelChgTkLike[iEle] = isoRawChgTkLike / ptCorr;

        vals_isoRelSumNeu[iEle] = isoRawSumNeu / ptCorr;
        vals_isoRelNeuSelfVeto[iEle] = isoRawNeuSelfVeto / ptCorr;
        vals_isoRelNeuPFphoVeto[iEle] = isoRawNeuPFphoVeto / ptCorr;

        vals_isoRelHybrid[iEle] = isoRawHybrid / ptCorr;

        // relative isolation using uncorrected pt
        vals_isoRelSumAllUncorrPt[iEle] = isoRawSumAll / pt;
        vals_isoRelAllSelfVetoUncorrPt[iEle] = isoRawAllSelfVeto / pt;
        vals_isoRelAllDzVetoUncorrPt[iEle] = isoRawAllDzVeto / pt;
        vals_isoRelAllPFeleVetoUncorrPt[iEle] = isoRawAllPFeleVeto / pt;
        vals_isoRelAllPFphoVetoUncorrPt[iEle] = isoRawAllPFphoVeto / pt;
        vals_isoRelAllPFegmVetoUncorrPt[iEle] = isoRawAllPFegmVeto / pt;

        vals_isoRelSumChgUncorrPt[iEle] = isoRawSumChg / pt;
        vals_isoRelChgSelfVetoUncorrPt[iEle] = isoRawChgSelfVeto / pt;
        vals_isoRelChgDzVetoUncorrPt[iEle] = isoRawChgDzVeto / pt;
        vals_isoRelChgPFeleVetoUncorrPt[iEle] = isoRawChgPFeleVeto / pt;

        vals_isoRelChgTkLikeUncorrPt[iEle] = isoRawChgTkLike / pt;

        vals_isoRelSumNeuUncorrPt[iEle] = isoRawSumNeu / pt;
        vals_isoRelNeuSelfVetoUncorrPt[iEle] = isoRawNeuSelfVeto / pt;
        vals_isoRelNeuPFphoVetoUncorrPt[iEle] = isoRawNeuPFphoVeto / pt;
        
        vals_isoRelHybridUncorrPt[iEle] = isoRawHybrid / pt;
    }

    // multiplicity branches
    out->addColumn<float>("nPFall", vals_nPFall, "number of PF candidates in the event");
    out->addColumn<float>("nPFallDr0p3", vals_nPFallDr0p3, "number of PF candidates within dR < 0.3");
    out->addColumn<float>("nPFallSelfVeto", vals_nPFallSelfVeto, "number of PF candidates within dR < 0.3 after self veto");
    out->addColumn<float>("nPFallDz", vals_nPFallDz, "number of PF candidates within dR < 0.3 after self veto and dz veto");
    out->addColumn<float>("nPFallEleVeto", vals_nPFallEleVeto, "number of PF candidates after PF electron veto");
    out->addColumn<float>("nPFallPhoVeto", vals_nPFallPhoVeto, "number of PF candidates after PF photon veto");
    out->addColumn<float>("nPFallEgmVeto", vals_nPFallEgmVeto, "number of PF candidates after PF EGM veto");

    out->addColumn<float>("nPFchg", vals_nPFchg, "number of charged PF candidates in the event");
    out->addColumn<float>("nPFchgDr0p3", vals_nPFchgDr0p3, "number of charged PF candidates within dR < 0.3");
    out->addColumn<float>("nPFchgSelfVeto", vals_nPFchgSelfVeto, "number of charged PF candidates after self veto");
    out->addColumn<float>("nPFchgDz", vals_nPFchgDz, "number of charged PF candidates after dz veto");
    out->addColumn<float>("nPFchgEleVeto", vals_nPFchgEleVeto, "number of charged PF candidates after PF electron veto");

    out->addColumn<float>("nPFneu", vals_nPFneu, "number of neutral PF candidates in the event");
    out->addColumn<float>("nPFneuDr0p3", vals_nPFneuDr0p3, "number of neutral PF candidates within dR < 0.3");
    out->addColumn<float>("nPFneuSelfVeto", vals_nPFneuSelfVeto, "number of neutral PF candidates after self veto");
    out->addColumn<float>("nPFneuPhoVeto", vals_nPFneuPhoVeto, "number of neutral PF candidates after PF photon veto");

    // raw isolation branches: all PF
    out->addColumn<float>("customPfIsoRawSumAll", vals_isoRawSumAll, "custom PF iso raw, all PF, no veto");
    out->addColumn<float>("customPfIsoRawAllSelfVeto", vals_isoRawAllSelfVeto, "custom PF iso raw, all PF, self veto");
    out->addColumn<float>("customPfIsoRawAllDzVeto", vals_isoRawAllDzVeto, "custom PF iso raw, all PF, dz veto");
    out->addColumn<float>("customPfIsoRawAllPFeleVeto", vals_isoRawAllPFeleVeto, "custom PF iso raw, all PF, PF electron veto");
    out->addColumn<float>("customPfIsoRawAllPFphoVeto", vals_isoRawAllPFphoVeto, "custom PF iso raw, all PF, PF photon veto");
    out->addColumn<float>("customPfIsoRawAllPFegmVeto", vals_isoRawAllPFegmVeto, "custom PF iso raw, all PF, PF EGM veto");

    // relative isolation branches: all PF, corrected pt
    out->addColumn<float>("customPfIsoRelSumAll", vals_isoRelSumAll, "custom PF iso relative, all PF, no veto, corrected pt");
    out->addColumn<float>("customPfIsoRelAllSelfVeto", vals_isoRelAllSelfVeto, "custom PF iso relative, all PF, self veto, corrected pt");
    out->addColumn<float>("customPfIsoRelAllDzVeto", vals_isoRelAllDzVeto, "custom PF iso relative, all PF, dz veto, corrected pt");
    out->addColumn<float>("customPfIsoRelAllPFeleVeto", vals_isoRelAllPFeleVeto, "custom PF iso relative, all PF, PF electron veto, corrected pt");
    out->addColumn<float>("customPfIsoRelAllPFphoVeto", vals_isoRelAllPFphoVeto, "custom PF iso relative, all PF, PF photon veto, corrected pt");
    out->addColumn<float>("customPfIsoRelAllPFegmVeto", vals_isoRelAllPFegmVeto, "custom PF iso relative, all PF, PF EGM veto, corrected pt");

    // relative isolation branches: all PF, uncorrected pt
    out->addColumn<float>("customPfIsoRelSumAllUncorrPt", vals_isoRelSumAllUncorrPt, "custom PF iso relative, all PF, no veto, uncorrected pt");
    out->addColumn<float>("customPfIsoRelAllSelfVetoUncorrPt", vals_isoRelAllSelfVetoUncorrPt, "custom PF iso relative, all PF, self veto, uncorrected pt");
    out->addColumn<float>("customPfIsoRelAllDzVetoUncorrPt", vals_isoRelAllDzVetoUncorrPt, "custom PF iso relative, all PF, dz veto, uncorrected pt");
    out->addColumn<float>("customPfIsoRelAllPFeleVetoUncorrPt", vals_isoRelAllPFeleVetoUncorrPt, "custom PF iso relative, all PF, PF electron veto, uncorrected pt");
    out->addColumn<float>("customPfIsoRelAllPFphoVetoUncorrPt", vals_isoRelAllPFphoVetoUncorrPt, "custom PF iso relative, all PF, PF photon veto, uncorrected pt");
    out->addColumn<float>("customPfIsoRelAllPFegmVetoUncorrPt", vals_isoRelAllPFegmVetoUncorrPt, "custom PF iso relative, all PF, PF EGM veto, uncorrected pt");

    // charged-only branches
    out->addColumn<float>("customPfIsoRawSumChg", vals_isoRawSumChg, "custom PF iso raw, charged PF, no veto");
    out->addColumn<float>("customPfIsoRawChgSelfVeto", vals_isoRawChgSelfVeto, "custom PF iso raw, charged PF, self veto");
    out->addColumn<float>("customPfIsoRawChgDzVeto", vals_isoRawChgDzVeto, "custom PF iso raw, charged PF, dz veto");
    out->addColumn<float>("customPfIsoRawChgPFeleVeto", vals_isoRawChgPFeleVeto, "custom PF iso raw, charged PF, PF electron veto");

    out->addColumn<float>("customPfIsoRelSumChg", vals_isoRelSumChg, "custom PF iso relative, charged PF, no veto, corrected pt");
    out->addColumn<float>("customPfIsoRelChgSelfVeto", vals_isoRelChgSelfVeto, "custom PF iso relative, charged PF, self veto, corrected pt");
    out->addColumn<float>("customPfIsoRelChgDzVeto", vals_isoRelChgDzVeto, "custom PF iso relative, charged PF, dz veto, corrected pt");
    out->addColumn<float>("customPfIsoRelChgPFeleVeto", vals_isoRelChgPFeleVeto, "custom PF iso relative, charged PF, PF electron veto, corrected pt");

    out->addColumn<float>("customPfIsoRelSumChgUncorrPt", vals_isoRelSumChgUncorrPt, "custom PF iso relative, charged PF, no veto, uncorrected pt");
    out->addColumn<float>("customPfIsoRelChgSelfVetoUncorrPt", vals_isoRelChgSelfVetoUncorrPt, "custom PF iso relative, charged PF, self veto, uncorrected pt");
    out->addColumn<float>("customPfIsoRelChgDzVetoUncorrPt", vals_isoRelChgDzVetoUncorrPt, "custom PF iso relative, charged PF, dz veto, uncorrected pt");
    out->addColumn<float>("customPfIsoRelChgPFeleVetoUncorrPt", vals_isoRelChgPFeleVetoUncorrPt, "custom PF iso relative, charged PF, PF electron veto, uncorrected pt");

    // tkIso-like
    out->addColumn<float>("customPfIsoRawChgTkLike", vals_isoRawChgTkLike, "charged PF isolation with tkIso-like selections");
    out->addColumn<float>("customPfIsoRelChgTkLike", vals_isoRelChgTkLike, "charged PF isolation with tkIso-like selections, corrected pt");
    out->addColumn<float>("customPfIsoRelChgTkLikeUncorrPt", vals_isoRelChgTkLikeUncorrPt, "charged PF isolation with tkIso-like selections, uncorrected pt");

    // neutral-only branches
    out->addColumn<float>("customPfIsoRawSumNeu", vals_isoRawSumNeu, "custom PF iso raw, neutral PF, no veto");
    out->addColumn<float>("customPfIsoRawNeuSelfVeto", vals_isoRawNeuSelfVeto, "custom PF iso raw, neutral PF, self veto");
    out->addColumn<float>("customPfIsoRawNeuPFphoVeto", vals_isoRawNeuPFphoVeto, "custom PF iso raw, neutral PF, PF photon veto");

    out->addColumn<float>("customPfIsoRelSumNeu", vals_isoRelSumNeu, "custom PF iso relative, neutral PF, no veto, corrected pt");
    out->addColumn<float>("customPfIsoRelNeuSelfVeto", vals_isoRelNeuSelfVeto, "custom PF iso relative, neutral PF, self veto, corrected pt");
    out->addColumn<float>("customPfIsoRelNeuPFphoVeto", vals_isoRelNeuPFphoVeto, "custom PF iso relative, neutral PF, PF photon veto, corrected pt");

    out->addColumn<float>("customPfIsoRelSumNeuUncorrPt", vals_isoRelSumNeuUncorrPt, "custom PF iso relative, neutral PF, no veto, uncorrected pt");
    out->addColumn<float>("customPfIsoRelNeuSelfVetoUncorrPt", vals_isoRelNeuSelfVetoUncorrPt, "custom PF iso relative, neutral PF, self veto, uncorrected pt");
    out->addColumn<float>("customPfIsoRelNeuPFphoVetoUncorrPt", vals_isoRelNeuPFphoVetoUncorrPt, "custom PF iso relative, neutral PF, PF photon veto, uncorrected pt");

    out->addColumn<float>("customPfIsoRawHybrid", vals_isoRawHybrid, "PF isolation with brems veto");
    out->addColumn<float>("customPfIsoRelHybrid", vals_isoRelHybrid, "PF isolation with brems veto, corrected pt");
    out->addColumn<float>("customPfIsoRelHybridUncorrPt", vals_isoRelHybridUncorrPt, "PF isolation with brems veto, uncorrected pt");

    // save to the event branches
    iEvent.put(std::move(out));
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1TkEleL2IsoTableProducer);
