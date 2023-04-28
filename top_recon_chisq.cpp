#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>
#include <algorithm>
#include <string>

#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TH1F.h"

// DECAY PROCESS
/* p p -> t t~
 * t -> b W+
 * W+ -> q q~
 * t~ -> b~ W-
 * W- -> q q~ */

// EFFICIENCY DEFINITIONS
/* 3 efficiencies: eff_1, eff_2, eff_event
 * 3 numerators (e_1, e_2, e_event), 2 denominators (e_1, e_2)
 * Denominators:
 *  subscript corresponds to number of matched (matched to truth particle) tops
 *  e_1 means only one top quark was matched, e_2 means both top quarks were matched
 * Numerators:
 *  subscript corresponds to number of top quarks reconstructed correctly using chisq method
 *  e_1 means one top quark was reconstructed correctly (for an e_1 event), e_2 means one top quark was reconstructed correctly (for an e_2 event), e_event means both top quarks were reconstructed correctly (for an e_2 event) */

TFile *f = TFile::Open("/Users/divitrawal/MG5_aMC_v3_4_0/TOPRECON/Events/run_01/tag_1_delphes_events.root","read");
TTree *t = (TTree*)f -> Get("Delphes");

TBranch *bJetBTag = t -> GetBranch("Jet.BTag");
TBranch *bJetPT = t -> GetBranch("Jet.PT");
TBranch *bJetEta = t -> GetBranch("Jet.Eta");
TBranch *bJetPhi = t -> GetBranch("Jet.Phi");
TBranch *bJetM = t -> GetBranch("Jet.Mass");
TBranch *bJetNum = t -> GetBranch("Jet_size");
TBranch *bParticlePT = t -> GetBranch("Particle.PT");
TBranch *bParticleEta = t -> GetBranch("Particle.Eta");
TBranch *bParticlePhi = t -> GetBranch("Particle.Phi");
TBranch *bParticleM = t -> GetBranch("Particle.Mass");
TBranch *bParticleID = t -> GetBranch("Particle.PID");
TBranch *bParticleStatus = t -> GetBranch("Particle.Status");
TBranch *bParticleMother = t -> GetBranch("Particle.M1");

// returns an ordered list of particles formatted as [[b,q,q~],[b~,q',q'~]]
// note that each particle is a tuple of its four-vector and btag
std::vector<std::vector<std::tuple<TLorentzVector, bool>>> get_particles()
{
    std::vector<int> unordered_particles; // holds all final state particle indexes
    for (int p = 0; p < t -> GetLeaf("Particle_size") -> GetValue(); ++p)
    {
        if (bParticleStatus -> GetLeaf("Particle.Status") -> GetValue(p) == 23) // is final state particle
        {
            unordered_particles.push_back(p);
        }
    }
    std::vector<std::tuple<TLorentzVector, bool>> t1;
    std::vector<std::tuple<TLorentzVector, bool>> t2;
    for (int q : unordered_particles) // add bottom quark to t1
    {
        if (bParticleID -> GetLeaf("Particle.PID") -> GetValue(q) == 5) // is bottom quark
        {
            TLorentzVector particle;
            particle.SetPtEtaPhiM(bParticlePT -> GetLeaf("Particle.PT") -> GetValue(q),
                                  bParticleEta -> GetLeaf("Particle.Eta") -> GetValue(q),
                                  bParticlePhi -> GetLeaf("Particle.Phi") -> GetValue(q),
                                  bParticleM -> GetLeaf("Particle.Mass") -> GetValue(q));
            std::tuple<TLorentzVector, bool> particle_tuple = {particle, true};
            t1.push_back(particle_tuple);
        }
    }
    for (int q : unordered_particles) // add q and q~ to t1
    {
        if (bParticleID -> GetLeaf("Particle.PID") -> GetValue(bParticleMother -> GetLeaf("Particle.M1") -> GetValue(q)) == 24) // originated from W+ boson
        {
            TLorentzVector particle;
            particle.SetPtEtaPhiM(bParticlePT -> GetLeaf("Particle.PT") -> GetValue(q),
                                  bParticleEta -> GetLeaf("Particle.Eta") -> GetValue(q),
                                  bParticlePhi -> GetLeaf("Particle.Phi") -> GetValue(q),
                                  bParticleM -> GetLeaf("Particle.Mass") -> GetValue(q));
            std::tuple<TLorentzVector, bool> particle_tuple = {particle, false};
            t1.push_back(particle_tuple);
        }
    }
    for (int q : unordered_particles)
    {
        if (bParticleID -> GetLeaf("Particle.PID") -> GetValue(q) == -5) // is anti-bottom quark
        {
            TLorentzVector particle;
            particle.SetPtEtaPhiM(bParticlePT -> GetLeaf("Particle.PT") -> GetValue(q),
                                  bParticleEta -> GetLeaf("Particle.Eta") -> GetValue(q),
                                  bParticlePhi -> GetLeaf("Particle.Phi") -> GetValue(q),
                                  bParticleM -> GetLeaf("Particle.Mass") -> GetValue(q));
            std::tuple<TLorentzVector, bool> particle_tuple = {particle, true};
            t2.push_back(particle_tuple);
        }
    }
    for (int q : unordered_particles)
    {
        if (bParticleID -> GetLeaf("Particle.PID") -> GetValue(bParticleMother -> GetLeaf("Particle.M1") -> GetValue(q)) == -24) // originated from W- boson
        {
            TLorentzVector particle;
            particle.SetPtEtaPhiM(bParticlePT -> GetLeaf("Particle.PT") -> GetValue(q),
                                  bParticleEta -> GetLeaf("Particle.Eta") -> GetValue(q),
                                  bParticlePhi -> GetLeaf("Particle.Phi") -> GetValue(q),
                                  bParticleM -> GetLeaf("Particle.Mass") -> GetValue(q));
            std::tuple<TLorentzVector, bool> particle_tuple = {particle, false};
            t2.push_back(particle_tuple);
        }
    }
    std::vector<std::vector<std::tuple<TLorentzVector, bool>>> ordered_particles = {t1,t2};
    for (std::vector<std::tuple<TLorentzVector, bool>> t : ordered_particles)
    {
        std::cout << "TOP QUARK" << std::endl;
        for (std::tuple<TLorentzVector, bool> p : t)
        {
            std::cout << "M: " << std::get<0>(p).M() << " PT: " << std::get<0>(p).Pt() << " Btag: " << std::get<1>(p) << std::endl;
        }
    }

    return ordered_particles;
}

// returns vector containing number of jets, number of btagged jets, and list of tuples (each vector contains four-vector and btag)
std::tuple<int, int, std::vector<std::tuple<TLorentzVector, bool>>> get_jets()
{
    std::vector<std::tuple<TLorentzVector, bool>> unordered_jets;
    int num_btagged = 0;
    for (int j = 0; j < bJetNum -> GetLeaf("Jet_size") -> GetValue(); ++j)
    {
        TLorentzVector jet;
        jet.SetPtEtaPhiM(bJetPT -> GetLeaf("Jet.PT") -> GetValue(j),
                         bJetEta -> GetLeaf("Jet.Eta") -> GetValue(j),
                         bJetPhi -> GetLeaf("Jet.Phi") -> GetValue(j),
                         bJetM -> GetLeaf("Jet.Mass") -> GetValue(j));
        bool btag = bJetBTag -> GetLeaf("Jet.BTag") -> GetValue(j);
        if (btag == 1)
        {
            num_btagged++;
        }
        std::tuple<TLorentzVector, bool> t = std::make_tuple(jet, btag);
        unordered_jets.push_back(t);
    }
    return std::make_tuple(unordered_jets.size(), num_btagged, unordered_jets);
}

// returns index of element in vector
int indexOf(std::vector<TLorentzVector> v, TLorentzVector element)
{
    auto idx = std::find(v.begin(), v.end(), element);
    return std::distance(v.begin(), idx);
}

bool contains(std::vector<std::tuple<TLorentzVector, bool>> list, std::tuple<TLorentzVector, bool> element)
{
    if (std::find(list.begin(), list.end(), element) != list.end()) 
    {
        return true;
    } 
    return false;
}

/*
// returns list of jets within angular distance 0.4 for each particle, structured as [[[(four-vector,btag),(...),(...)],[(...),(...),(...)]]]
std::vector<std::vector<std::vector<std::tuple<TLorentzVector, bool>>>> match(std::vector<std::vector<TLorentzVector>> ordered_particles, std::vector<std::tuple<TLorentzVector, bool>> unordered_jets)
{
    std::vector<std::vector<std::vector<std::tuple<TLorentzVector, bool>>>> ordered_jets;
    for (std::vector<TLorentzVector> quark : ordered_particles)
    {
        for (TLorentzVector particle : quark)
        {
            std::vector<double> dR_values;
            std::vector<std::tuple<TLorentzVector, bool>> dR_jets;
            if (indexOf(quark, particle) == 0) // particle is bottom or anti-bottom quark
            {
                for (std::tuple<TLorentzVector, bool> jet: unordered_jets)
                {
                    double del_r = particle.deltaR(std::get<0>(jet));
                    if (del_r < 0.4 && std::get<1>(jet) == true) // within angular distance 0.4 and is btagged
                    {
                        dR_values.push_back(del_r);
                        dR_jets.push_back(jet);
                    }
                }
            }
            else
            {
                for (std::tuple<TLorentzVector, bool> jet: unordered_jets)
                {
                    double del_r = particle.deltaR(std::get<0>(jet));
                    if (del_r < 0.4) // within angular distance 0.4
                    {
                        dR_values.push_back(del_r);
                        dR_jets.push_back(jet);
                    }
                }
            }
            if (dR_values.size() == 0)
            {
                TLorentzVector v;
                ordered_jets.push_back(v);
            }
            else
            {
                double min_dR = std::min_element(dR_values.begin(), dR_values.end());
                int index = std::distance(dR_values.begin(), std::find(dR_values.begin(), dR_values.end(), min_dR));
                ordered_jets.push_back(dR_jets[index]);
            }
        }







            for (std::tuple<TLorentzVector, bool> jet : unordered_jets)
            {
                double del_r = particle.DeltaR(std::get<0>(jet))
                if (del_r < 0.4)
                {
                    dR_values.push_back(del_r);
                    dR_jets.push_back(jet);
                }
            }
            if (dR_values.size() == 0)
            {
                TLorentzVector v;
                ordered_jets.push_back(std::make_tuple(v, false));
            }
        }

    }

*/

/*
std::vector<std::vector<std::vector<std::tuple<TLorentzVector, bool>>>> all_corresponding_jets
{
    for (std::vector<TLorentzVector> quark : ordered_particles)
    {
        std::vector<std::vector<std::tuple<TLorentzVector, bool>>> corresponding_jets;
        for (TLorentzVector p : quark)
        {
            std::vector<std::tuple<TLorentzVector, bool>> matched_jets;
            if (indexOf(quark,p) == 0) // is bottom quark
            {
                for (std::tuple<TLorentzVector, bool> j : unordered_jets)
                {
                    if (p.DeltaR(std::get<0>(j)) < 0.4 && std::get<1>(j)) // jet is within angular distance of 0.4 and btagged
                    {
                        matched_jets.push_back(j);
                    }
                }
            }
            else
            {
                for(std::tuple<TLorentzVector, bool> j : unordered_jets)
                {
                    if (p.DeltaR(std::get<0>(j)) < 0.4)
                    {
                        matched_jets.push_back(j);
                    }
                }
            }
            corresponding_jets.push_back(matched_jets);
        }
        all_corresponding_jets.push_back(corresponding_jets);
    }
    return all_corresponding_jets;
}
*/

// returns list of ordered jets formatted as [[b,q,q~],[b~,q',q'~]]
std::vector<std::vector<std::tuple<TLorentzVector, bool>>> match(std::vector<std::vector<std::tuple<TLorentzVector, bool>>> ordered_particles, std::vector<std::tuple<TLorentzVector, bool>> unordered_jets)
{
    std::vector<std::vector<std::tuple<TLorentzVector, bool>>> ordered_jets;
    std::vector<std::tuple<TLorentzVector, bool>> used_jets;
    for (std::vector<std::tuple<TLorentzVector, bool>> quark : ordered_particles)
    {
        std::vector<std::tuple<TLorentzVector, bool>> top_quark_jets;
        for (std::tuple<TLorentzVector, bool> parton : quark)
        {   
            if(std::get<1>(parton)) // particle is b or b~
            {
                double min_dR = 99999999;
                std::tuple<TLorentzVector, bool> min_dR_jet;
                for (std::tuple<TLorentzVector, bool> jet : unordered_jets)
                {
                    if (!contains(used_jets, jet) && std::get<1>(jet)) // jet has not been used yet and is btagged
                    {
                        double del_R = std::get<0>(parton).DeltaR(std::get<0>(jet));
                        if (del_R < min_dR)
                        {
                            min_dR = del_R;
                            min_dR_jet = jet;
                        }
                    }
                }
                if (min_dR < 0.4) // angular distance less than 0.4
                {
                    top_quark_jets.push_back(min_dR_jet);
                    used_jets.push_back(min_dR_jet);
                }
                else // add dummy jet
                {
                    TLorentzVector v;
                    v.SetPtEtaPhiM(0,0,0,0);
                    std::tuple<TLorentzVector, bool> v_tuple = {v, false};
                    top_quark_jets.push_back(v_tuple);
                }
            }
            else // particle is not b or b~
            {
                double min_dR = 99999999;
                std::tuple<TLorentzVector, bool> min_dR_jet;
                for (std::tuple<TLorentzVector, bool> jet : unordered_jets)
                {
                    if (!contains(used_jets, jet) && !std::get<1>(jet)) // jet has not been used yet and is not btagged
                    {
                        double del_R = std::get<0>(parton).DeltaR(std::get<0>(jet));
                        if (del_R < min_dR)
                        {
                            min_dR = del_R;
                            min_dR_jet = jet;
                        }
                    }
                }
                if (min_dR < 0.4) // angular distance less than 0.4
                {
                    top_quark_jets.push_back(min_dR_jet);
                    used_jets.push_back(min_dR_jet);
                }
                else // add dummy jet
                {
                    TLorentzVector v;
                    v.SetPtEtaPhiM(0,0,0,0);
                    std::tuple<TLorentzVector, bool> v_tuple = {v, false};
                    top_quark_jets.push_back(v_tuple);
                }
            }
        }
        ordered_jets.push_back(top_quark_jets);
    }
    for (std::vector<std::tuple<TLorentzVector, bool>> q : ordered_jets)
    {
        std::cout << "MATCHED QUARK" << std::endl;
        for (std::tuple<TLorentzVector, bool> j : q)
        {
            std::cout << "M: " << std::get<0>(j).M() << " PT: " << std::get<0>(j).Pt() << " Btag: " << std::get<1>(j) << std::endl;
        }
    }
    return ordered_jets;
}

// returns which top quarks are unambiguously matched
std::tuple<bool, bool, bool> unambiguously_matched(std::vector<std::vector<std::tuple<TLorentzVector, bool>>> ordered_jets)
{
    std::vector<bool> um_vector;
    TLorentzVector dummy_jet;
    dummy_jet.SetPtEtaPhiM(0,0,0,0);
    for (std::vector<std::tuple<TLorentzVector, bool>> quark : ordered_jets)
    {
        bool all_matched = true;
        for (std::tuple<TLorentzVector, bool> jet : quark)
        {
            if (std::get<0>(jet) == dummy_jet)
            {
                all_matched = false;
                break;
            }
        }
        um_vector.push_back(all_matched);
    }
    bool um_1 = false;
    bool um_2 = false;
    bool um_event = false;
    if (um_vector[0])
    {
        um_1 = true;
    }
    if (um_vector[1])
    {
        um_2 = true;
    }
    if (um_1 && um_2)
    {
        um_event = true;
    }
    return std::make_tuple(um_1, um_2, um_event);
}

/*// returns which top quarks are unambiguously matched
// NOTE: if top1 and top2 are individually unambiguously matched but there is >= 1 jet overlap, prefers top1 to be unambiguously matched
std::tuple<bool, bool, bool> unambiguously_matched(std::vector<std::vector<std::vector<std::tuple<TLorentzVector, bool>>>> all_matched_jets)
{
    std::vector<std::vector<std::tuple<TLorentzVector, bool>>> t1 = all_matched_jets[0];
    std::vector<std::vector<std::tuple<TLorentzVector, bool>>> t2 = all_matched_jets[1];
    bool unambiguously_matched_event = false;
    bool unambiguously_matched_1 = true;
    bool unambiguously_matched_2 = true;
    std::vector<std::tuple<TLorentzVector, bool>> used_jets;
    for (std::vector<std::tuple<TLorentzVector, bool>> matched_jets : t1)
    {
        if (unambiguously_matched_1)
        {
            if (matched_jets.size() != 1) // more than one matched jet for particle
            {
                unambiguously_matched_1 = false;
            }
            else if (matched_jets.size() == 1 && std::find(used_jets.begin(), used_jets.end(), matched_jets[0]) != used_jets.end()) // jet has already been used in top1
            {
                unambiguously_matched_1 = false;
            }
            else //only one matched jet that has not been used
            {
                used_jets.push_back(matched_jets[0]);
            }
        }
    }
    for (std::vector<std::tuple<TLorentzVector, bool>> matched_jets : t2)
    {
        if (unambiguously_matched_2)
        {
            if (matched_jets.size() != 1) // more than one matched jet for particle
            {
                unambiguously_matched_2 = false;
            }
            else if (matched_jets.size() == 1 && std::find(used_jets.begin(), used_jets.end(), matched_jets[0]) != used_jets.end()) // jet has already been used in top1 or top2
            {
                unambiguously_matched_2 = false;
            }
            else //only one matched jet that has not been used
            {
                used_jets.push_back(matched_jets[0]);
            }
        }
    }
    if (unambiguously_matched_1 && unambiguously_matched_2)
    {
        unambiguously_matched_event = true;
    }
    return std::make_tuple(unambiguously_matched_1, unambiguously_matched_2, unambiguously_matched_event);
}*/

// returns chi^2 value for ordered list of four-vectors
double chisq(std::vector<std::tuple<TLorentzVector, bool>> jets)
{   
        TLorentzVector jet1 = std::get<0>(jets[0]);
        TLorentzVector jet2 = std::get<0>(jets[1]);
        TLorentzVector jet3 = std::get<0>(jets[2]);
        TLorentzVector jet4 = std::get<0>(jets[3]);
        TLorentzVector jet5 = std::get<0>(jets[4]);
        TLorentzVector jet6 = std::get<0>(jets[5]);

        double recon_top_1 = (jet1 + jet2 + jet3).M();
        double recon_top_2 = (jet4 + jet5 + jet6).M();
        double sigma_W = 12.3;
        double sigma_diff_between_tops = 26.3;
        double recon_W_1 = (jet2 + jet3).M();
        double recon_W_2 = (jet4 + jet5).M();
        double known_W_mass = 81.3;

        double term_1 = (std::pow(recon_top_1 - recon_top_2, 2))/std::pow(sigma_diff_between_tops, 2);
        double term_2 = (std::pow(recon_W_1 - known_W_mass, 2))/std::pow(sigma_W, 2);
        double term_3 = (std::pow(recon_W_2 - known_W_mass, 2))/std::pow(sigma_W, 2);
        double chisq = term_1 + term_2 + term_3;
        return chisq;
}

// compares transverse momentum of two four-vectors
bool compare_pt(const std::tuple<TLorentzVector, bool>& a, const std::tuple<TLorentzVector, bool>& b) {
  return std::get<0>(a).Pt() < std::get<0>(b).Pt();
}

// returns vector of TLorentzVectors with mimimum chi^2 value
std::vector<TLorentzVector> min_chisq(std::vector<std::tuple<TLorentzVector, bool>> jets)
{
    double min_chisq_val = 9999999;
    std::vector<TLorentzVector> min_chisq_jets;
    do
    {
        if (std::get<1>(jets[0]) && std::get<1>(jets[3]))
        {
            double perm_chisq = chisq(jets);
            if (perm_chisq < min_chisq_val)
            {
                min_chisq_val = perm_chisq;
                min_chisq_jets.clear();
                for (std::tuple<TLorentzVector, bool> j : jets)
                {
                    min_chisq_jets.push_back(std::get<0>(j));
                }
            }
        }
    }
    while (std::next_permutation(jets.begin(), jets.end(), compare_pt));

    while (std::next_permutation(jets.begin(), jets.end(), compare_pt))
    {   
        double perm_chisq = chisq(jets);
        if (perm_chisq != -1 && perm_chisq < min_chisq_val)
        {
            min_chisq_val = perm_chisq;
            min_chisq_jets.clear();
            for (std::tuple<TLorentzVector, bool> j : jets)
            {
                min_chisq_jets.push_back(std::get<0>(j));
            }
        }
    }

    return min_chisq_jets;
}

// returns whether top quark is predited correctly by chisq method
bool match_eval(std::vector<std::tuple<TLorentzVector, bool>> ordered_jets_btagged, std::vector<TLorentzVector> chisq_jets)
{
    std::vector<TLorentzVector> ordered_jets;
    for (std::tuple<TLorentzVector, bool> jet : ordered_jets_btagged)
    {
        ordered_jets.push_back(std::get<0>(jet));
    }
    bool matched = true;
    for (int i = 0; i < 3; ++i)
    {
        if (ordered_jets[i] != chisq_jets[i])
        {
            matched = false;
            break;
        }
    }
    return matched;
}


/*
// returns vector of booleans indicating if the minimized chi^2 value matches the truth jets
std::vector<bool> match_eval(std::vector<std::tuple<TLorentzVector, bool>> matched_jets, std::vector<TLorentzVector> chisq_jets)
{
    std::vector<bool> match;
    if (std::get<0>(matched_jets[0]) == chisq_jets[0])
    {
        match.push_back(true);
    }
    else
    {
        match.push_back(false);
    }
    if (((std::get<0>(matched_jets[1]) == chisq_jets[1]) && (std::get<0>(matched_jets[2]) == chisq_jets[2])) ^ ((std::get<0>(matched_jets[2]) == chisq_jets[1]) && (std::get<0>(matched_jets[1]) == chisq_jets[2])))
    {
        match.push_back(true);
        match.push_back(true);
    }
    else
    {
        match.push_back(false);
        match.push_back(false);
    }
        if (std::get<0>(matched_jets[3]) == chisq_jets[3])
    {
        match.push_back(true);
    }
    else
    {
        match.push_back(false);
    }
    if (((std::get<0>(matched_jets[4]) == chisq_jets[4]) && (std::get<0>(matched_jets[5]) == chisq_jets[5])) ^ ((std::get<0>(matched_jets[5]) == chisq_jets[4]) && (std::get<0>(matched_jets[4]) == chisq_jets[5])))
    {
        match.push_back(true);
        match.push_back(true);
    }
    else
    {
        match.push_back(false);
        match.push_back(false);
    }
    return match;
}
*/

int main()
{
    int e_1_denom = 0;
    int e_2_denom = 0;
    int e_1_ct = 0;
    int e_2_ct = 0;
    int e_event_ct = 0;
    TFile* file = new TFile("top_recon_chisq_results.root", "RECREATE");
    TH1F* e1 = new TH1F("1_identifiable_top", "Masses of Reconstructed Tops for Events with 1 Identifiable Top Quark", 100, 0, 1000);
    e1 -> GetXaxis() -> SetTitle("Mass (GeV)");
    e1 -> GetYaxis() -> SetTitle("Frequency");
    TH1F* e2 = new TH1F("2_identifiable_tops", "Masses of Reconstructed Tops for Events with 2 Identifiable Top Quarks", 100, 0, 1000);
    e2 -> GetXaxis() -> SetTitle("Mass (GeV)");
    e2 -> GetYaxis() -> SetTitle("Frequency");
    TH1F* all = new TH1F("all_reconstructed_tops", "All Reconstructed Top Quarks", 100, 0, 1000);
    all -> GetXaxis() -> SetTitle("Mass (GeV)");
    all -> GetYaxis() -> SetTitle("Frequency");
    for (int event_num = 0; event_num < t -> GetEntries(); ++event_num)
    {
        std::cout << event_num << std::endl;
        t -> GetEntry(event_num);
        if (bJetNum -> GetLeaf("Jet_size") -> GetValue() < 6) // ignore events with < 6 jets
        {
            std::cout << "ERR: less than 6 jets" << std::endl;
            continue;
        }
        std::vector<std::vector<std::tuple<TLorentzVector, bool>>> ordered_particles = get_particles();
        std::tuple<int, int, std::vector<std::tuple<TLorentzVector, bool>>> get_jet_output = get_jets();
        auto [num_jets, num_btagged, unordered_jets] = get_jet_output;
        std::cout << "num jets " << num_jets << std::endl;
        std::cout << "num btagged " << num_btagged << std::endl;
        if (num_btagged < 2) // ignore events with < 2 btagged jets
        {
            std::cout << "ERR: less than 2 btagged jets" << std::endl;
            continue;
        }
        std::vector<std::vector<std::tuple<TLorentzVector, bool>>> ordered_jets = match(ordered_particles, unordered_jets);
        std::cout << "matched" << std::endl;
        std::tuple<bool, bool, bool> um_output = unambiguously_matched(ordered_jets);
        std::cout << "um" << std::endl;
        auto[unambiguously_matched_1, unambiguously_matched_2, unambiguously_matched_event] = um_output;
        std::vector<TLorentzVector> chisq_jets_unformatted = min_chisq(unordered_jets);
        std::cout << "min chisq" << std::endl;
        std::cout << chisq_jets_unformatted.size() << std::endl;
        for (TLorentzVector jet : chisq_jets_unformatted)
        {
            std::cout << jet.M() << " " << jet.Pt() << std::endl;
        }
        std::vector<std::vector<TLorentzVector>> chisq_jets = {{chisq_jets_unformatted[0], chisq_jets_unformatted[1], chisq_jets_unformatted[2]}, {chisq_jets_unformatted[3], chisq_jets_unformatted[4], chisq_jets_unformatted[5]}};
        std::cout << "chisq jets created" << std::endl;
        if (unambiguously_matched_event)
        {
            std::cout << "um event" << std::endl;
            e_2_denom ++;
            int num_correctly_predicted_tops = 0;
            if (match_eval(ordered_jets[0], chisq_jets[0])) // first top is correctly predicted
            {
                num_correctly_predicted_tops++;
            }
            if (match_eval(ordered_jets[1], chisq_jets[1])) // second top is correctly predicted
            {
                num_correctly_predicted_tops++;
            }
            if (num_correctly_predicted_tops == 1)
            {
                e_2_ct ++;
            } 
            else if (num_correctly_predicted_tops == 2)
            {
                e_event_ct ++;
            }
        }
        else if (unambiguously_matched_1 ^ unambiguously_matched_2) // only one top matched
        {
            std::cout << "um 1 or 2" << std::endl;
            e_1_denom ++;
            if (unambiguously_matched_1)
            {
                if (match_eval(ordered_jets[0], chisq_jets[0]))
                {
                    e_1_ct ++;
                }
            }
            else if (unambiguously_matched_2)
            {
                if (match_eval(ordered_jets[1], chisq_jets[1]))
                {
                    e_1_ct ++;
                }
            }
        }


        /*
        if (unambiguously_matched_event) // both top quarks are correctly matched
        {
            e_2_denom ++;
            std::vector<TLorentzVector> chisq_jets = min_chisq(unordered_jets);
            std::vector<std::tuple<TLorentzVector, bool>> ordered_jets;
            for (std::vector<std::vector<std::tuple<TLorentzVector, bool>>> top : all_matched_jets)
            {
                for (std::vector<std::tuple<TLorentzVector, bool>> j : top)
                {
                    ordered_jets.push_back(j[0]);
                }
            }
            std::vector<bool> matched = match_eval(ordered_jets, chisq_jets);
            int matched_tops = 0;
            if (matched[0] && matched[1] && matched[2])
            {
                matched_tops ++;
                TLorentzVector recon_top = (chisq_jets[0] + chisq_jets[1] + chisq_jets[2]);
                e2 -> Fill(recon_top.M());
                all -> Fill(recon_top.M());
            }
            if (matched[3] && matched[4] && matched[5])
            {
                matched_tops ++;
                TLorentzVector recon_top = (chisq_jets[3] + chisq_jets[4] + chisq_jets[5]);
                e2 -> Fill(recon_top.M());
                all -> Fill(recon_top.M());
            }
            if (matched_tops == 1)
            {
                e_2_ct ++;
            }
            else if (matched_tops == 2)
            {
                e_event_ct ++;
            }
        }
        else if (unambiguously_matched_1 ^ unambiguously_matched_2) // only one top quark is correctly matched
        {
            e_1_denom ++;
            std::vector<TLorentzVector> chisq_jets = min_chisq(unordered_jets);
            std::vector<std::tuple<TLorentzVector, bool>> ordered_jets;
            if (unambiguously_matched_1) // first top quark is correctly matched
            {
                for (std::vector<std::tuple<TLorentzVector, bool>> j : all_matched_jets[0])
                {
                    ordered_jets.push_back(j[0]);
                }
                for (int i = 0; i < 3; ++i)
                {
                    TLorentzVector v;
                    v.SetPtEtaPhiM(0, 0, 0, 0);
                    if (i == 0)
                    {
                        ordered_jets.push_back(std::make_tuple(v,true));
                    }
                    else
                    {
                        ordered_jets.push_back(std::make_tuple(v,false));
                    }
                }
                std::vector<bool> matched = match_eval(ordered_jets, chisq_jets);
                int matched_tops = 0;
                if (matched[0] && matched[1] && matched[2])
                {
                    TLorentzVector recon_top = (chisq_jets[0] + chisq_jets[1] + chisq_jets[2]);
                    e1 -> Fill(recon_top.M());
                    all -> Fill(recon_top.M());
                    matched_tops ++;
                }
                if (matched_tops == 1)
                {
                    e_1_ct ++;
                }
            }
            else if (unambiguously_matched_2) // second top quark is correctly matched
            {
                for (int i = 0; i < 3; ++i)
                {
                    TLorentzVector v;
                    v.SetPtEtaPhiM(0, 0, 0, 0);
                    if (i == 0)
                    {
                        ordered_jets.push_back(std::make_tuple(v,true));
                    }
                    else
                    {
                        ordered_jets.push_back(std::make_tuple(v,false));
                    }
                }
                for (std::vector<std::tuple<TLorentzVector, bool>> j : all_matched_jets[1])
                {
                    ordered_jets.push_back(j[0]);
                }
                std::vector<bool> matched = match_eval(ordered_jets, chisq_jets);
                int matched_tops = 0;
                if (matched[3] && matched[4] && matched[5])
                {
                    TLorentzVector recon_top = (chisq_jets[3] + chisq_jets[4] + chisq_jets[5]);
                    e1 -> Fill(recon_top.M());
                    all -> Fill(recon_top.M());
                    matched_tops ++;
                }
                if (matched_tops == 1)
                {
                    e_1_ct ++;
                }
            }
        }*/
        std::cout << "EVENT " << event_num <<std::endl;
        std::cout << "e_1 Denominator: " << e_1_denom << std::endl;
        std::cout << "e_2 Denominator: " << e_2_denom << std::endl;
        if (e_1_denom != 0)
        {
            std::cout << "e_top_1 Efficiency: " << e_1_ct/e_1_denom << std::endl;
        }
        if (e_2_denom != 0)
        {
            std::cout << "e_top_2 Efficiency: " << e_2_ct/e_2_denom << std::endl;
            std::cout << "e_event Efficiency: " << e_event_ct/e_2_denom << std::endl;
        }
    }
    e1 -> Write();
    e2 -> Write();
    all -> Write();
    file -> Close();

    std::cout << "e_1 Denominator: " << e_1_denom << std::endl;
    std::cout << "e_2 Denominator: " << e_2_denom << std::endl;
    std::cout << "e_top_1 Efficiency: " << e_1_ct/e_1_denom << std::endl;
    std::cout << "e_top_2 Efficiency: " << e_2_ct/e_2_denom << std::endl;
    std::cout << "e_event Efficiency: " << e_event_ct/e_2_denom << std::endl;
    return 0;
}