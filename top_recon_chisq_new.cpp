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
#include "Tleaf.h"
#include "TLeaf.h"
#include "TH1F.h"

using namespace std;

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

TLeaf *JetBTag = t -> GetBranch("Jet.BTag") -> GetLeaf("Jet.BTag");
TLeaf *JetPT = t -> GetBranch("et.PT") -> GetLeaf("Jet.PT");
TLeaf *JetEta = t -> GetBranch("Jet.Eta") -> GetLeaf("Jet.Eta");
TLeaf *JetPhi = t -> GetBranch("Jet.Phi") -> GetLeaf("Jet.Phi");
TLeaf *JetM = t -> GetBranch("Jet.Mass") -> GetLeaf("Jet.Mass");
TLeaf *JetNum = t -> GetBranch("Jet_size") -> GetLeaf("Jet_size");
TLeaf *ParticlePT = t -> GetBranch("Particle.PT") -> GetLeaf("Particle.PT");
TLeaf *ParticleEta = t -> GetBranch("Particle.Eta") -> GetLeaf("Particle.Eta");
TLeaf *ParticlePhi = t -> GetBranch("Particle.Phi") -> GetLeaf("Particle.Phi");
TLeaf *ParticleM = t -> GetBranch("Particle.Mass") -> GetLeaf("Particle.Mass");
TLeaf *ParticleID = t -> GetBranch("Particle.PID") -> GetLeaf("Particle.PID");
TLeaf *ParticleStatus = t -> GetBranch("Particle.Status") -> GetLeaf("Particle.Status");
TLeaf *ParticleMother = t -> GetBranch("Particle.M1") -> GetLeaf("Particle.M1");

// returns an ordered list of particles formatted as [[b,q,q~],[b~,q',q'~]]
// note that each particle is a tuple of its four-vector and btag
vector<vector<tuple<TLorentzVector, bool>>> get_particles()
{
    vector<int> unordered_particles; // holds all final state particle indexes
    for (int p = 0; p < t -> GetLeaf("Particle_size") -> GetValue(); ++p)
    {
        if (ParticleStatus -> GetValue(p) == 23) // is final state particle
        {
            unordered_particles.push_back(p);
        }
    }
    vector<tuple<TLorentzVector, bool>> t1;
    vector<tuple<TLorentzVector, bool>> t2;
    for (int q : unordered_particles) // add bottom quark to t1
    {
        if (ParticleID -> GetValue(q) == 5) // is bottom quark
        {
            TLorentzVector particle;
            particle.SetPtEtaPhiM(ParticlePT -> GetValue(q),
                                  ParticleEta -> GetValue(q),
                                  ParticlePhi -> GetValue(q),
                                  ParticleM -> GetValue(q));
            tuple<TLorentzVector, bool> particle_tuple = {particle, true};
            t1.push_back(particle_tuple);
        }
    }
    for (int q : unordered_particles) // add q and q~ to t1
    {
        if (ParticleID -> GetValue(ParticleMother -> GetValue(q)) == 24) // originated from W+ boson
        {
            TLorentzVector particle;
            particle.SetPtEtaPhiM(ParticlePT -> GetValue(q),
                                  ParticleEta -> GetValue(q),
                                  ParticlePhi -> GetValue(q),
                                  ParticleM -> GetValue(q));
            tuple<TLorentzVector, bool> particle_tuple = {particle, false};
            t1.push_back(particle_tuple);
        }
    }
    for (int q : unordered_particles)
    {
        if (ParticleID ->  GetValue(q) == -5) // is anti-bottom quark
        {
            TLorentzVector particle;
            particle.SetPtEtaPhiM(ParticlePT -> GetValue(q),
                                  ParticleEta -> GetValue(q),
                                  ParticlePhi -> GetValue(q),
                                  ParticleM -> GetValue(q));
            tuple<TLorentzVector, bool> particle_tuple = {particle, true};
            t2.push_back(particle_tuple);
        }
    }
    for (int q : unordered_particles)
    {
        if (ParticleID -> GetValue(ParticleMother -> GetValue(q)) == -24) // originated from W- boson
        {
            TLorentzVector particle;
            particle.SetPtEtaPhiM(ParticlePT -> GetValue(q),
                                  ParticleEta -> GetValue(q),
                                  ParticlePhi -> GetValue(q),
                                  ParticleM -> GetValue(q));
            tuple<TLorentzVector, bool> particle_tuple = {particle, false};
            t2.push_back(particle_tuple);
        }
    }
    vector<vector<tuple<TLorentzVector, bool>>> ordered_particles = {t1,t2};
    for (vector<tuple<TLorentzVector, bool>> t : ordered_particles)
    {
        cout << "TOP QUARK" << endl;
        for (tuple<TLorentzVector, bool> p : t)
        {
            cout << "M: " << get<0>(p).M() << " PT: " << get<0>(p).Pt() << " Btag: " << get<1>(p) << endl;
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



int main()
{
    for (int i = 0; i< t -> GetEntries(); ++i)
    {
        t -> GetEntry(i);
        vector<vector<tuple<TLorentzVector, bool>>> particles = get_particles();
    }
    return 0;
}