import uproot
import awkward
from awkward import Array
import uproot_methods.classes.TLorentzVector as LV
import numpy as np

# DECAY PROCESS
# p p -> t t~
# t -> b W+
# W+ -> q q~
# t~ -> b~ W-
# W- -> q q~

# EFFICIENCY DEFINITIONS
# 3 efficiencies: eff_1, eff_2, eff_event
# 3 numerators (e_1, e_2, e_event), 2 denominators (e_1, e_2)
# Denominators:
#  subscript corresponds to number of matched (matched to truth particle) tops
#  e_1 means only one top quark was matched, e_2 means both top quarks were matched
# Numerators:
#  subscript corresponds to number of top quarks reconstructed correctly using chisq method
#  e_1 means one top quark was reconstructed correctly (for an e_1 event), e_2 means one top quark was reconstructed correctly (for an e_2 event), e_event means both top quarks were reconstructed correctly (for an e_2 event)

file = uproot.open("/Users/divitrawal/MG5_aMC_v3_4_0/TOPRECON/Events/run_01/tag_1_delphes_events.root")
tree = file["Delphes"]

bJetBTag = tree["Jet.BTag"].array()
print(bJetBTag)