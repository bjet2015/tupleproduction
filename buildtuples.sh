#root -l -q buildtuplemc.C\(\"PbPb\",\"qcd\",\"akPu4PFJetAnalyzer\"\) &
root -l -q buildtuplemc.C\(\"pp\",\"qcd\",\"ak4PFJetAnalyzer\"\) &



#root -l -q -b buildtupledata.C\(\"pp\",\"ak4PFJetAnalyzer\"\) &

#root -l -q -b buildtupledata.C\(\"PbPb\",\"akPu4PFJetAnalyzer\"\) &
#root -l -q -b buildtupledata.C\(\"PbPbBJet\",\"akPu4PFJetAnalyzer\"\) &