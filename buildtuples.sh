#root -l -q buildtuplemc.C\(\"PbPb\",\"qcd\"\) &
root -l -q buildtuplemc.C\(\"pp\",\"qcd\"\) &
#root -l -q buildtuplemc.C\(\"PbPb\",\"bjet\"\) &
#root -l -q buildtuplemc.C\(\"pp\",\"bjet\"\) &

root -l -q buildtupledata.C\(\"pp_PFLowPt\",\"ak4PFJetAnalyzer\"\) &
