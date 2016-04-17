#Build FCR (and BFA)
root -l <<< ".L buildtupledata.C+"
root -l <<< ".L buildtuplemc.C+"
#pp
#root -l -q -b buildtuplemc.C+\(\"mcppbfcak4PF\"\) >out_mcppbfcak4PF &
#PbPb
root -l -q -b buildtuplemc.C+\(\"mcPbbfcakPu4PF\"\) >out_mcPbbfcakPu4PF &
root -l -q -b buildtuplemc.C+\(\"mcPbbfcakPu3PF\"\) >out_mcPbbfcakPu3PF &
#root -l -q -b buildtuplemc.C+\(\"mcPbbfcakCs4PF\"\) >out_mcPbbfcakCs4PF &
#root -l -q -b buildtuplemc.C+\(\"mcPbbfcakCs3PF\"\) >out_mcPbbfcakCs3PF &

#root -l -q -b buildtuplemc.C+\(\"mcPbbfcakPu4Calo\"\) >out_mcPbbfcakPu4Calo &
#root -l -q -b buildtuplemc.C+\(\"mcPbbfcakPu3Calo\"\) >out_mcPbbfcakPu3Calo &
#root -l -q -b buildtuplemc.C+\(\"mcPbbfcakCs4Calo\"\) >out_mcPbbfcakCs4Calo &
#root -l -q -b buildtuplemc.C+\(\"mcPbbfcakCs3Calo\"\) >out_mcPbbfcakCs3Calo &