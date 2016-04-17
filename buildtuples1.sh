root -l <<< ".L buildtupledata.C+"
root -l <<< ".L buildtuplemc.C+"

#Build all data tuples
root -l -q -b buildtupledata.C+\(\"dtppjpfak4PF\"\) >out_dtppjpfak4PF &
root -l -q -b buildtupledata.C+\(\"dtppjpfak3PF\"\) >out_dtppjpfak3PF &

root -l -q -b buildtupledata.C+\(\"dtPbj60akPu4PF\"\) >out_dtPbj60akPu4PF &
root -l -q -b buildtupledata.C+\(\"dtPbbjtakPu4PF\"\) >out_dtPbbjtakPu4PF &
root -l -q -b buildtupledata.C+\(\"dtPbj60akPu3PF\"\) >out_dtPbj60akPu3PF &
root -l -q -b buildtupledata.C+\(\"dtPbbjtakPu3PF\"\) >out_dtPbbjtakPu3PF &

root -l -q -b buildtupledata.C+\(\"dtPbj60akCs4PF\"\) >out_dtPbj60akCs4PF &
root -l -q -b buildtupledata.C+\(\"dtPbbjtakCs4PF\"\) >out_dtPbbjtakCs4PF &
root -l -q -b buildtupledata.C+\(\"dtPbj60akCs3PF\"\) >out_dtPbj60akCs3PF &
root -l -q -b buildtupledata.C+\(\"dtPbbjtakCs3PF\"\) >out_dtPbbjtakCs3PF &
