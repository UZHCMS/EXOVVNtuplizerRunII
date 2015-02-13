from ROOT import *

numberOfEvents=1
options="precision=3,colsize=7"

f1=TFile("../../dijetTree_test.root")
tree1=f1.Get("dijets").Get("events")
branches=[b.GetName() for b in tree1.GetListOfBranches()]
print branches
branches_selected=[
#'runNo', 'lumi',
'evtNo', 
'rho', 
#'nJetsAK8',
'jetPtAK8', 'jetEtaAK8', 'jetPhiAK8', 'jetMassAK8', 'jetEnergyAK8',
'jetJecAK8', 'jetAreaAK8',
'idLAK8',# 'idTAK8',
'jetChfAK8', 'jetNhfAK8', 'jetPhfAK8', 'jetMufAK8', 'jetElfAK8',
'jetMassPrunedAK8', 'jetMassSoftDropAK8',
'jetTau1AK8', 'jetTau2AK8', 'jetTau3AK8']
tree1.GetPlayer().SetScanRedirect(True)
tree1.GetPlayer().SetScanFileName("tree1.txt")
tree1.Scan(":".join(branches_selected),"",options,numberOfEvents)

f2=TFile("flatTuple.root")
tree2=f2.Get("ntuplizer").Get("tree")
branches=[b.GetName() for b in tree2.GetListOfBranches()]
print branches
branches_selected=[
#'EVENT_run', 'EVENT_lumiBlock',
'EVENT_event',
'rho',
#'njetsAK8',
'jetAK8_pt', 'jetAK8_eta', 'jetAK8_phi', 'jetAK8_mass', 'jetAK8_e',
'jetAK8_jec', 'jetAK8_area',
'jetAK8_IDLoose',
'jetAK8_chf', 'jetAK8_nhf', 'jetAK8_phf', 'jetAK8_muf', 'jetAK8_emf',
'jetAK8_prunedmassUnCorr', 'jetAK8_softdropmassCorr',
'jetAK8_tau1', 'jetAK8_tau2', 'jetAK8_tau3', #'jetAK8_tau21', 
]
tree2.GetPlayer().SetScanRedirect(True)
tree2.GetPlayer().SetScanFileName("tree2.txt")
tree2.Scan(":".join(branches_selected),"",options,numberOfEvents)

t1=open("tree1.txt")
t2=open("tree2.txt")

l1=t1.readlines()
l2=t2.readlines()

i1=0
i2=0
while i1<len(l1) and i2<len(l2):
  index1=l1[i1].split("*")[2].strip()
  index2=l2[i2].split("*")[2].strip()
  try:
    if float(index1)>float(index2):
     print "tree2 has less entries than tree1:\n",l1[i1]
     i1+=1
     continue
    if float(index1)<float(index2):
     print "tree1 has less entries than tree2:\n", l2[i2]
     i2+=1
     continue
  except: pass
  if l1[i1]!=l2[i2]:
   print "tree1 and tree2 differ:\n",l1[i1],"\n",l2[i2]
   pass
  i1+=1
  i2+=1
