#!/bin/bash

#./had_tau /home/nukul/EHEP_Projects/Summer_2019/GJets/runlist_gjets.txt /home/nukul/EHEP_Projects/Summer_2019/GJets/gjets_had_tau.root a

#./had_tau /home/nukul/EHEP_Projects/Summer_2019/QCD/runlist_qcd.txt /home/nukul/EHEP_Projects/Summer_2019/QCD/qcd_had_tau.root  a

./had_tau /home/nukul/EHEP_Projects/Summer_2019/TTGJets/runlist_ttgjets.txt /home/nukul/EHEP_Projects/Summer_2019/TTGJets/ttgjets_had_tau.root  a

./had_tau /home/nukul/EHEP_Projects/Summer_2019/TTJets/runlist_ttjets.txt /home/nukul/EHEP_Projects/Summer_2019/TTJets/ttjets_had_tau.root  a

./had_tau /home/nukul/EHEP_Projects/Summer_2019/WGJets_monophotons/runlist_wgjets_monophoton.txt /home/nukul/EHEP_Projects/Summer_2019/WGJets_monophotons/wgjets_monophoton_had_tau.root  a

./had_tau /home/nukul/EHEP_Projects/Summer_2019/WGJets_to_Lnu/runlist_wgjets_lnu.txt /home/nukul/EHEP_Projects/Summer_2019/WGJets_to_Lnu/wgjets_lnu_had_tau.root a

#./had_tau /home/nukul/EHEP_Projects/Summer_2019/ZJets_to_nunu/runlist_zjets_nunu.txt /home/nukul/EHEP_Projects/Summer_2019/ZJets_to_nunu/zjets_nunu_had_tau.root  a

./had_tau /home/nukul/EHEP_Projects/Summer_2019/WGJets/runlist_wgjets.txt /home/nukul/EHEP_Projects/Summer_2019/WGJets/wgjets_had_tau.root a


root /home/nukul/EHEP_Projects/Summer_2019/TTGJets/ttgjets_had_tau.root /home/nukul/EHEP_Projects/Summer_2019/TTJets/ttjets_had_tau.root /home/nukul/EHEP_Projects/Summer_2019/WGJets_monophotons/wgjets_monophoton_had_tau.root /home/nukul/EHEP_Projects/Summer_2019/WGJets_to_Lnu/wgjets_lnu_had_tau.root /home/nukul/EHEP_Projects/Summer_2019/WGJets/wgjets_had_tau.root


#/home/nukul/EHEP_Projects/Summer_2019/GJets/gjets_had_tau.root
#/home/nukul/EHEP_Projects/Summer_2019/ZJets_to_nunu/zjets_nunu_had_tau.root
# /home/nukul/EHEP_Projects/Summer_2019/QCD/qcd_had_tau.root 
