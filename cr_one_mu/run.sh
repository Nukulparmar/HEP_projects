#!/bin/bash

./one_mu /home/nukul/EHEP_Projects/Summer_2019/GJets/runlist_gjets.txt /home/nukul/EHEP_Projects/Summer_2019/GJets/gjets_one_mu.root a

./one_mu /home/nukul/EHEP_Projects/Summer_2019/QCD/runlist_qcd.txt /home/nukul/EHEP_Projects/Summer_2019/QCD/qcd_one_mu.root  a

./one_mu /home/nukul/EHEP_Projects/Summer_2019/TTGJets/runlist_ttgjets.txt /home/nukul/EHEP_Projects/Summer_2019/TTGJets/ttgjets_one_mu.root  a

./one_mu /home/nukul/EHEP_Projects/Summer_2019/TTJets/runlist_ttjets.txt /home/nukul/EHEP_Projects/Summer_2019/TTJets/ttjets_one_mu.root  a

./one_mu /home/nukul/EHEP_Projects/Summer_2019/WGJets_monophotons/runlist_wgjets_monophoton.txt /home/nukul/EHEP_Projects/Summer_2019/WGJets_monophotons/wgjets_monophoton_one_mu.root  a

./one_mu /home/nukul/EHEP_Projects/Summer_2019/WGJets_to_Lnu/runlist_wgjets_lnu.txt /home/nukul/EHEP_Projects/Summer_2019/WGJets_to_Lnu/wgjets_lnu_one_mu.root a

./one_mu /home/nukul/EHEP_Projects/Summer_2019/ZJets_to_nunu/runlist_zjets_nunu.txt /home/nukul/EHEP_Projects/Summer_2019/ZJets_to_nunu/zjets_nunu_one_mu.root  a


root /home/nukul/EHEP_Projects/Summer_2019/GJets/gjets_one_mu.root /home/nukul/EHEP_Projects/Summer_2019/QCD/qcd_one_mu.root /home/nukul/EHEP_Projects/Summer_2019/TTGJets/ttgjets_one_mu.root /home/nukul/EHEP_Projects/Summer_2019/TTJets/ttjets_one_mu.root /home/nukul/EHEP_Projects/Summer_2019/WGJets_monophotons/wgjets_monophoton_one_mu.root /home/nukul/EHEP_Projects/Summer_2019/WGJets_to_Lnu/wgjets_lnu_one_mu.root /home/nukul/EHEP_Projects/Summer_2019/ZJets_to_nunu/zjets_nunu_one_mu.root
