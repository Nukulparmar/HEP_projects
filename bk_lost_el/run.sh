echo "TTGJets"
./lost_el /home/work/nparmar/runlist/runlist_TTGJets.txt /home/work/nparmar/git_workspace/HEP_projects/output_files/TTGJets_lost_el.root a
echo "TTJets"
./lost_el /home/work/nparmar/runlist/runlist_ttjets.txt /home/work/nparmar/git_workspace/HEP_projects/output_files/ttjets_lost_el.root a
echo "WGJets"
./lost_el /home/work/nparmar/runlist/runlist_WGJets.txt /home/work/nparmar/git_workspace/HEP_projects/output_files/WGJets_lost_el.root a
echo "WJets"
./lost_el /home/work/nparmar/runlist/runlist_WJets.txt /home/work/nparmar/git_workspace/HEP_projects/output_files/WJets_lost_el.root a

root /home/work/nparmar/git_workspace/HEP_projects/output_files/*.root

