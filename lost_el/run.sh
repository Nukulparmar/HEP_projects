#!/bin/bash

echo "Runing on 2016"

./lost_el /home/work/nparmar/runlist/runlist_TTGJets.txt /home/work/nparmar/analysis/lost_el/output_files/TTGJets_lost_el.root TTGJets 2016

./lost_el /home/work/nparmar/runlist/runlist_WGJets.txt /home/work/nparmar/analysis/lost_el/output_files/WGJets_lost_el.root WGJets 2016

./lost_el /home/work/nparmar/runlist/runlist_WJets.txt /home/work/nparmar/analysis/lost_el/output_files/WJets_lost_el.root WJets 2016

./lost_el /home/work/nparmar/runlist/runlist_TTJets.txt /home/work/nparmar/analysis/lost_el/output_files/TTJets_lost_el.root TTJets 2016

echo "hadd command running for 2016"
hadd -f /home/work/nparmar/analysis/lost_el/output_files/TTWGamma_lost_el_2016.root /home/work/nparmar/analysis/lost_el/output_files/TTGJets_lost_el.root /home/work/nparmar/analysis/lost_el/output_files/WGJets_lost_el.root
hadd -f /home/work/nparmar/analysis/lost_el/output_files/All_TTW_lost_el_2016.root /home/work/nparmar/analysis/lost_el/output_files/TTGJets_lost_el.root /home/work/nparmar/analysis/lost_el/output_files/WGJets_lost_el.root /home/work/nparmar/analysis/lost_el/output_files/WJets_lost_el.root /home/work/nparmar/analysis/lost_el/output_files/TTJets_lost_el.root
# ./lost_el /home/work/nparmar/runlist/runlist_both.txt /home/work/nparmar/analysis/lost_el/output_files/both_lost_el.root 2016

# ./lost_el /home/work/nparmar/runlist/runlist_all.txt /home/work/nparmar/analysis/lost_el/output_files/all_lost_el.root 2016

echo "Runing on 2017"
./lost_el /home/work/nparmar/runlist/runlist_ttgamma_2017.txt /home/work/nparmar/analysis/lost_el/output_files/TTGJets_lost_el_2017.root TTGJets 2017

./lost_el /home/work/nparmar/runlist/runlist_WGJets_2017.txt /home/work/nparmar/analysis/lost_el/output_files/WGJets_lost_el_2017.root WGJets 2017

./lost_el /home/work/nparmar/runlist/runlist_WJets_2017.txt /home/work/nparmar/analysis/lost_el/output_files/WJets_lost_el_2017.root WJets 2017

./lost_el /home/work/nparmar/runlist/runlist_TTJets_2017.txt /home/work/nparmar/analysis/lost_el/output_files/TTJets_lost_el_2017.root TTJets 2017

echo "hadd command running for 2017"
hadd -f /home/work/nparmar/analysis/lost_el/output_files/TTWGamma_lost_el_2017.root /home/work/nparmar/analysis/lost_el/output_files/TTGJets_lost_el_2017.root /home/work/nparmar/analysis/lost_el/output_files/WGJets_lost_el_2017.root
hadd -f /home/work/nparmar/analysis/lost_el/output_files/All_TTW_lost_el_2017.root /home/work/nparmar/analysis/lost_el/output_files/TTGJets_lost_el_2017.root /home/work/nparmar/analysis/lost_el/output_files/WGJets_lost_el_2017.root /home/work/nparmar/analysis/lost_el/output_files/WJets_lost_el_2017.root /home/work/nparmar/analysis/lost_el/output_files/TTJets_lost_el_2017.root

# ./lost_el /home/work/nparmar/runlist/runlist_both_2017.txt /home/work/nparmar/analysis/lost_el/output_files/both_lost_el_2017.root 2017

# ./lost_el /home/work/nparmar/runlist/runlist_all_2017.txt /home/work/nparmar/analysis/lost_el/output_files/all_lost_el_2017.root 2017

echo "Runing on 2018"
./lost_el /home/work/nparmar/runlist/runlist_ttgjets_2018.txt /home/work/nparmar/analysis/lost_el/output_files/TTGJets_lost_el_2018.root TTGJets 2018

./lost_el /home/work/nparmar/runlist/runlist_WJets_2018.txt /home/work/nparmar/analysis/lost_el/output_files/WJets_lost_el_2018.root WJets 2018

./lost_el /home/work/nparmar/runlist/runlist_WGJets_2018.txt /home/work/nparmar/analysis/lost_el/output_files/WGJets_lost_el_2018.root WGJets 2018

./lost_el /home/work/nparmar/runlist/runlist_TTJets_2018.txt /home/work/nparmar/analysis/lost_el/output_files/TTJets_lost_el_2018.root TTJets 2018

echo "hadd command running for 2018"
hadd -f /home/work/nparmar/analysis/lost_el/output_files/TTWGamma_lost_el_2018.root /home/work/nparmar/analysis/lost_el/output_files/TTGJets_lost_el_2018.root 
hadd -f /home/work/nparmar/analysis/lost_el/output_files/All_TTW_lost_el_2018.root /home/work/nparmar/analysis/lost_el/output_files/TTGJets_lost_el_2018.root /home/work/nparmar/analysis/lost_el/output_files/WJets_lost_el_2018.root /home/work/nparmar/analysis/lost_el/output_files/TTJets_lost_el_2018.root /home/work/nparmar/analysis/lost_el/output_files/WGJets_lost_el_2018.root


# root /home/work/nparmar/analysis/lost_el/output_files/TTGJets_lost_el.root /home/work/nparmar/analysis/lost_el/output_files/WGJets_lost_el.root /home/work/nparmar/analysis/lost_el/output_files/both_lost_el.root

read -p "Enter the name of the new folder:" folder
mkdir /home/work/nparmar/analysis/lost_el/$folder
read -p "Enter some details about the current run:" details
echo "$details" > /home/work/nparmar/analysis/lost_el/$folder/readme.txt

# root -l -q 'TF_lost_el_all_region.C("both",2016,"/home/work/nparmar/analysis/lost_el/'$folder'/lost_electron_2016.pdf ","/home/work/nparmar/analysis/lost_el/'$folder'/Closure_2016.pdf)' > /home/work/nparmar/analysis/lost_el/$folder/tf_lost_el_2016.txt
# root -l -q 'TF_lost_el_all_region.C("both",2017,"/home/work/nparmar/analysis/lost_el/'$folder'/lost_electron_2017.pdf","/home/work/nparmar/analysis/lost_el/'$folder'/Closure_2017.pdf")' > /home/work/nparmar/analysis/lost_el/$folder/tf_lost_el_2017.txt
# root -l -q 'TF_lost_el_all_region.C("both",2018,"/home/work/nparmar/analysis/lost_el/'$folder'/lost_electron_2018.pdf","/home/work/nparmar/analysis/lost_el/'$folder'/Closure_2018.pdf")' > /home/work/nparmar/analysis/lost_el/$folder/tf_lost_el_2018.txt
# root -l -q 'TF_lost_el_all_region.C("ttgjets",2016,"/home/work/nparmar/analysis/lost_el/'$folder'/ttgjets_lost_el_2016.pdf")' > /home/work/nparmar/analysis/lost_el/$folder/tf_ttgjets_2016.txt
# root -l -q 'TF_lost_el_all_region.C("ttgjets",2017,"/home/work/nparmar/analysis/lost_el/'$folder'/ttgjets_lost_el_2017.pdf")' > /home/work/nparmar/analysis/lost_el/$folder/tf_ttgjets_2017.txt
# root -l -q 'TF_lost_el_all_region.C("ttgjets",2018,"/home/work/nparmar/analysis/lost_el/'$folder'/ttgjets_lost_el_2018.pdf")' > /home/work/nparmar/analysis/lost_el/$folder/tf_ttgjets_2018.txt
# root -l -q 'TF_lost_el_all_region.C("wgjets",2016,"/home/work/nparmar/analysis/lost_el/'$folder'/wgjets_lost_el_2016.pdf")' > /home/work/nparmar/analysis/lost_el/$folder/tf_wgjets_2016.txt 
# root -l -q 'TF_lost_el_all_region.C("wgjets",2017,"/home/work/nparmar/analysis/lost_el/'$folder'/wgjets_lost_el_2017.pdf")' > /home/work/nparmar/analysis/lost_el/$folder/tf_wgjets_2017.txt
root -l -q 'TF_lost_el_all_region.C("all",2016,"/home/work/nparmar/analysis/lost_el/'$folder'/lost_electron_2016_all.pdf","/home/work/nparmar/analysis/lost_el/'$folder'/Closure_2016_all.pdf")' > /home/work/nparmar/analysis/lost_el/$folder/tf_lost_el_all_2016.txt
root -l -q 'TF_lost_el_all_region.C("all",2017,"/home/work/nparmar/analysis/lost_el/'$folder'/lost_electron_2017_all.pdf","/home/work/nparmar/analysis/lost_el/'$folder'/Closure_2017_all.pdf")' > /home/work/nparmar/analysis/lost_el/$folder/tf_lost_el_all_2017.txt
root -l -q 'TF_lost_el_all_region.C("all",2018,"/home/work/nparmar/analysis/lost_el/'$folder'/lost_electron_2018_all.pdf","/home/work/nparmar/analysis/lost_el/'$folder'/Closure_2018_all.pdf")' > /home/work/nparmar/analysis/lost_el/$folder/tf_lost_el_all_2018.txt
