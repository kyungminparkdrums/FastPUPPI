import os

#indir = '/eos/cms/store/group/cmst3/group/l1tr/kypark/REPO_L1HGcalID/OldID/'
indir = '/eos/cms/store/cmst3/group/l1tr/kypark/BestMET/AR_Prep_April24_tuneWP_more/'

f_ttbar = f'{indir}/TT_PU200/FP/perfNano.root'
f_nugun = f'{indir}/NuGunAllEta_PU200/FP/perfNano.root'

outdir = '/eos/user/k/kypark/www/L1T_HGC_TunePFwp/NewID/'

cmd_jecs = f'python3 scripts/makeJecs.py {f_ttbar} -A -o {outdir}/jecs.root'
#cmd_puppi = f'python3 scripts/jetHtSuite.py {f_ttbar} {f_nugun} {outdir}/METcentral -j {outdir}/jecs.root -w l1pfpu_metnoref -v metCentral --eta 2.4'
#cmd_puppi = f'python3 scripts/jetHtSuite.py {f_ttbar} {f_nugun} {outdir}/MET -j {outdir}/jecs.root -w l1pfpu_metnoref -v met --eta 5.0'
#cmd_puppi = f'python3 scripts/jetHtSuite.py {f_ttbar} {f_nugun} {outdir}/JET -j {outdir}/jecs.root -w l1pfpu_jetnoref -v jet4 --eta 5.0'
cmd_puppi = f'python3 scripts/jetHtSuite.py {f_ttbar} {f_nugun} {outdir}/HTcentral -j {outdir}/jecs.root -w l1pfpu_jetnoref -v ht --eta 2.4'
#cmd_puppi = f'python3 scripts/jetHtSuite.py {f_ttbar} {f_nugun} {outdir}/HT -j {outdir}/jecs.root -w l1pfpu_jetnoref -v ht --eta 5.0'

#print(cmd_jecs)
#os.system(cmd_jecs)

print(cmd_puppi)
os.system(cmd_puppi)
