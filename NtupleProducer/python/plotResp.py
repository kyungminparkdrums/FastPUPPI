import os

f_ttbar = '/eos/cms/store/group/cmst3/group/l1tr/kypark/BESTPUPPI/AR2025/TT_PU200/FP/perfTuple.root'

outdir = '/eos/user/k/kypark/www/L1BESTPUPPI/LatestGreatest/'

cmd_resp = f'python3 scripts/respPlots.py {f_ttbar} {outdir} -w l1pf -p jet'

print(cmd_resp)
os.system(cmd_resp)
