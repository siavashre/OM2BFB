import os
a = ['BOKU', 'CaSki', 'Detroit562', 'HCS2', 'HN120Met','HN120Pri','HN148_Met', 'HN148_Pri', 'HN159_Met', 'HN159_Pri', 'HN160Met', 'HN160Pri', 'OSC20']
for i in a:
    cmd = 'mkdir /home/sraeisid/test/{name}'.format(name = i)
    print(cmd)
    os.system(cmd)
    cmd = "python3 /home/sraeisid/OM2BFB/runOM2BFB.py  -r /nucleus/projects/sraeisid/BFB/Daniel/{name}/copynumber/cnv_rcmap_exp.txt -c /home/sraeisid/OM2BFB/hg38_centro.txt -n {name} -o /home/sraeisid/test/{name}/ -s /nucleus/projects/sraeisid/BFB/Daniel/{name}/merged_smaps/exp_refineFinal1_merged_filter_inversions.smap -x /nucleus/projects/sraeisid/BFB/Daniel/{name}/merged_smaps/exp_refineFinal1_merged.xmap -cmap /nucleus/projects/sraeisid/BFB/Daniel/{name}/merged_smaps/exp_refineFinal1_merged_q.cmap -fol /nucleus/projects/sraeisid/BFB/Daniel/{name}/refine1_ExperimentLabel/  -bfbfinder /home/sraeisid/libraries/BFBFinder/BFBFinder/build/libs/BFBFinder.jar".format(name = i)
    print(cmd)
    os.system(cmd)
