import os
import sys
import subprocess

patientdir = sys.argv[1]
markerdir = sys.argv[2]
bednames = sys.argv[3]
outputdir = sys.argv[4]

patient_fns = os.listdir(patientdir)
marker_fns = file(os.path.join(markerdir,bednames),'r').readlines()

stem = lambda fn:  "_".join(fn.split(".")[:-1])

for marker_fn in marker_fns:
    for patient_fn in patient_fns:
        pat_fq = os.path.join(patientdir,patient_fn)
        mar_fq = os.path.join(markerdir,marker_fn).strip()
        out_fn = stem(marker_fn) +"__" + stem(patient_fn) + ".bed"
        out_fq = os.path.join(outputdir,out_fn)
        cmd = " ".join(['bedtools','intersect','-a',mar_fq,'-b',pat_fq,'>',out_fq])
        print(cmd)
        if not os.path.exists(out_fq):
            subprocess.call(cmd, shell=True)
