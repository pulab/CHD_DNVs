import os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-d","--directory", help = "input directory")
parser.add_argument("-o","--output", help = "output file")
args = parser.parse_args()

directory = args.directory
outputfilename = args.output



linevalue = lambda line: int(line[2]) - int(line[1])

filenames = os.listdir(directory)

# Prepare the data frame by finding all the annotation names
#markers, patients = list(zip(*(filename[:-4].split("__") for filename in filenames)))

data = {}

patients = set()
markers = set()

for fn in filenames:
    print(fn)
    marker, patient = fn[:-4].split("__")
    patients.add(patient)
    markers.add(marker)
    fh = open(os.path.join(directory,fn),'r')
    lines = fh.readlines()
    fh.close()
    total = 0
    for line in lines:
#        print(line)
        total = total + linevalue(line.split())
      data[(patient,marker)]=total

markers = sorted(list(markers))
patients = sorted(list(patients))


outfile = open(outputfilename,'w')

outfile.write('patient\t'+ '\t'.join(markers) + '\n')
for patient in patients:
    values = (data.get((patient,marker),"NA") for marker in markers)
    outfile.write(patient+'\t' + '\t'.join(map(str,values)) + '\n')

outfile.close()
