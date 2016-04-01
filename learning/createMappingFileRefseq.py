speciesFile = "tmp.txt"
#speciesFile = "bergeys_data/datasetSpeciesRefSeq.txt"
#gemd = "bergeys_data/genome_metadata"
bioprojects = "/home/kyra/Arbeit/repositories/traitar-model/learning/bergeys_data/bioprojects_20140115.txt"
gemd = "/home/kyra/Arbeit/repositories/traitar-model/learning/bergeys_data/genome_metadata"
#out = "bergeys_data/map_v2.txt"
out = "/home/kyra/Dropbox/Arbeit/traitar/paperMaptmp.txt"
out2 = "/home/kyra/Dropbox/Arbeit/traitar/refseqpr2ncbi.txt"

refseqs = {}
#refseqs = []
id2faa = {}
op = []

# parse input
with open(speciesFile) as file:
    for line in file:
        tmp = line.replace("\n", "").split("\t")
        #refseqs.extend([el.strip() for el in tmp[1].split(",")])
        refseqs[tmp[0]] = [el.strip() for el in tmp[1].split(",")]
        #print refseqs
#refseqs = set(refseqs)

# create a dictionary mapping refseq project ids to ncbi ids and strain names
with open(bioprojects) as file:
    for line in file:
        tmp = line.split("\t")
        id2faa[tmp[3]] = [tmp[0], tmp[1]]

# create a dictionary to map ncbi_id to filename
gemddic = {}
with open(gemd) as file:
    for line in file:
        tmp = line.split("\t")
        gemddic[tmp[3]] = tmp[0]

mis = 0
of = open(out2,"w")
with open(out, "w") as file:
    file.write("sample_file_name\tsample_name\n")
    for sp in refseqs:
        for rs in refseqs[sp]:
            t = "?"
            try:
                tmp = id2faa[rs]
                try:
                    t = gemddic[tmp[1]]
                    file.write(t + ".RefSeq.faa\t" + tmp[1] + "\n")
                    of.write(rs + "\t" + tmp[1] + "\n")
                except KeyError:
                    print "no filename found in genome_metadata. NCBI: ", tmp[1]

            except KeyError:
                mis += 1
                #print sp

print "not found in file: ", mis
print len(refseqs) - mis