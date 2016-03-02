speciesFile = "bergeys_data/datasetSpeciesRefSeq.txt"
gemd = "bergeys_data/genome_metadata"
out = "bergeys_data/map.txt"

#refseqs = {}
refseqs = []
id2faa = {}
op = []

# parse input
with open(speciesFile) as file:
    for line in file:
        tmp = line.replace("\n", "").split("\t")
        refseqs.extend([el.strip() for el in tmp[1].split(",")])
       # refseqs[tmp[0]] = [el.strip() for el in tmp[1].split(",")]
        #print refseqs
refseqs = set(refseqs)

# create a dictionary
with open(gemd) as file:
    for line in file:
        tmp = line.split("\t")
        id2faa[tmp[3]] = [tmp[0], tmp[1]]

mis = 0
with open(out, "w") as file:
    file.write("sample_file_name\tsample_name\n")
    for sp in refseqs:
 #       try:
  #          sp_name = id2faa[sp][1]
   #     except KeyError:
    #        mis += 1
     #       print "Main species not found: ", sp, "Some lines skipped"
      #      continue

       # for rs in refseqs[sp]:
            #try:
             #   tmp = id2faa[rs]
              #  file.write(tmp[0] + ".RefSeq.faa\t" + sp_name + "\n")
            #except KeyError:
             #   mis += 1
        try:
            tmp = id2faa[sp]
            file.write(tmp[0] + ".RefSeq.faa\t" + tmp[1] + "\n")
        except KeyError:
            mis += 1

print "not found in file: ", mis
