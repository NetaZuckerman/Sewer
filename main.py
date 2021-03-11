import pysam

file = pysam.AlignmentFile('s7510_S483_L002_001.mapped.sorted.bam', "rb")
iter = file.pileup()
mystring = ''
mylist = []
for pileupcol in file.pileup():
    mystring = ''
    for pileupread in pileupcol.pileups:
        if not pileupread.is_del and not pileupread.is_refskip:
            mystring = mystring + pileupread.alignment.query_sequence[pileupread.query_position]

    mystring = "".join(sorted(list(mystring)))
    mylist.append((pileupcol.pos+1, mystring))
    for (k, v) in mylist:
        print("{}: {}".format(k, len(v)))


for (k, v) in mylist:
    print("{}: {}".format(k, len(v)))
