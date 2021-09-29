
#read in each line of ASO position file
file = open("RemovedMissing")
positiondata = file.readlines()
file.close()

#loop through position data, open appropriate ct file, define the ASO position, count up basepairs
for line in positiondata:
	transcript,start,end,asoid = line.split()	
	start = int(start)
	end = int(end)
	bigstart = start
	bigend = end
	print(transcript,asoid,start,end)
	ctfile1name = transcript + ".ct"
	#print(ctfile1name)
	ctfile1 = open(ctfile1name)
	header = ctfile1.readline()
	structuredata = ctfile1.readlines()
	ctfile1.close()
	print(header)
	length = 0
	bp = 0
	un = 0
	resultfile = "mod_" + transcript + "_" + asoid + ".ct"
	f = open(resultfile, "a")
	print(header, end ='', file=f)
	for basepair in structuredata:
		num,base,val1,val2,pairedness,val4 = basepair.split()
		num = int(num)
		newpairedness = 0
		pairedness = int(pairedness)
		if bigend >= num >= bigstart:
			print(num, base, val1, val2, newpairedness, val4, file=f)
		elif bigend >= pairedness >= bigstart: 
			print(num, base, val1, val2, newpairedness, val4, file=f)
		else:	
			print(num, base, val1, val2, pairedness, val4, file=f)
	f.close()
