
#read in each line of ASO position file
file = open("RemovedMissing")
positiondata = file.readlines()
file.close()

#loop through position data, open appropriate ct file, define the ASO position, count up basepairs
for line in positiondata:
	#print(line)
	transcript,start,end,asoid = line.split()	
	start = int(start)
	end = int(end)
	bigstart = start - 10
	bigend = end + 10
	#print(start,end)
	ctfile1name = "mod_" + transcript + "_" + asoid + ".ct"
	#print(ctfile1name)
	ctfile1 = open(ctfile1name)
	structuredata = ctfile1.readlines()[1:]
	ctfile1.close()
	length = 0
	bp = 0
	un = 0
	for basepair in structuredata:
		num,base,val1,val2,pairedness,val4 = basepair.split()
		num = int(num)
		pairedness = int(pairedness)
		if bigend >= num >= bigstart:
			#print(start, end, num, pairedness)
			length += 1
			if pairedness == 0:
				un += 1
			else:
				bp += 1
	print(asoid, transcript, start, end, length, un, bp)
