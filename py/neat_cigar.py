import re

class CigarString:
	def __init__(self, stringIn=None, listIn=None):

		if stringIn == None and listIn == None:
			print '\nError: CigarString object not initialized.\n'
			exit(1)

		self.cigarData = []

		if stringIn != None:
			self.joinCigar(j_stringIn=stringIn)

		if listIn != None:
			self.joinCigar(j_listIn=listIn)


	def stringToList(self, s):

		cigarDat = []
		letters = re.split(r"\d+",s)[1:]
		numbers = [int(n) for n in re.findall(r"\d+",s)]
		dReserve = 0
		for i in xrange(len(letters)):
			if letters[i] == 'D':
				dReserve = numbers[i]
			if letters[i] == 'M' or letters[i] == 'I':
				if dReserve:
					cigarDat += ['D'*dReserve+letters[i]] + [letters[i]]*(int(numbers[i])-1)
				else:
					cigarDat += [letters[i]]*int(numbers[i])
				dReserve = 0
		return cigarDat


	def listToString(self, l):

		symbols      = ''
		currentSym   = l[0]
		currentCount = 1
		if 'D' in currentSym:
			currentSym   = currentSym[-1]
		for k in xrange(1,len(l)):
			nextSym = l[k]
			if len(nextSym) == 1 and nextSym == currentSym:
				currentCount += 1
			else:
				symbols += str(currentCount) + currentSym
				if 'D' in nextSym:
					symbols += str(nextSym.count('D')) + 'D'
					currentSym   = nextSym[-1]
				else:
					currentSym   = nextSym
				currentCount = 1
		symbols += str(currentCount) + currentSym
		return symbols

	def getList(self):

		return self.cigarData


	def getString(self):

		return self.listToString(self.cigarData)


	def joinCigar(self, j_stringIn=None, j_listIn=None):

		if j_stringIn == None and j_listIn == None:
			print '\nError: Invalid join operation in CigarString\n'
			exit(1)

		if j_stringIn != None:
			self.cigarData += self.stringToList(j_stringIn)

		if j_listIn != None:
			self.cigarData += j_listIn


	def insertCigarElement(self, pos, i_stringIn=None, i_listIn=None):

		if i_stringIn == None and i_listIn == None:
			print '\nError: Invalid insertion operation in CigarString\n'
			exit(1)

		if pos < 0 or pos >= len(self.cigarData):
			print '\nError: Invalid insertion position in CigarString\n'
			exit(1)

		if i_stringIn != None:
			self.cigarData = self.cigarData[:pos] + self.stringToList(i_stringIn) + self.cigarData[pos:]

		if i_listIn != None:
			self.cigarData = self.cigarData[:pos] + i_listIn + self.cigarData[pos:]


if __name__ == '__main__':
	print 'testing CigarString class...'

	str1 = '50M10D7I23M'
	str2 = '10I25M'
	iPos = 20
	myCigar  = CigarString(stringIn=str1)
	myCigar.insertCigarElement(iPos,i_stringIn=str2)
	print str1,'+',str2,'[inserted at position',str(iPos)+']','=',myCigar.getString()

