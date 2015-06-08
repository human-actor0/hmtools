#!/usr/bin/env python
import sys,getopt

usage = """
Tool: demultiplex orignal fastq 
Version: v0.1
Author: Hyunmin Kim (Hyun.Kim@ucdenver.edu)
Output: demultiplexed fastq matching barcodes are chopped
Usage:  %s [options] <fastq> 
	<fastq>: stdin is reserved for the standard input

	[options]:
	 -3 <barcode>: 3'end barcode 
	 -5 <barcode>: 5'end barcode 
""" % (sys.argv[0])


def main(argv):
	threebar=None 
	fivebar=None
	fastq = 'stdin'

	try:
		opts,args = getopt.getopt(argv,'h5:3:')
		if len(args) > 1:
			fastq = args[0]
	except getopt.GetoptError, err:
		print usageo; sys.exit(-1)

	for o,a in opts:	
		if o=='-h':
			print usage; sys.exit(-1)
		elif o=='-3': threebar = a		
		elif o=='-5': fivebar = a		
		else:
			raise getopt.GetoptError, "unknown option"

	if fastq == 'stdin': fh = sys.stdin
	else:fh = open(fastq,'r')
	buf = []
	while 1:
		line = fh.readline().strip()
		if line =='': break
		buf.append(line)
		if len(buf) == 4:
			# handle buf
			#print '\n'.join(buf)
			s = buf[1]
			q = buf[3]
			if fivebar != None:
				le = len(fivebar)
				if fivebar != s[:le]:
					buf = []
					continue
				buf[0] += '__F_'+s[:le]+'_'+q[:le]
				buf[1] =buf[1][le:]
				buf[3] =buf[3][le:]
			if threebar != None:
				le = len(threebar)
				if threebar != s[-le:]:
					buf = []
					continue
				buf[0] += '__T_'+s[-le:]+'_'+q[-le:]
				buf[1] = buf[1][:-le]
				buf[3] = buf[3][:-le]
			print '\n'.join(buf)
			buf = []

	if fastq != 'stdin': fh.close()
		


if __name__=='__main__':
	main(sys.argv[1:])	

