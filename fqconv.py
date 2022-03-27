#!/usr/bin/env python

#######################################
# methylCtools fqconv
# v1.1.0
# 26 march 2022
#
# volker hovestadt
# developed at the german cancer research center, 2011-2015
# methylctools@hovestadt.bio
#
#
# converts illumina BS sequencing reads to 3-letter alphabet.
# reads must be in fastq format. supports gzip-compressed input files.
# to read an uncompressed file from stdin, use "-".
# to write to stdout, use "-".
#
# for SE: reads are C to T converted (define input file using -1 argument).
# for PE: first end reads are C to T converted, second-end reads are G to A converted (-2 argument).
# define both ends to create an interleaved output file (bwa-mem compatible).
#
# read ids are trimmed to the first occurrence of # or space character.
# conversion positions are appended to read id (in hexadecimal).
# if the resulting read id becomes longer than 250 chars (violates SAM
# format specifications, could happen with very long reads, the read is not converted.


def mod_fqconv(sysargv):
	import sys
	import argparse
	import re
	import gzip
	import datetime
	def nicetime(): return datetime.datetime.now().strftime("[fqconv %Y-%m-%d %H:%M:%S]")
	
	
	#######################################
	# arguments, filehandles

	parser = argparse.ArgumentParser(prog="methylCtools fqconv", version="0.9.4", description="converts bisulfite sequencing reads to fully converted state")
	parser.add_argument("-s", "--silent", dest="qf", action="store_false", help="do not show status messages")

	groupinput = parser.add_argument_group("input files, required (fastq format, \"-\" for stdin, supports gzip compressed files)")
	groupinput.add_argument("-1", dest="inFQ1", metavar="reads1.fq", action="store", help="first read in pair or single end")
	groupinput.add_argument("-2", dest="inFQ2", metavar="reads2.fq", action="store", help="second read in pair")

	groupoutput = parser.add_argument_group("output file, will be created (fastq format)")
	groupoutput.add_argument("outFQ", metavar="reads.conv.fq", action="store", default=False, help="converted sequencing reads, \"-\" for stdout")

	args = parser.parse_args(sysargv)
		
	try:
		if args.inFQ1:
			f1 = True
			if args.inFQ1.split(".")[-1].lower() in ["gz", "gzip"]: in1 = gzip.open(args.inFQ1, "rb")
			else:
				if args.inFQ1 == "-": in1 = sys.stdin
				else: in1 = open(args.inFQ1, "r")
		else: f1 = False
		
		if args.inFQ2:
			f2 = True
			if args.inFQ2.split(".")[-1].lower() in ["gz", "gzip"]: in2 = gzip.open(args.inFQ2, "rb")
			else:
				if args.inFQ2 == "-": in2 = sys.stdin
				else: in2 = open(args.inFQ2, "r")
		else: f2 = False
				
		if args.outFQ == "-": out1 = sys.stdout
		else: out1 = open(args.outFQ, "w")
		
	except IOError as strerror:
		sys.exit("%s error: %s" % (nicetime(), strerror))
	
	if not f1 and not f2: sys.exit("%s error: -1 or/and -2 must be defined.", nicetime())

	if args.qf: sys.stderr.write("%s command: %s\n" % (nicetime(), " ".join(sys.argv)))
	

	#######################################
	# main

	cfrom1, cto1 = "C", "T"
	cfrom2, cto2 = "G", "A"

	n = 0
	l1, l2 = "", ""
	id1, id2 = "", ""
	seq1, seq2 = "", ""
	qual1, qual2 = "", ""
        
	if f1 and args.qf: sys.stderr.write("%s start: converting reads %s to %s\n" % (nicetime(), cfrom1, cto1))
	if f2 and args.qf: sys.stderr.write("%s start: converting reads %s to %s\n" % (nicetime(), cfrom2, cto2))
	c = [0, 0]																# [reads, bases converted]
	wc = 0																	# warning count
	
	o12 = ""
	while True:
		if f1: l1 = in1.readline().rstrip()									# read line
		if f2: l2 = in2.readline().rstrip()
		
		if f1:																# check if file ends
			if not l1:
				if n: sys.exit("%s error: %s is truncated" % (nicetime(), args.inFQ1))
				elif l2: sys.exit("%s error: %s is truncated" % (nicetime(), args.inFQ1))
				else: break
		if f2:
			if not l2:
				if n: sys.exit("%s error: %s is truncated" % (nicetime(), args.inFQ2))
				elif l1: sys.exit("%s error: %s is truncated" % (nicetime(), args.inFQ2))
				else: break
				
		n += 1
		if n == 1:															# line 1: name
			if f1 and f2: c[0] += 2
			else: c[0] += 1
			
			if f1: id = l1.split()[0].split("#")[0]							# BWA only uses id up to first space, clip # (not picard-compatible?)
			if not f1: id = l2.split()[0].split("#")[0]						# assume same name for paired-end reads
			
		elif n == 2:														# line 2: sequence
			if f1: seq1 = l1
			if f2: seq2 = l2
																			# (line 3: quality name, replace to "+")		
		elif n == 4:														# line 4: quality
			if f1: qual1 = l1
			if f2: qual2 = l2

			if f1:															# convert sequence and append positions to id (in hexadecimal)
				id += "."
				c[1] += seq1.count(cfrom1)
				id += hex(int(re.sub("\D", "0", re.sub(cfrom1, "1", seq1)), 2))[2:-1]
				if id[-1] == ".": id += "0"
				seq1 = re.sub(cfrom1, cto1, seq1)
			if f2:
				id += "."
				c[1] += seq2.count(cfrom2)
				id += hex(int(re.sub("\D", "0", re.sub(cfrom2, "1", seq2)), 2))[2:-1]
				if id[-1] == ".": id += "0"
				seq2 = re.sub(cfrom2, cto2, seq2)	
				
			if len(id) >= 250:												# if id too long (SAM allows 255 characters)
				id = id.split(".")[0] 
				if f1: id += ".0"
				if f2: id += ".0"
				if wc < 100: sys.stderr.write("%s warning: %s is not converted (%s)\n" % (nicetime(), id, seq1))
				if wc == 99: sys.stderr.write("%s warning: only showing 100 warnings\n" % nicetime())
				wc += 1
			else:
				if f1: seq1 = seq1.replace(cfrom1, cto1)					# convert
				if f2: seq2 = seq2.replace(cfrom2, cto2)
		
			if f1: o12 += ("%s\n%s\n+\n%s\n" % (id, seq1, qual1))			# append
			if f2: o12 += ("%s\n%s\n+\n%s\n" % (id, seq2, qual2))
			n = 0															# reset counter

			if c[0]%10**6 == 0:
				out1.write(o12)												# write 100k reads at once
				o12 = ""
				if args.qf: sys.stderr.write("%s status: %i reads processed\n" % (nicetime(), c[0]))

	out1.write(o12)
		

	#######################################
	# end

	if f1: in1.close()
	if f2: in2.close()
	out1.close()

	if args.qf: sys.stderr.write("%s end: %i reads processed, %i bases converted\n" % (nicetime(), c[0], c[1]))


if __name__ == "__main__":
    import sys
    mod_fqconv(sys.argv[1:])


