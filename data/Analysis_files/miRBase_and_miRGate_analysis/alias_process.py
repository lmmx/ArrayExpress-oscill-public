#!/usr/bin/env python2.7

import sys
aliases=sys.stdin.read().split("\n")
for linenum in xrange(len(aliases)):
	line = aliases[linenum]
	fields = line.split("\t")
	mir_accession = fields[0]
	if len(fields) > 1:
		mir_names = filter(None, fields[1].split(";"))
		num_names = len(mir_names)
		for name_num in xrange(len(mir_names)):
			print mir_accession+"\t"+mir_names[name_num]
