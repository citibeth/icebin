# Parses ModelE constants from Fortran source file and writes a Python
# file to generate them.  That python file is then used to create a
# netCDF file of the constants.
#
# See: write_modele_const.py
#
# by Robert Fischer        November 27, 2013

import re

# Name of the constants file to parse
constants_mod_f90 = '../model/shared/Constants_mod.F90'
# Name of output file
output_fname = 'modele_const.py'

print 'Reading ModelE constants from %s and writing to %s' % (constants_mod_f90, output_fname)

def remove_comments(fin) :
	for line in fin :
		pt = line.find('!')
		if pt >= 0 : yield line[:pt]
		else : yield line

def longlines(rawdata):
	lines = []
	for line in rawdata :
		if line is None : continue
		line = line.rstrip()
		if line.endswith("&") :
			lines.append(line[:-1].rstrip())
		else :
			lines.append(line)
			ret = " ".join(lines)
			ret = ret.replace('&', ' ')
			yield ret
			lines = []
	if len(lines) > 0:
		ret = " ".join(lines)
		ret = ret.replace('&', ' ')
		yield ret

def splitcomma(str) :
	ret = []
	i=0
	brk = 0
	parenlevel = 0
	while i < len(str) :
		if str[i] == '(' : parenlevel += 1
		if str[i] == ')' : parenlevel -= 1
		if str[i] == ',' and parenlevel == 0 :
			ret.append(str[brk:i])
			brk = i+1
		i += 1
	ret.append(str[brk:])
	return ret

# Convert Fortran to standard format floating point constants
expRE = re.compile(r'(d)-?[0-9]+')
def d2e(str) :
	lstr = list(str)
	for match in expRE.finditer(str) :
		lstr[match.start(1)] = 'e'
	ret = "".join(lstr)
	return ret

#restr = 
#
#\s*((.+?)\s*=\s*(.*))(,(/+?)\s*=\s*(.*))*'

ignore_vars = set(['nan'])

fconstRE = re.compile(r'\s+(real\*8|integer)(\s*,\s*parameter)?\s*::(.*)', re.IGNORECASE)
vardefRE = re.compile(r'\s*(.+?)\s*=\s*(.*?)\s*', re.IGNORECASE)

out = open(output_fname, 'w')
defvars = []
for line in longlines(remove_comments(open(constants_mod_f90, 'r'))) :
	if line.lstrip().startswith('contains') : break
	match = fconstRE.match(line)
	if match is None : continue
	type = match.group(1)
	for vardef in splitcomma(match.group(3)) :
		eq = vardef.split('=')
		name = eq[0].strip().lower()
		if name in ignore_vars : continue
		val = d2e(eq[1]).strip().lower()
		out.write('%s=%s\n' % (name, val))
		defvars.append(name)

# Define the set of vars at the end
out.write('defvars=%s\n' % str(defvars))

out.close()
