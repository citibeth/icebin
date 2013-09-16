import re
import os.path
import sys

license = [
"GLINT2: A Coupling Library for Ice Models and GCMs",
"",
"This program is free software: you can redistribute it and/or modify",
"it under the terms of the GNU General Public License as published by",
"the Free Software Foundation, either version 3 of the License, or",
"(at your option) any later version.",
"",
"This program is distributed in the hope that it will be useful,",
"but WITHOUT ANY WARRANTY; without even the implied warranty of",
"MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the",
"GNU General Public License for more details.",
"",
"You should have received a copy of the GNU General Public License",
"along with this program.  If not, see <http://www.gnu.org/licenses/>."
]

licenses = {}

ret = []
for line in license :
	ret.append('# ')
	ret.append(line)
	ret.append('\n')
ret.append('\n')
license_py = "".join(ret)
re_py = re.compile(r".*?\n# along with this program.  If not, see <http://www.gnu.org/licenses/>.(\n\n)", re.MULTILINE|re.DOTALL)
licenses['.py'] = (license_py, re_py)

ret = []
for line in license :
	ret.append('! ')
	ret.append(line)
	ret.append('\n')
ret.append('\n')
license_f = "".join(ret)
re_f = re.compile(r".*?\n! along with this program.  If not, see <http://www.gnu.org/licenses/>.(\n\n)", re.MULTILINE|re.DOTALL)
licenses['.f90'] = (license_f, re_f)
licenses['.f'] = (license_f, re_f)
licenses['.F90'] = (license_f, re_f)
licenses['.F'] = (license_f, re_f)

ret = ['/*\n']
for line in license :
	ret.append(' * ')
	ret.append(line)
	ret.append('\n')
ret.append('*/\n')
ret.append('\n')
re_c = re.compile(r".*?\n * along with this program.  If not, see <http://www.gnu.org/licenses/>.\n*/(\n\n)", re.MULTILINE|re.DOTALL)
license_c = "".join(ret)
licenses['.c'] = (license_c, re_c)
licenses['.h'] = (license_c, re_c)
licenses['.cpp'] = (license_c, re_c)
licenses['.hpp'] = (license_c, re_c)

# -------------------------------------
def relicense_file(fname) :
	# Determine extension
	ext = os.path.splitext(fname)[1]
	(mylicense, myre) = licenses[ext]

	# Read entire contents
	with open(fname) as file :
		txt = file.read()

	match = myre.match(txt)
	if match is None :
		body = txt
	else :
		body = txt[match.end():]
	newtxt = "".join([mylicense, body])

	if newtxt == txt :
		print '%s: Already has license' % fname
	else :
		print '%s: Adding license' % fname
		os.rename(fname, fname + '.bak')
		out = open(fname, 'w')
		out.write(newtxt)
		out.close()


for f in sys.argv[1:] :
	relicense_file(f)
