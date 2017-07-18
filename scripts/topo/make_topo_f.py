from giss import ioutil
import subprocess
import os

# Find one input file, and like that as a Z directory
def link_zdir():
    for path in os.environ['MODELE_FILE_PATH'].split(os.pathsep):
        if os.path.exists(os.path.join(path, 'Z2MX2M.NGDC')):
            os.symlink(path, 'Z')
            return 
try:
    link_zdir()
except FileExistsError:
    pass


os.environ['IC'] = os.getcwd()
cmd = [ioutil.search_file('make_topo_f', os.environ['PATH'])]

print('cmd', cmd)

with open('make_topo_f.out', 'w') as fout:
    with open('make_topo_f.err', 'w') as ferr:
        subprocess.call(cmd, stdout=fout, stderr=ferr)

cmd = ['giss2nc', 'Z1QX1F.BS1', 'Z1QX1F.BS1.nc']
subprocess.check_call(cmd)

cmd = ['giss2nc', 'Z1QX1N.BS1', 'Z1QX1N.BS1.nc']
subprocess.check_call(cmd)

cmd = ['giss2nc', 'int_f', 'int_f.nc']
subprocess.check_call(cmd)


