import make_PWfile_obj
from sys import argv
if len(argv) < 2:
    print('Usage: get_calculation [Path to file]')
    exit()
else:
    file = make_PWfile_obj.make_PWfile_obj(argv[1])

print(file.calculation)
