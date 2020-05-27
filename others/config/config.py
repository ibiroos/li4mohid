import json

path_in = 'f'
path_exe_in = 'c:\\pedrof_exe'

with open('config.json', 'r') as cfg_file:
    try:
        data = json.loads(cfg_file.read())
    except:
        data = {}

    if 'path' not in data:
        print('1. data = ', data)
        data['path'] = path_in
        print('2. data = ', data)
    elif data['path']:
        print('4.data = ', data)
        data['path'] = path_in

    if 'path_exe' not in data:
        print('1. data = ', data)
        data['path_exe'] = path_exe_in
        print('2. data = ', data)
    else:
        data['path_exe'] = path_exe_in
        print('4.data = ', data)
    json_data = json.dumps(data)
with open('config.json', 'w') as cfg_file:
    cfg_file.write(json_data)

