import os.path
import json

file = 'config.json'
chaves = ['path', 'path_exe', 'clave3']

if not os.path.isfile(file):
    print(f'el fichero {file} no existe')
    open('config.json', 'w').close()

with open('config.json', 'r') as cfg_file:
    try:
        data = json.loads(cfg_file.read())

    except:
        data = {}

    for chave in chaves:
        if chave not in data:
            data[chave] = input(f"Necesito un valor para {chave}:")
            print(f"Gracias. Ahora {chave} vale {data[chave]}")
        else:
            respuesta = input(f"Quieres cambiar el valor de {chave}? (s/other)")
            if respuesta == "s":
                data[chave] = input(f"Necesito un valor para {chave}:")
                print(f"Gracias. Ahora {chave} vale {data[chave]}")

    json_data = json.dumps(data)
with open('config.json', 'w') as cfg_file:
    cfg_file.write(json_data)

