import json

# Opening JSON file
with open('thredds_config.json') as thredds_cfg:
    data = json.loads(thredds_cfg.read())

    print(data["grids"][i]["type"] for i in list(data["grids"].keys()))
    print(data["grids"]["noia"]["template"])

    hydro_models = [m for m in data["grids"].keys() if data["grids"][m]["type"] == "hydro"]
    wind_models = [m for m in data["grids"].keys() if data["grids"][m]["type"] == "wind"]

    model = 'vigo'
    url = data["grids"][model]["catalog"]
    print(url)




