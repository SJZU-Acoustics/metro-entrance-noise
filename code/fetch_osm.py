"""Fetch the major road network of the study area from the OpenStreetMap Overpass
API and cache it as one CSV of polylines (data/osm_roads_shenyang.csv) for the
point maps in Figs 1 and 5. The cache shipped with the repository was fetched on
2026-09-12; rerun only to refresh it. Data (c) OpenStreetMap contributors, ODbL. Run once;
the figure scripts read the cache and draw nothing if it is absent."""
import csv, json, pathlib, sys, time, urllib.parse, urllib.request
OUT = pathlib.Path(__file__).resolve().parents[1] / "data" / "osm_roads_shenyang.csv"  # repository data/
S, W, N, E = 41.62, 123.20, 41.99, 123.58   # study-area bbox (WGS84), all 111 stations inside
Q = f"""[out:json][timeout:180];
(
  way["highway"~"^(motorway|trunk|primary)$"]({S},{W},{N},{E});
);
out geom;"""
MIRRORS = ["https://overpass-api.de/api/interpreter",
           "https://overpass.kumi.systems/api/interpreter",
           "https://lz4.overpass-api.de/api/interpreter"]
j = None
for url in MIRRORS:
    try:
        req = urllib.request.Request(url, data=urllib.parse.urlencode({"data": Q}).encode(),
                                     headers={"User-Agent": "p11-metro-map/1.0"})
        j = json.loads(urllib.request.urlopen(req, timeout=240).read()); print("fetched from", url); break
    except Exception as e:
        print("mirror failed:", url, e); time.sleep(3)
if j is None:
    sys.exit("all Overpass mirrors failed")
rows = []
for el in j["elements"]:
    if el["type"] != "way":
        continue
    t = el.get("tags", {})
    for i, g in enumerate(el.get("geometry", [])):
        rows.append((el["id"], i, t.get("highway", ""), g["lon"], g["lat"]))
OUT.parent.mkdir(exist_ok=True)
with open(OUT, "w", newline="") as f:
    w = csv.writer(f); w.writerow(["line_id", "seq", "highway", "lng", "lat"]); w.writerows(rows)
from collections import Counter
print("ways", sum(1 for e in j["elements"] if e["type"] == "way"), "rows", len(rows), Counter(r[2] for r in rows))
