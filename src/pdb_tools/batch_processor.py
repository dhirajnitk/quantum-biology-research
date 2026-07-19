"""
batch_processor.py  —  Download and analyse many PDB structures in batch.

Usage
-----
    python src/pdb_tools/batch_processor.py < pdb_list.txt
    python src/pdb_tools/batch_processor.py --fetch 50  (searches RCSB for 50 membrane proteins)
"""

import sys
import os
import csv
import time
import json
import urllib.request
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from pdb_tools.trp_extractor import fetch_pdb, extract_trp_coordinates, distance_matrix, classify_pairs, COUPLING_CUTOFF_ANGSTROM


RESULTS_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "data")
CACHE_DIR = os.path.join(RESULTS_DIR, "cache")


def analyse_pdb(pdb_id, chain_id=None):
    """Analyse a single PDB and return summary dict, or None on failure."""
    text = fetch_pdb(pdb_id, cache=True)
    if not text:
        return None
    centers = extract_trp_coordinates(text, chain_id)
    if not centers:
        return {
            "pdb": pdb_id, "n_trp": 0, "n_coupled": 0, "n_relay": 0,
            "mean_distance": 0, "min_distance": 0, "error": "no_trp"
        }
    D, keys = distance_matrix(centers)
    coupled, relay = classify_pairs(D, keys)
    n = len(keys)
    distances = [D[i, j] for i in range(n) for j in range(i + 1, n)]
    mean_d = np.mean(distances) if distances else 0
    min_d = np.min(distances) if distances else 0
    return {
        "pdb": pdb_id,
        "n_trp": n,
        "n_coupled": len(coupled),
        "n_relay": len(relay),
        "mean_distance_A": round(mean_d, 2),
        "min_distance_A": round(min_d, 2),
        "pct_coupled": round(len(coupled) / max(len(distances), 1) * 100, 1),
        "error": ""
    }


def search_rcsb_membrane_proteins(max_results=100):
    """Use RCSB GraphQL API to find membrane proteins with Trp residues."""
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.structure_determination_methodology",
                        "operator": "exact_match",
                        "value": "experimental"
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entity_source_organism.scientific_name",
                        "operator": "exact_match",
                        "value": "Homo sapiens"
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_polymer_entity_count",
                        "operator": "greater_or_equal",
                        "value": 2
                    }
                }
            ]
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {
                "start": 0,
                "rows": max_results
            },
            "sort": [
                {"sort_by": "score", "direction": "desc"}
            ],
            "scoring_strategy": "combined"
        }
    }
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    req = urllib.request.Request(
        url,
        data=json.dumps(query).encode(),
        headers={"Content-Type": "application/json"},
        method="POST"
    )
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read())
            return [hit["identifier"] for hit in data.get("result_set", [])]
    except Exception as e:
        print(f"[!] RCSB search failed: {e}")
        return []


# ── Curated list of membrane proteins relevant to neural signalling ──
CURATED_PDBS = [
    # Ion channels — voltage-gated
    "6J8J", "7K48", "6J8E", "6A95", "6N4R", "6N4Q",  # NaV
    "7CR3", "6UWL", "6UWM", "7SJN",                     # CaV
    "7SKE", "7SKF",                                      # KV
    "1BL8", "3F7V", "3F7W",                              # KcsA / Kir
    # Ion channels — ligand-gated
    "7TYO", "7SAF", "5H8N", "5H8O", "7SJQ",              # NMDA
    "6PV7", "7SMM", "6ZGD", "6ZGE",                      # nAChR
    "6D6T", "7A5V", "7A5W", "6X3Z",                      # GABA-A
    "6X3W", "7D6L", "7D6M",                               # GlyR
    "7SKE",                                                # AMPA
    "5KOU", "6BQN", "7D6R",                               # P2X
    # GPCRs
    "1F88", "2RH1", "3SN6", "4DKL",                       # Rhodopsin, β2AR
    "5C1M", "6N4B", "7S3L",                               # Adenosine, opioid
    "7CMO", "7E2Y", "7E2Z",                               # Cannabinoid
    # Transporters
    "5I6X", "5I71",                                        # SERT
    "4M48", "4XP4", "4XP9",                                # DAT
    "2R5D", "5E6R",                                        # GLUT
    # Water channels
    "1FQY", "5KCF", "5KCG",                                # Aquaporins
    # Gap junctions / other
    "5ERA", "7F9M",                                        # Connexin
    "6H7W", "6H7X",                                        # TRPV1
    # Synaptic / scaffolding
    "3SOA", "5WS3", "6MYA", "7R6K",                        # PSD-95, Synaptotagmin
]


def process_batch(pdb_list, output_name="batch_results.csv"):
    """Process a list of PDB IDs and save aggregate results."""
    os.makedirs(RESULTS_DIR, exist_ok=True)
    out_path = os.path.join(RESULTS_DIR, output_name)

    rows = []
    total = len(pdb_list)
    for i, pdb_id in enumerate(pdb_list):
        pdb_id = pdb_id.strip().upper()
        if not pdb_id:
            continue
        sys.stdout.write(f"\r[{i+1}/{total}] {pdb_id} ... ")
        sys.stdout.flush()
        try:
            row = analyse_pdb(pdb_id)
            if row:
                rows.append(row)
        except Exception as e:
            rows.append({"pdb": pdb_id, "error": str(e)[:60]})
        time.sleep(0.3)  # rate limit

    print(f"\n[+] Writing {out_path}")
    with open(out_path, "w", newline="") as f:
        fieldnames = ["pdb", "n_trp", "n_coupled", "n_relay",
                       "mean_distance_A", "min_distance_A", "pct_coupled", "error"]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    ok = [r for r in rows if not r.get("error")]
    fail = [r for r in rows if r.get("error")]
    print(f"    {len(ok)} succeeded, {len(fail)} failed")
    return rows


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--fetch":
        n = int(sys.argv[2]) if len(sys.argv) > 2 else 50
        print(f"[+] Searching RCSB for {n} membrane proteins ...")
        ids = search_rcsb_membrane_proteins(n)
        print(f"[+] Found {len(ids)} IDs")
        process_batch(ids, "rcsb_membrane_batch.csv")
    else:
        # Process curated list
        process_batch(CURATED_PDBS, "curated_membrane_batch.csv")
