#!/usr/bin/env python3

import sys
from json import load
from typing import Annotated
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import urlopen

from mcp.server import MCPServer
from pydantic import Field


# Create the MCP server.
mcp = MCPServer("apoholo-workshop")


@mcp.tool(title="AhojDB: Search apo/holo protein structures")
def search_apoholo(
    pdb_ids: Annotated[
        str,
        Field(
            description=(
                "Optional. Comma-separated PDB IDs to look up in AhojDB. "
                "Example: 1a73,2hhb"
            ),
        ),
    ] = "",
    uniprot_ids: Annotated[
        str,
        Field(
            description="Optional. Comma-separated UniProt accessions. Example: Q94702",
        ),
    ] = "",
    ligands: Annotated[
        str,
        Field(
            description=(
                "Optional. Comma-separated ligand codes (PDB 3-letter chemical "
                "component IDs). Example: ZN,HEM"
            ),
        ),
    ] = "",
    pdb_limit: Annotated[
        int,
        Field(
            description=(
                "Optional. Maximum number of apo/holo PDB IDs to list per entry. "
                "The exact total counts are always returned regardless of this "
                "limit. Only raise it if the user explicitly wants the full list. "
                "Default: 10."
            ),
            ge=1,
        ),
    ] = 10,
) -> dict:
    """
    Searches AhojDB (apoholo.cz) for precomputed apo (ligand-free) and holo
    (ligand-bound) forms of protein structures.

    Provide at least one of pdb_ids, uniprot_ids, or ligands. Each matching
    entry describes one binding pocket, with counts of the apo and holo
    structures found for it and a sample of the PDB IDs of those structures
    (up to pdb_limit each; the exact totals are always reported).
    """
    print(
        "Tool search_apoholo called with parameters: "
        f"pdb_ids={pdb_ids}, uniprot_ids={uniprot_ids}, ligands={ligands}, "
        f"pdb_limit={pdb_limit}",
        file=sys.stderr,
    )

    if not (pdb_ids.strip() or uniprot_ids.strip() or ligands.strip()):
        return {"error": "Provide at least one of pdb_ids, uniprot_ids, or ligands."}

    parameters = urlencode(
        {
            "pdb_ids": pdb_ids.strip(),
            "uniprot_ids": uniprot_ids.strip(),
            "ligands": ligands.strip(),
        }
    )
    url = f"https://apoholo.cz/api/db/search?{parameters}"

    try:
        with urlopen(url, timeout=30) as response:
            result = load(response)
    except HTTPError as error:
        return {"error": f"Could not search AhojDB (HTTP {error.code})."}
    except URLError:
        return {"error": "Could not reach AhojDB."}

    entries = []
    for entry in result.get("entries", []):
        apo = sorted(set(entry.get("found_apo_pdbids") or []))
        holo = sorted(set(entry.get("found_holo_pdbids") or []))
        entries.append(
            {
                "entry_key": entry.get("entry_key"),
                "query": entry.get("query"),
                "pdb_id": entry.get("target_pdb_id"),
                "ligand": entry.get("target_ligand"),
                "uniprot_ids": entry.get("target_uniprot_ids"),
                "assignment": entry.get("target_apoholo_assignment"),  # A=apo, H=holo
                "resolution": entry.get("target_resolution"),
                "num_apo_pdbs": entry.get("num_apo_pdbids"),
                "num_holo_pdbs": entry.get("num_holo_pdbids"),
                "apo_pdbs_sample": apo[:pdb_limit],
                "holo_pdbs_sample": holo[:pdb_limit],
                "pdbs_truncated": len(apo) > pdb_limit or len(holo) > pdb_limit,
            }
        )

    return {
        "query": {"pdb_ids": pdb_ids, "uniprot_ids": uniprot_ids, "ligands": ligands},
        "num_entries": len(entries),
        "entries": entries,
    }


if __name__ == "__main__":
    try:
        mcp.run(
            transport="streamable-http",
            host="0.0.0.0",
            port=8000,
            stateless_http=True,
        )
    except KeyboardInterrupt:
        pass
