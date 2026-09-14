#!/usr/bin/env python3

import os
import re
import sys
from json import load
from typing import Annotated
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlencode
from urllib.request import urlopen

from mcp.server import MCPServer
from pydantic import Field


# Create the MCP server.
mcp = MCPServer("apoholo-workshop")


def log_call(tool: str, **params: object) -> None:
    """One line per tool invocation on stderr (captured by journald as a service)."""
    args = ", ".join(f"{key}={value}" for key, value in params.items())
    print(f"tool {tool} called: {args}", file=sys.stderr, flush=True)


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

    Use this to START an investigation from a bare PDB ID, UniProt accession,
    or ligand code that has not yet been looked up in this conversation. It
    returns entry_key handles (form: PDB-CHAIN-LIGAND-POS, e.g. 1a01-A-HEM-142)
    that the drill-down tools (get_pocket_residues, find_conformers) then use.
    Do NOT use this to get more detail on a structure already listed in a
    previous result — use a drill-down tool with its entry_key instead.

    AhojDB has no free-text name search. To search by a protein NAME (e.g.
    'hemoglobin'), first call resolve_protein to turn it into UniProt
    accessions, then pass those as uniprot_ids here.

    Provide at least one of pdb_ids, uniprot_ids, or ligands. Two output modes:
    - By uniprot_ids only (a whole protein): a single accession can map to
      hundreds of structures, so results are summarized by bound ligand
      (mode='uniprot_ligand_summary') rather than listed row by row.
    - By pdb_ids and/or ligands: one entry per matching binding pocket, with
      apo/holo counts and a sample of PDB IDs (up to pdb_limit each; exact
      totals always reported).
    """
    log_call(
        "search_apoholo",
        pdb_ids=pdb_ids,
        uniprot_ids=uniprot_ids,
        ligands=ligands,
        pdb_limit=pdb_limit,
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

    entries_raw = result.get("entries", [])

    # Protein-level query (by UniProt accession only): one accession can map to
    # hundreds of entries -- the same pocket repeated across many PDBs -- so
    # summarize by bound ligand instead of returning every row.
    if uniprot_ids.strip() and not pdb_ids.strip():
        groups: dict = {}
        for entry in entries_raw:
            ligand = entry.get("target_ligand")
            group = groups.setdefault(
                ligand, {"pdb_ids": set(), "positions": set(), "entry_keys": []}
            )
            group["pdb_ids"].add(entry.get("target_pdb_id"))
            parts = (entry.get("entry_key") or "").split("-")
            if len(parts) >= 4:
                group["positions"].add(parts[-1])  # residue position of the ligand
            group["entry_keys"].append(entry.get("entry_key"))

        ligand_groups = sorted(
            (
                {
                    "ligand": ligand,
                    "num_structures": len(group["entry_keys"]),
                    "num_pdbs": len(group["pdb_ids"]),
                    "num_positions": len(group["positions"]),
                    "sample_entry_keys": group["entry_keys"][:pdb_limit],
                }
                for ligand, group in groups.items()
            ),
            key=lambda g: g["num_structures"],
            reverse=True,
        )

        return {
            "query": {"uniprot_ids": uniprot_ids},
            "mode": "uniprot_ligand_summary",
            "total_entries": len(entries_raw),
            "total_pdbs": len({e.get("target_pdb_id") for e in entries_raw}),
            "num_ligands": len(ligand_groups),
            "note": (
                "Results for a whole UniProt accession, grouped by bound ligand. "
                "A ligand group can span several binding pockets and residue "
                "positions (num_positions) -- it is NOT necessarily one pocket, "
                "and residue numbering also varies across structures. Use "
                "sample_entry_keys with the drill-down tools, or search_apoholo "
                "by pdb_ids, to inspect specific pockets."
            ),
            "ligand_groups": ligand_groups[:pdb_limit],
            "ligands_truncated": len(ligand_groups) > pdb_limit,
        }

    entries = []
    for entry in entries_raw:
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


@mcp.tool(title="UniProt: Resolve a protein name to accessions")
def resolve_protein(
    name: Annotated[
        str,
        Field(
            description=(
                "Required. A protein name, gene name, or synonym to look up, "
                "e.g. 'hemoglobin', 'cytochrome P450', 'HBB'."
            ),
        ),
    ],
    organism_id: Annotated[
        int,
        Field(
            description=(
                "Optional. NCBI taxonomy id to restrict the search. "
                "Default: 9606 (human)."
            ),
            ge=1,
        ),
    ] = 9606,
    limit: Annotated[
        int,
        Field(
            description="Optional. Maximum number of candidate accessions. Default: 10.",
            ge=1,
        ),
    ] = 10,
) -> dict:
    """
    Resolves a protein name (or gene/synonym) to candidate UniProt accessions
    for use as the uniprot_ids input of search_apoholo.

    AhojDB is keyed by PDB structures and UniProt accessions, not free-text
    names, so use this FIRST to turn a name like 'hemoglobin' into accessions.
    A name is often ambiguous: it can map to several accessions -- protein
    subunits (hemoglobin alpha/beta/...) or isoforms (many cytochrome P450s).
    Choose the relevant one(s), ideally confirming with the user, then search.

    Results are restricted to reviewed (Swiss-Prot) entries of the given
    organism (human by default). Tell the user that the human default was
    applied and can be changed via organism_id.
    """
    log_call("resolve_protein", name=name, organism_id=organism_id, limit=limit)

    name = name.strip()
    if not name:
        return {"error": "name must not be empty."}

    query = f"{name} AND reviewed:true AND organism_id:{organism_id}"
    parameters = urlencode(
        {
            "query": query,
            "fields": "accession,id,protein_name,gene_names,organism_name",
            "format": "json",
            "size": limit,
        }
    )
    url = f"https://rest.uniprot.org/uniprotkb/search?{parameters}"

    try:
        with urlopen(url, timeout=30) as response:
            data = load(response)
    except HTTPError as error:
        return {"error": f"Could not query UniProt (HTTP {error.code})."}
    except URLError:
        return {"error": "Could not reach UniProt."}

    candidates = []
    for record in data.get("results", []):
        description = record.get("proteinDescription", {})
        protein_name = (
            description.get("recommendedName", {})
            .get("fullName", {})
            .get("value")
        )
        if not protein_name:
            submitted = description.get("submissionNames") or []
            if submitted:
                protein_name = submitted[0].get("fullName", {}).get("value")
        genes = [
            gene.get("geneName", {}).get("value")
            for gene in record.get("genes", [])
            if gene.get("geneName")
        ]
        candidates.append(
            {
                "accession": record.get("primaryAccession"),
                "protein_name": protein_name,
                "genes": genes,
                "organism": record.get("organism", {}).get("scientificName"),
            }
        )

    return {
        "name": name,
        "organism_id": organism_id,
        "note": (
            f"Searched reviewed UniProt entries for organism {organism_id} "
            "(9606 = human, the default). Tell the user the organism filter was "
            "applied and can be changed. A name can map to several accessions "
            "(subunits or isoforms); pick the relevant one(s) before searching."
        ),
        "num_candidates": len(candidates),
        "candidates": candidates,
    }


@mcp.tool(title="AhojDB: Get binding-pocket residues")
def get_pocket_residues(
    entry_key: Annotated[
        str,
        Field(
            description=(
                "Required. An AhojDB entry key from a search_apoholo result, "
                "e.g. 1a01-A-HEM-142. Identifies one binding pocket."
            ),
        ),
    ],
    pdb_id: Annotated[
        str,
        Field(
            description=(
                "Optional. A specific apo or holo PDB ID from that entry to get "
                "the pocket residues for. Leave empty to get the residues of the "
                "query pocket itself."
            ),
        ),
    ] = "",
) -> dict:
    """
    Returns the binding-site (pocket) residues for one AhojDB entry.

    Use this ONLY to drill into an entry you already found with search_apoholo.
    It requires an entry_key (e.g. 1a01-A-HEM-142), never a bare 4-character PDB
    ID; to inspect one specific apo/holo structure from the entry, pass it as
    pdb_id.

    The residues are parsed from the entry's processing log. With no pdb_id the
    residues of the query pocket are returned (residue name + number). With a
    pdb_id, a query->conformer mapping is returned: each query-pocket residue
    paired with its counterpart in that structure (matched by UniProt residue
    number), flagging residues that are unobserved in the conformer.
    """
    log_call("get_pocket_residues", entry_key=entry_key, pdb_id=pdb_id)

    entry_key = entry_key.strip()
    if not entry_key:
        return {"error": "entry_key must not be empty."}

    url = f"https://apoholo.cz/api/db/entry/{quote(entry_key, safe='')}/log"
    try:
        with urlopen(url, timeout=30) as response:
            log = response.read().decode("utf-8", "replace")
    except HTTPError as error:
        return {"error": f"Could not retrieve entry log (HTTP {error.code})."}
    except URLError:
        return {"error": "Could not reach AhojDB."}

    # Query-pocket residues: the first residue list in the protein
    # binding-site section, e.g. "A_HEM_142: [21] ['A_A_MET_32', ...]".
    match = re.search(
        r"binding_site_residues_prot\)\s*\n[^\n]*?:\s*\[\d+\]\s*(\[[^\]]*\])", log
    )
    if not match:
        return {"error": f"No pocket residues found for entry {entry_key}."}

    query_residues = []
    for token in re.findall(r"'([^']+)'", match.group(1)):
        parts = token.rsplit("_", 2)  # e.g. A_A_MET_32 -> ['A_A', 'MET', '32']
        if len(parts) < 3 or not parts[-1].isdigit():
            continue
        name, number = parts[-2], int(parts[-1])
        query_residues.append({"name": name, "number": number})

    pdb_id = pdb_id.strip()
    if not pdb_id:
        return {
            "entry_key": entry_key,
            "pdb_id": entry_key.split("-", 1)[0],
            "is_query_pocket": True,
            "count": len(query_residues),
            "residues": query_residues,
        }

    # Specific structure: residue numbers from its PyMOL selection line,
    # e.g. "2w72 and ( (chain A and resi 32+39+42+...".
    target = re.search(
        rf"{re.escape(pdb_id)} and \( \(chain \w+ and resi ([\d+]+)", log
    )
    if not target:
        return {
            "error": (
                f"No pocket residues found for {pdb_id} in entry {entry_key}. "
                "Is that PDB one of the entry's apo/holo structures?"
            )
        }

    conformer_numbers = {int(n) for n in target.group(1).split("+") if n}

    residue_mapping = []
    observed = 0
    for residue in query_residues:
        present = residue["number"] in conformer_numbers
        if present:
            observed += 1
        residue_mapping.append(
            {
                "query": {"name": residue["name"], "number": residue["number"]},
                # Conformer shares the query's UniProt residue numbering, so a
                # present residue is the same name+number; null means it is
                # unobserved/missing in this structure.
                "conformer": (
                    {"name": residue["name"], "number": residue["number"]}
                    if present
                    else None
                ),
                "observed_in_conformer": present,
            }
        )

    return {
        "entry_key": entry_key,
        "query_pdb_id": entry_key.split("-", 1)[0],
        "conformer_pdb_id": pdb_id,
        "is_query_pocket": False,
        "num_query_residues": len(query_residues),
        "num_observed_in_conformer": observed,
        "num_missing_in_conformer": len(query_residues) - observed,
        "mapping_note": (
            "Each row maps a query-pocket residue to its counterpart in the "
            "conformer, matched by UniProt residue number. conformer=null "
            "(observed_in_conformer=false) means that binding residue is "
            "unobserved/missing in the conformer. Present this to the user as a "
            "query->conformer mapping, not just a flat list of residues."
        ),
        "residue_mapping": residue_mapping,
    }


@mcp.tool(title="AhojDB: Find conformers by binding-residue coverage")
def find_conformers(
    entry_key: Annotated[
        str,
        Field(
            description=(
                "Required. An AhojDB entry key from a search_apoholo result, "
                "e.g. 1a01-A-HEM-142. The entry is already specific to one "
                "binding pocket and ligand."
            ),
        ),
    ],
    min_observed_percent: Annotated[
        float,
        Field(
            description=(
                "Optional. Keep only conformers where at least this percentage of "
                "the binding-site residues are structurally observed (resolved). "
                "Default: 80.0."
            ),
            ge=0.0,
            le=100.0,
        ),
    ] = 80.0,
    state: Annotated[
        str,
        Field(
            description=(
                "Optional. Which conformers to consider: 'holo' (bind the ligand, "
                "default), 'apo' (ligand-free), or 'all'."
            ),
        ),
    ] = "holo",
    limit: Annotated[
        int,
        Field(
            description="Optional. Maximum number of conformers to return. Default: 10.",
            ge=1,
        ),
    ] = 10,
) -> dict:
    """
    Finds conformers for one AhojDB pocket that share its ligand-binding site and
    whose binding-site residues are sufficiently observed.

    Use this ONLY to drill into an entry you already found with search_apoholo.
    It requires an entry_key (e.g. 1a01-A-HEM-142), never a bare 4-character PDB
    ID.

    Conformers are filtered to those with at least min_observed_percent of their
    binding-site residues resolved, then sorted best-first (highest coverage,
    then most similar). The applied threshold is reported so it can be explained
    and adjusted.
    """
    log_call(
        "find_conformers",
        entry_key=entry_key,
        min_observed_percent=min_observed_percent,
        state=state,
        limit=limit,
    )

    entry_key = entry_key.strip()
    if not entry_key:
        return {"error": "entry_key must not be empty."}

    state = state.strip().lower()
    if state not in ("holo", "apo", "all"):
        return {"error": "state must be one of: holo, apo, all."}

    url = f"https://apoholo.cz/api/db/entry/{quote(entry_key, safe='')}/query-result"
    try:
        with urlopen(url, timeout=60) as response:
            data = load(response)
    except HTTPError as error:
        return {"error": f"Could not retrieve entry (HTTP {error.code})."}
    except URLError:
        return {"error": "Could not reach AhojDB."}

    structures = []
    if state in ("holo", "all"):
        structures += data.get("found_holo") or []
    if state in ("apo", "all"):
        structures += data.get("found_apo") or []

    matched = []
    for s in structures:
        num = s.get("mapped_binding_residues_num") or 0
        observed = s.get("mapped_binding_residues_observed") or 0
        observed_percent = round(observed / num * 100, 1) if num else 0.0
        if observed_percent < min_observed_percent:
            continue
        matched.append(
            {
                "pdb_id": s.get("pdb_id"),
                "chains": s.get("chains"),
                "assignment": s.get("apoholo_assignment"),  # A=apo, H=holo
                "observed_percent": observed_percent,
                "binding_residues_observed": observed,
                "binding_residues_total": num,
                "rmsd": s.get("rmsd"),
                "tm_score": s.get("tm_score"),
                "resolution": s.get("resolution"),
                "ligands": s.get("ligands"),
            }
        )

    matched.sort(
        key=lambda c: (c["observed_percent"], c["tm_score"] or 0), reverse=True
    )

    return {
        "entry_key": entry_key,
        "state": state,
        "applied_filter": {
            "min_observed_percent": min_observed_percent,
            "note": (
                "Results are filtered to conformers with at least "
                f"{min_observed_percent}% of their binding-site residues observed. "
                "Tell the user this threshold was applied and that it can be "
                "changed via the min_observed_percent parameter."
            ),
        },
        "num_considered": len(structures),
        "num_matched": len(matched),
        "results": matched[:limit],
        "truncated": len(matched) > limit,
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
