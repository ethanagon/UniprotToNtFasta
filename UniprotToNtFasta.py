from Bio import Entrez, SeqIO
from io import StringIO
import time
import requests
from requests.adapters import HTTPAdapter, Retry
import re


Entrez.email = "ethanagon22@gmail.com"
TEST = True

input_file = "idsOut_PF00520.txt"
output_file = "ntseqs_PF00520.fa"
species_map_out = "PF00520_map.smap"
species_list_out = "PF00520_species_list.txt"

if TEST:
    input_file = "test_list.txt"
    output_file = "ntseqs_small_test.fa"
    species_map_out = "small_map_test.smap"
    species_list_out = "small_sp_list_test.txt"

re_next_link = re.compile(r'<(.+)>; rel="next"')
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

UNIPROT_API_BASE = "https://rest.uniprot.org"
BATCH_SIZE = 500
UNIPROT_BATCH_SIZE = 10000

def submit_uniprot_mapping(ids, from_db="UniProtKB_AC-ID", to_db="EMBL-GenBank-DDBJ_CDS"):
    response = requests.post(
        f"{UNIPROT_API_BASE}/idmapping/run",
        data={"from": from_db, "to": to_db, "ids": ",".join(ids)}
    )
    response.raise_for_status()
    return response.json()["jobId"]

def wait_for_uniprot_job(job_id):
    secs = 0
    status_url = f"{UNIPROT_API_BASE}/idmapping/status/{job_id}"
    while True:
        response = requests.get(status_url)
        response.raise_for_status()
        data = response.json()
        if "results" in data or data.get("jobStatus") == "FINISHED":
            return  # Job done
        elif data.get("jobStatus") == "FAILED":
            raise Exception("UniProt ID mapping job failed.")
        print(f"Waiting {secs} seconds for UniProt mapping job to finish...", end="\r")
        time.sleep(2)
        secs += 2
        
def get_next_link(headers):
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)
        
def get_page(page_url):
        while page_url:
            response = session.get(page_url)
            response.raise_for_status()
            total = response.headers["x-total-results"]
            yield response, total
            page_url = get_next_link(response.headers)
        
def fetch_uniprot_results(job_id):
    print("\nFetching UniProt mapping results...")

    mapping = {}
    url = f"{UNIPROT_API_BASE}/idmapping/results/{job_id}?format=json&size=500"
    for page, total in get_page(url):
        page_json = page.json()
        for entry in page_json.get("results", []):
            from_id = entry["from"]
            to_id = entry["to"]
            mapping.setdefault(from_id, to_id)
        print(f"processed {len(mapping)} ids in batch")
    return mapping

def fetch_cds_from_ena(ENA_mapping, fmt="embl"):
    ENA_ids = list(ENA_mapping.values())
    id_str = ",".join(ENA_ids)
    nuc_to_uniprot = {v:k for k,v in ENA_mapping.items()} #mapping bijective from how we built it
    
    all_filtered_records = []
    n_batches = 1
    for batch in [ENA_ids[i:i+BATCH_SIZE] for i in range(0, len(ENA_ids), BATCH_SIZE)]:
        id_str = ",".join(batch)
        url = f"https://www.ebi.ac.uk/ena/browser/api/{fmt}/{id_str}"
        print(f"Requesting sequences for IDs {(n_batches-1) * BATCH_SIZE}-{((n_batches-1) * BATCH_SIZE) + len(batch)}")
        response = requests.get(url)
        n_batches += 1
    
        if not response.ok:
            print(f"Failed to fetch records for batch: {response.status_code}")
            

        handle = StringIO(response.text)
        records = list(SeqIO.parse(handle, fmt))
        filtered_records = [record for record in records if any(ch != 'N' for ch in record.seq)]
        for record in filtered_records:
            organism = record.annotations.get('organism')
            species_name = extract_binomial_name(organism)
            if(species_name):
                    nuc_id = record.id
                    up_id = nuc_to_uniprot.get(nuc_id, "UnknownUniProt")
                    record.id = f"{up_id}|{record.id}"
                    record.description = species_name
                    all_filtered_records.append(record)
    return all_filtered_records

def extract_binomial_name(species):
    pattern = r"\b([A-Z][a-z]+) (?!sp\b|sp\.\b)([a-z]+)\b"
    match = re.search(pattern, species.strip())
    if match:
        genus, species = match.groups()
        return f"{genus} {species}"
    return False
    
# create .smap file for Spimap
def write_species_mapping(organisms):
    with open(species_map_out, 'w') as file:
        for organism in organisms:
            file.write(f"*{organism}    {organism}\n")

# Create easily copied species list for Timetree to spit out nwk file
def write_species_list(organisms, list_out="species_list.txt"):
    with open (list_out, 'w') as file:
        for organism in organisms:
            file.write(f"{organism}\n")

#Doing this separately cuz easiest CDS from refseq is fasta only 
def get_organism_names_from_uniprot(uniprot_ids): 
    print("Getting organism names from UniProt for refseq-only sequences")
    query = " OR ".join(f"(accession:{uid})" for uid in uniprot_ids)
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": query,
        "fields": "accession,organism_name",
        "format": "json",
        "size": len(uniprot_ids)
    }

    try:
        response = requests.get(url, params=params)
        response.raise_for_status
        results = response.json()['results']

    except Exception as e:
        print(f"Batch UniProt query failed: {e}")
        return {}

    species_mapping = {}
    key_errors = []
    for result in results:
        uni_id = result.get('primaryAccession','unkownkUniprotID')
        try:
            species_name = result['organism']['scientificName']
            species_mapping.update({uni_id: species_name})
        except KeyError:
            species_mapping.update({uni_id: "unknownSpecies"})
            key_errors.append(uni_id)
    if(len(key_errors) > 0):
        print(f"Could not get species for refseq-only UniProt ids {', '.join(key_errors)}")
    return species_mapping

def fetch_cds_from_refseq(mapping):
    all_filtered_records = []
    flat_nuc_ids = []
    nuc_to_uniprot = {}
    
    
    for up_id, nuc in mapping.items():
        nuc_base = nuc.split('.')[0]  # Trim version
        flat_nuc_ids.append(nuc_base)
        nuc_to_uniprot[nuc_base] = up_id
    n_batches = 1
    for batch in [flat_nuc_ids[i:i+BATCH_SIZE] for i in range(0, len(flat_nuc_ids), BATCH_SIZE)]:
        print(f"Requesting sequences for IDs {(n_batches-1) * BATCH_SIZE}-{((n_batches-1) * BATCH_SIZE) + len(batch)}")
        try:
            handle = Entrez.efetch(db="nucleotide", id=batch, rettype="fasta_cds_na", retmode="text")
            records = list(SeqIO.parse(handle, "fasta"))

        except Exception as e:
            print(f"Error fetching FASTA for batch {batch}: {e}")
            
        filtered_records = [record for record in records if any(ch != 'N' for ch in record.seq)]
        uniprots = [nuc_to_uniprot[nuc_base] for nuc_base in batch]
        accession_to_organism = get_organism_names_from_uniprot(uniprots)
        for record in filtered_records:
            messy_id = record.id
            nuc_id =  messy_id[messy_id.find('|')+1 : messy_id.find('.')]
            species_name = accession_to_organism[nuc_to_uniprot[nuc_id]]
            if(species_name and species_name != "unknownSpecies"):
                    up_id = nuc_to_uniprot.get(nuc_id, "UnknownUniProt")
                    record.id = f"{up_id}|{nuc_id}"
                    record.description = species_name
                    all_filtered_records.append(record)
        time.sleep(0.4)
        n_batches += 1
    return all_filtered_records

def get_uniprot_mapping(ids, from_db="UniProtKB_AC-ID", to_db="EMBL-GenBank-DDBJ_CDS"):
    job_id = submit_uniprot_mapping(ids, from_db=from_db, to_db=to_db)
    wait_for_uniprot_job(job_id)
    mapping = fetch_uniprot_results(job_id)

    if not mapping:
        print("No mappings found.")
        return
    
    return mapping


def main():
    with open(input_file) as f:
        uniprot_ids = [line.strip() for line in f if line.strip()]
        
    ENA_mapping = {}        
    n_ENA_batches = 1
    for batch in [uniprot_ids[i:i+UNIPROT_BATCH_SIZE] for i in range(0,len(uniprot_ids), UNIPROT_BATCH_SIZE)]:
        ENA_mapping.update(get_uniprot_mapping(batch))
        print(f"Finished {n_ENA_batches} batches of Uniprot IDs")
        print(f"Submitted {min(n_ENA_batches*UNIPROT_BATCH_SIZE, len(uniprot_ids))} total Uniprot IDs for mapping")
        print(f"Mapped {len(ENA_mapping.values())} Uniprot IDs to ENA ids")
        n_ENA_batches += 1

    
    
    ENA_records = fetch_cds_from_ena(ENA_mapping)
    
    remaining_ids = list(set(uniprot_ids) - set(ENA_mapping.keys())) 
    print(f"{len(remaining_ids)} Uniprot ids not found in ENA. Querying Refseq.")

    #Write a list of what didn't come through in ENA, because it's refseq that keeps fucking up
    write_species_list(remaining_ids, list_out="ENA_missed_ids.txt")
    #Write ENA results to separate file for my sanity and minimizing reruns
    refseq_mapping = {} 
    n_refseq_batches = 1   
    for batch in [remaining_ids[i:i+UNIPROT_BATCH_SIZE] for i in range(0,len(remaining_ids), UNIPROT_BATCH_SIZE)]:
        refseq_mapping.update(get_uniprot_mapping(batch, to_db="RefSeq_Nucleotide"))
        print(f"Finished {n_refseq_batches} batches of Uniprot IDs")
        print(f"Submitted {min(n_refseq_batches*UNIPROT_BATCH_SIZE, len(remaining_ids))} total Uniprot IDs for mapping")
        print(f"Mapped {len(refseq_mapping.values())} Uniprot IDs to RefSeq ids")
        n_refseq_batches += 1

    
    print("Fetching CDS records from refseq")
    refseq_records = fetch_cds_from_refseq(refseq_mapping)
    all_records = refseq_records + ENA_records 

    organisms = {record.description for record in all_records}


    print(f"writing species mapping to {species_map_out}")
    write_species_mapping(organisms)
    print(f"writing species list to {species_list_out}")
    write_species_list(organisms, list_out=species_list_out)

    SeqIO.write(all_records, output_file, "fasta")
    print(f"Wrote {len(all_records)} sequences to {output_file}")

if __name__ == "__main__":
    main()
