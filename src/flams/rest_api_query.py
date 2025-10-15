import requests

def query_uniprot_rest(ptm_keyword, taxon_id, limit=10):
    url = "https://rest.uniprot.org/uniprotkb/search"

    params = {
        "query": f"reviewed:true AND organism_id:{taxon_id} AND ft_mod_res:{ptm_keyword}",
        "fields": "accession,ft_mod_res",
        "format": "tsv",
        "size": limit
    }

    # response = requests.get(url, params=params)
    
    # if response.status_code == 200:
    #     print("Query successful! Showing results:\n")
    #     print(response.text)
    # else:
    #     print(f"Query failed. Status code: {response.status_code}")
    #     print(response.text)

if __name__ == "__main__":
    ptm = input("Enter PTM keyword (e.g. acetyl, phospho): ").strip().lower()
    taxon = input("Enter NCBI Taxonomy ID (default = 9606 for human): ").strip()
    limit = input("Enter number of results to retrieve (default = 10): ").strip()

    if not taxon:
        taxon = "9606"
    if not limit.isdigit():
        limit = 10
    else:
        limit = int(limit)

    query_uniprot_rest(ptm_keyword=ptm, taxon_id=taxon, limit=limit)