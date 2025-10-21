import requests
#from collections import defaultdict
from ratelimit import limits, sleep_and_retry

def fetch_target_data(compound):
    # Placeholder function to integrate Open Targets API
    return {"gene": "TP53", "disease": "Cancer", "evidence_score": 0.95}

class KEGGAPIClient:
    """Client for interacting with the KEGG API."""
    
    BASE_URL = "https://rest.kegg.jp"
    RATE_LIMIT_CALLS = 3
    RATE_LIMIT_PERIOD = 1
    
    def __init__(self, timeout: float = 10.0):
        self.timeout = timeout
        self.session = requests.Session()

    # Set rate limit for KEGG API calls ("Please limit your API calls up to 3 times per second. Otherwise, your access will be blocked.")
    @sleep_and_retry
    @limits(calls=RATE_LIMIT_CALLS, period=RATE_LIMIT_PERIOD)
    def _make_request(self, endpoint: str) -> requests.Response:
        """Make a rate-limited request to KEGG API."""
        url = f"{self.BASE_URL}/{endpoint}"
        response = self.session.get(url, timeout=self.timeout)
        self._validate_response(response)  # Validates and raises if needed
        return response
    
    def _validate_response(self, response: requests.Response, api_name: str = "KEGG API") -> None:
        """
        Validate API response status code and raise appropriate errors.
        """
        if response.status_code == 400:
            raise ValueError(
                f"{api_name} returned Status Code 400: Bad request "
                "(syntax error, wrong database name, etc.)"
            )
        elif response.status_code == 404:
            raise ValueError(
                f"{api_name} returned Status Code 404: Not found "
                "(e.g., requesting amino acid sequence for RNA)"
            )
        elif response.status_code != 200:
            raise ValueError(
                f"{api_name} returned Status Code {response.status_code}"
            )
    
    def find_kegg_id(self, gene_target: str, organism: str = "hsa") -> str:
        """
        Get the KEGG gene ID for a given gene symbol (e.g., BRCA1, TP53).
        Ensures we pick the entry whose primary name exactly matches the symbol.
        """
        try:
            response = self._make_request(f"find/genes/{gene_target}")
            if not response.ok or not response.text.strip():
                return None

#            candidates = []
            for line in response.text.strip().splitlines():
                # Can shorten here maybe (all hsa are at the start of the queries)
                if not line.startswith(f"{organism}:"):
                    continue
                parts = line.split("\t", 1)
                if len(parts) < 2:
                    continue
                gene_id, desc = parts
                gene_names = desc.split(";")[0]  # take the first part before the semicolon
                primary_name = gene_names.split(",")[0].strip().upper()

#                # Save candidate info
#                candidates.append((gene_id, primary_name, gene_names))

                # Return if exact match
                if primary_name == gene_target.upper():
                    return gene_id  # e.g. "hsa:672"

            # fallback: partial match (if no exact)
            # CAN REMOVE THIS? RAISE ERROR IF NO EXACT MATCH TO KEGG ID
#           if candidates:
#                return candidates[0][0]
            return None
        
        except ValueError as e:
            raise ValueError(f"Failed to find KEGG id for {gene_target}: {e}")
        except requests.exceptions.Timeout:
            raise TimeoutError(f"Request timed out after {self.timeout}s")
        except requests.exceptions.RequestException as e:
            raise ConnectionError(f"Failed to connect to KEGG API: {str(e)}")

    def fetch_pathway_all(self, kegg_id: str) -> list[str]:
        """Fetch pathway codes associated with KEGG ID."""
        try:
            response = self._make_request(f"link/pathway/{kegg_id}")
            return self._parse_pathway_response(response.text)
        
        except ValueError as e:
            raise ValueError(f"Failed to fetch pathways for {kegg_id}: {e}")
        except requests.exceptions.Timeout:
            raise TimeoutError(f"Request timed out after {self.timeout}s")
        except requests.exceptions.RequestException as e:
            raise ConnectionError(f"Failed to connect to KEGG API: {str(e)}")
        
    def _parse_pathway_response(self, response_text: str) -> list[str]:
        """Parse the pathway response text into a list."""
        pathway_list = []

        for line in response_text.strip().splitlines():
            parts = line.split()
            if len(parts) == 2 and parts[1].startswith("path:"):
                path = parts[1].split("path:")[1]
                pathway_list.append(path)
        
        return pathway_list
    
if __name__ == "__main__":
    kegg_client = KEGGAPIClient()
    kegg_id = kegg_client.find_kegg_id("BRCA1")
    print(kegg_id)
    pathway_all = kegg_client.fetch_pathway_all(kegg_id)
    print(pathway_all)