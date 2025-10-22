#/usr/bin/env python3

import requests
from typing import List, Optional, Dict
from ratelimit import limits, sleep_and_retry
import ml_model  

def fetch_target_data(compound: str) -> Dict[str, object]:
    targets = ml_model.predict_targets(compound)
    # targets should be a dictionary, for example: {"gene": "TP53"}
    return targets 

class OpenTargetsClient:
    """Client for querying the Open Targets GraphQL API."""

    BASE_URL = "https://api.platform.opentargets.org/api/v4/graphql"

    def __init__(self, timeout: float = 15.0):
        self.timeout = timeout
        self.session = requests.Session()

    def run_query(self, query: str, variables: Optional[dict] = None) -> dict:
        """Run a GraphQL query and return JSON results."""
        response = self.session.post(
            self.BASE_URL,
            json={"query": query, "variables": variables},
            timeout=self.timeout,
        )
        if response.status_code != 200:
            raise ValueError(
                f"Open Targets API error {response.status_code}: {response.text}"
            )
        return response.json()

    def get_target_disease_evidence(self, gene_ensembl_id: str, disease_efo_id: str) -> List[dict]:
        """Retrieve evidence linking a target gene and a disease."""
        query = """
        query targetDiseaseEvidence($geneId: String!, $diseaseId: String!) {
          disease(efoId: $diseaseId) {
            id
            name
            evidences(
              datasourceIds: ["intogen"],
              ensemblIds: [$geneId]
            ) {
              count
              rows {
                disease { id name }
                diseaseFromSource
                target { id approvedSymbol }
                mutatedSamples {
                  functionalConsequence { id label }
                  numberSamplesTested
                  numberMutatedSamples
                }
                resourceScore
                significantDriverMethods
                cohortShortName
              }
            }
          }
        }
        """
        variables = {"geneId": gene_ensembl_id, "diseaseId": disease_efo_id}
        data = self.run_query(query, variables)
        evidences = (
            data.get("data", {})
            .get("disease", {})
            .get("evidences", {})
            .get("rows", [])
        )
        return evidences

class KEGGAPIClient:
    """Client for interacting with the KEGG REST API."""

    BASE_URL = "https://rest.kegg.jp"
    RATE_LIMIT_CALLS = 3
    RATE_LIMIT_PERIOD = 1

    def __init__(self, timeout: float = 10.0):
        self.timeout = timeout
        self.session = requests.Session()

    # Rate-limited requests: KEGG allows 3 requests/sec
    @sleep_and_retry
    @limits(calls=RATE_LIMIT_CALLS, period=RATE_LIMIT_PERIOD)
    def _make_request(self, endpoint: str) -> requests.Response:
        """Make a rate-limited request to the KEGG API."""
        url = f"{self.BASE_URL}/{endpoint}"
        response = self.session.get(url, timeout=self.timeout)
        self._validate_response(response)
        return response

    def _validate_response(self, response: requests.Response, api_name: str = "KEGG API") -> None:
        """Validate API response and raise errors if needed."""
        if response.status_code == 400:
            raise ValueError(f"{api_name} 400: Bad request (syntax error or wrong database).")
        elif response.status_code == 404:
            raise ValueError(f"{api_name} 404: Not found (invalid entry or type).")
        elif response.status_code != 200:
            raise ValueError(f"{api_name} returned unexpected status {response.status_code}.")

    def find_kegg_id(self, gene_target: str, organism: str = "hsa") -> Optional[str]:
        """
        Get the KEGG gene ID for a given gene symbol (e.g., BRCA1, TP53).
        Ensures we pick the entry whose primary name exactly matches the symbol.
        """
        try:
            response = self._make_request(f"find/genes/{gene_target}")
            if not response.ok or not response.text.strip():
                return None

            for line in response.text.strip().splitlines():
                if not line.startswith(f"{organism}:"):
                    continue
                parts = line.split("\t", 1)
                if len(parts) < 2:
                    continue
                gene_id, desc = parts
                gene_names = desc.split(";")[0]
                primary_name = gene_names.split(",")[0].strip().upper()

                if primary_name == gene_target.upper():
                    return gene_id  # e.g. "hsa:672"
            return None

        except ValueError as e:
            raise ValueError(f"Failed to find KEGG ID for {gene_target}: {e}")
        except requests.exceptions.Timeout:
            raise TimeoutError(f"Request timed out after {self.timeout}s")
        except requests.exceptions.RequestException as e:
            raise ConnectionError(f"Failed to connect to KEGG API: {str(e)}")

    def fetch_pathway_all(self, kegg_id: str) -> List[str]:
        """Fetch all pathway codes associated with a KEGG ID."""
        try:
            response = self._make_request(f"link/pathway/{kegg_id}")
            return self._parse_pathway_response(response.text)
        except ValueError as e:
            raise ValueError(f"Failed to fetch pathways for {kegg_id}: {e}")
        except requests.exceptions.Timeout:
            raise TimeoutError(f"Request timed out after {self.timeout}s")
        except requests.exceptions.RequestException as e:
            raise ConnectionError(f"Failed to connect to KEGG API: {str(e)}")

    def _parse_pathway_response(self, response_text: str) -> List[str]:
        """Parse the KEGG pathway link response text into a list of pathway codes."""
        pathway_list = []
        for line in response_text.strip().splitlines():
            parts = line.split()
            if len(parts) == 2 and parts[1].startswith("path:"):
                path = parts[1].split("path:")[1]
                pathway_list.append(path)
        return pathway_list


if __name__ == "__main__":
    compound = "aspirin"
    predicted_data = fetch_target_data(compound)
    # KEGG query example
    kegg_client = KEGGAPIClient()
    gene_symbol = predicted_data["gene"]

    kegg_id = kegg_client.find_kegg_id(gene_symbol)
    print(f"KEGG ID: {kegg_id}")

    if kegg_id:
        pathways = kegg_client.fetch_pathway_all(kegg_id)
        print(f"Pathways for {gene_symbol}: {pathways}")
    else:
        print("No KEGG ID found for this gene.")
