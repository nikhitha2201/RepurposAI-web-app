# RepurposAI - Drug & Target Discovery Studio

RepurposAI is a modular, open-source web app for drug repurposing and target discovery.  
It integrates ML, cheminformatics, and the Open Targets Platform.

## Getting Started

1. Install dependencies: `pip install -r requirements.txt`
2. Run the Streamlit app: `streamlit run app/main.py`
3. Explore modules: similarity search, target prediction, pathway mapping

Here’s a suggested structure:
```
RepurposAI/
├── README.md
├── requirements.txt
├── .gitignore
├── app/
│   ├── main.py                 # Streamlit main app
│   ├── dashboard.py            # Dashboard assembly module
│   ├── visualization.py        # Plotly/Seaborn visualization module
│   └── utils/
│       ├── api_integration.py  # Open Targets & KEGG API helpers
│       ├── similarity.py       # RDKit similarity functions
│       └── ml_model.py         # Placeholder for target prediction ML
├── data/                       # Sample data or placeholder CSVs
├── docs/                       # README, usage guide, hackathon slides
├── tests/                      # Basic test scripts for functions
└── notebooks/                  # Jupyter notebooks for testing ML/API modules
    └── example_notebook.ipynb
```
