# (Streamlit starter)
import streamlit as st
from app.dashboard import show_dashboard

st.set_page_config(page_title="RepurposAI", layout="wide")
st.title("RepurposAI - Drug & Target Discovery Studio")

# Call dashboard function (placeholder)
show_dashboard()
