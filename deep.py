import streamlit as st
from Bio.Seq import Seq
import random
import plotly.express as px
from collections import Counter

# --- App Branding ---
st.set_page_config(page_title="🧬 HelixSleuth", layout="wide", page_icon="🧬")

# Custom CSS for professional look
st.markdown("""
<style>
    .header { color: #2E86AB; font-family: 'Helvetica Neue', sans-serif; }
    .subheader { color: #5C6BC0; font-size: 1.1rem; }
    .stButton>button { background-color: #4CAF50; color: white; }
    .stTextInput>div>div>input { font-family: 'Courier New', monospace; }
    .footer { font-size: 0.8rem; color: #666; margin-top: 2rem; }
    .dna-display { font-family: 'Courier New', monospace; letter-spacing: 1px; }
    .highlight { background-color: #FFF59D; padding: 0 2px; }
</style>
""", unsafe_allow_html=True)

# --- Header ---
st.title("🧬 HelixSleuth")
st.markdown("""
<p class="subheader">Precision genomic analysis for research and diagnostics</p>
""", unsafe_allow_html=True)
st.divider()

# --- Sidebar (User Inputs) ---
with st.sidebar:
    st.header("⚙️ Simulation Parameters")
    dna_input = st.text_area("Enter DNA sequence (ATGC only):", "ATGCGTACGTACGTACGT", height=100).upper().strip()
    
    if st.button("Generate Random DNA"):
        dna_input = ''.join(random.choices('ATGC', k=random.randint(50, 100)))
    
    mutation_type = st.radio("Mutation Type:", 
                           ["Point", "Insertion", "Deletion"],
                           index=0,
                           help="Point: Single base change | Insertion: Add bases | Deletion: Remove bases")

# --- Core Functions (Optimized with Caching) ---
@st.cache_data
def validate_dna(sequence):
    return all(base in 'ATGC' for base in sequence)

@st.cache_data
def translate_dna(dna):
    return str(Seq(dna).translate(table=1, to_stop=True))

def apply_mutation(seq, mut_type, pos, change=None):
    if mut_type == "Point":
        return seq[:pos] + change + seq[pos+1:]
    elif mut_type == "Insertion":
        return seq[:pos] + change + seq[pos:]
    elif mut_type == "Deletion":
        return seq[:pos] + seq[pos+change:]

# --- Main Display ---
col1, col2 = st.columns([1, 1])

with col1:
    st.subheader("🔬 Original Sequence")
    if not validate_dna(dna_input):
        st.error("Invalid DNA sequence! Only A, T, G, C allowed.")
        st.stop()
    
    st.code(dna_input, language='text')
    
    if len(dna_input) >= 3:
        with st.expander("View Protein Translation"):
            protein = translate_dna(dna_input)
            st.write(f"Length: {len(protein)} amino acids")
            st.code(protein)

with col2:
    st.subheader("🧪 Mutation Controls")
    
    if mutation_type == "Point":
        pos = st.number_input("Position (0-based):", 0, len(dna_input)-1, 0)
        new_base = st.selectbox("New base:", ["A", "T", "G", "C"])
        mutated = apply_mutation(dna_input, "Point", pos, new_base)
        
    elif mutation_type == "Insertion":
        pos = st.number_input("Insert at position:", 0, len(dna_input), 0)
        insert_seq = st.text_input("Bases to insert:", "AT").upper()
        if validate_dna(insert_seq):
            mutated = apply_mutation(dna_input, "Insertion", pos, insert_seq)
        
    elif mutation_type == "Deletion":
        pos = st.number_input("Delete from position:", 0, len(dna_input)-1, 0)
        del_len = st.number_input("Delete length:", 1, len(dna_input)-pos, 1)
        mutated = apply_mutation(dna_input, "Deletion", pos, del_len)

# --- Analysis & Visualization ---
if 'mutated' in locals():
    st.divider()
    st.subheader("📊 Mutation Analysis")
    
    # Mutation Highlighting
    st.markdown("**DNA Comparison:**")
    col_a, col_b = st.columns(2)
    with col_a:
        st.markdown(f"Original ({len(dna_input)} bp)")
        st.code(dna_input)
    with col_b:
        st.markdown(f"Mutated ({len(mutated)} bp)")
        st.code(mutated)
    
    # Protein Impact
    if len(dna_input) >= 3:
        with st.expander("Protein Impact Analysis", expanded=True):
            orig_protein = translate_dna(dna_input)
            mut_protein = translate_dna(mutated)
            
            st.markdown(f"**Original Protein ({len(orig_protein)} aa):**")
            st.code(orig_protein)
            
            st.markdown(f"**Mutated Protein ({len(mut_protein)} aa):**")
            st.code(mut_protein)
            
            # Mutation statistics
            changes = sum(1 for a, b in zip(orig_protein, mut_protein) if a != b)
            st.metric("Amino Acid Changes", changes)
            
            # Visualization
            aa_counts = Counter(mut_protein)
            fig = px.bar(x=list(aa_counts.keys()), y=list(aa_counts.values()),
                         labels={'x': 'Amino Acid', 'y': 'Count'},
                         title="Mutated Protein Composition")
            st.plotly_chart(fig, use_container_width=True)

# --- Footer ---
st.divider()
st.markdown("""
<p class="footer">
    Developed with ❤️ by [Your Name] | <a href="mailto:your.email@example.com">Contact</a> | 
    Powered by Biopython & Streamlit
</p>
""", unsafe_allow_html=True)
