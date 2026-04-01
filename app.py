"""
NeBULA - Bioisostere Replacement Web App
For GM sir - Access from anywhere!
"""
import streamlit as st
import pandas as pd
import os
import sys
from datetime import datetime
import requests

# Add NeBULA to path
sys.path.insert(0, os.path.dirname(__file__))

st.set_page_config(page_title="NeBULA - Bioisostere Replacement", page_icon="🧬", layout="wide")

st.markdown("""
<style>
    .stApp { background: linear-gradient(135deg, #0d1b2a 0%, #1b263b 100%); }
    .main-header { color: #00d4ff; font-size: 2.5rem; font-weight: bold; }
    .result-card { background: rgba(0,212,255,0.1); border-radius: 15px; padding: 1.5rem; margin: 0.75rem 0; border: 1px solid #00d4ff; }
    .molecule-box { background: #1b263b; border-radius: 10px; padding: 1rem; margin: 0.5rem 0; font-family: monospace; word-break: break-all; }
    .success-box { background: rgba(0,255,136,0.15); border-left: 4px solid #00ff88; padding: 1rem; border-radius: 0 10px 10px 0; }
    .info-box { background: rgba(0,212,255,0.1); border-left: 4px solid #00d4ff; padding: 1rem; border-radius: 0 10px 10px 0; }
</style>
""", unsafe_allow_html=True)

st.markdown('<h1 class="main-header">🧬 NeBULA</h1>', unsafe_allow_html=True)
st.markdown('<p style="color:#88a4c4;">Next-Generation Bioisostere Utility for Lead Optimization</p>', unsafe_allow_html=True)
st.markdown("---")

# Import NeBULA functions
try:
    from reaction import count_heavy_atoms, get_heavy_atom_difference
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
    RDKIT_AVAILABLE = True
except Exception as e:
    RDKIT_AVAILABLE = False
    st.error(f"RDKit not available: {e}")

def get_pubchem_image(smiles, size=300):
    """Get 2D structure image from PubChem"""
    if not smiles:
        return None
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/PNG?image_size={size}"
        resp = requests.get(url, timeout=10)
        if resp.status_code == 200:
            return resp.content
    except:
        pass
    return None

def run_nebula(smiles_input):
    """Run NeBULA bioisostere replacement"""
    if not RDKIT_AVAILABLE:
        return None, "RDKit not available"
    
    try:
        # Load reactions
        reaction_file = os.path.join(os.path.dirname(__file__), 'Reaction_Fsp3-rich.csv')
        df_reactions = pd.read_csv(reaction_file, encoding='latin-1')
        
        # Parse input
        mol = Chem.MolFromSmiles(smiles_input)
        if mol is None:
            return None, "Invalid SMILES"
        
        Chem.SanitizeMol(mol)
        input_canonical = Chem.MolToSmiles(mol)
        input_heavy_atoms = count_heavy_atoms(mol)
        
        results = []
        seen_products = set()
        
        for idx, row in df_reactions.iterrows():
            smarts_reaction = row['smarts_reaction']
            try:
                reaction = AllChem.ReactionFromSmarts(smarts_reaction)
                reaction.Initialize()
                products = reaction.RunReactants((mol,))
            except:
                continue
            
            for prod_set in products:
                for prod in prod_set:
                    try:
                        prod_smiles = Chem.MolToSmiles(prod, isomericSmiles=True)
                        if not prod_smiles:
                            continue
                        try:
                            Chem.Kekulize(prod, clearAromaticFlags=True)
                            Chem.SanitizeMol(prod)
                        except:
                            continue
                        
                        prod_heavy_atoms = count_heavy_atoms(prod)
                        diff = get_heavy_atom_difference(smarts_reaction)
                        if abs(prod_heavy_atoms - input_heavy_atoms) <= abs(diff) and prod_smiles not in seen_products:
                            seen_products.add(prod_smiles)
                            results.append({
                                'Product_SMILES': prod_smiles,
                                'Reaction_SMARTS': smarts_reaction
                            })
                    except:
                        continue
        
        return results, None
        
    except Exception as e:
        return None, str(e)

# Sidebar
with st.sidebar:
    st.markdown("### 📊 About NeBULA")
    st.markdown("""
    **Bioisostere Replacement Tool**
    
    - 18,736+ reactions
    - Fsp3-rich transformations
    - Lead optimization
    - Drug design
    """)
    
    st.markdown("---")
    st.markdown("### 📖 How to Use")
    st.markdown("""
    1. Enter a SMILES string
    2. Click **Generate**
    3. View bioisostere replacements
    4. Click structures to enlarge
    """)
    
    st.markdown("---")
    st.markdown("### 🔬 Example SMILES")
    examples = [
        "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
        "c1ccccc1",  # Benzene
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O",  # Ibuprofen
        "CC(=O)Nc1ccc(cc1)O",  # Phenacetin
    ]
    for ex in examples:
        if st.button(ex[:30]+"...", key=f"ex_{ex}"):
            st.session_state.example_smiles = ex

# Main content
col1, col2 = st.columns([3, 1])

with col1:
    smiles_input = st.text_area(
        "🧪 Enter SMILES:",
        value=st.session_state.get('example_smiles', 'CC(=O)Oc1ccccc1C(=O)O'),
        height=100,
        placeholder="e.g., CC(=O)Oc1ccccc1C(=O)O"
    )

with col2:
    st.write("")
    st.write("")
    generate = st.button("🔄 Generate Bioisosteres", use_container_width=True)

st.markdown("---")

if generate and smiles_input:
    with st.spinner("Generating bioisostere replacements..."):
        results, error = run_nebula(smiles_input)
    
    if error:
        st.error(f"Error: {error}")
    elif results:
        st.success(f"✅ Generated {len(results)} bioisostere replacements!")
        
        # Show original
        st.markdown("### 📦 Original Molecule")
        col_orig_img, col_orig_info = st.columns([1, 2])
        
        with col_orig_img:
            img = get_pubchem_image(smiles_input, 300)
            if img:
                st.image(img, caption="Original Structure", width=250)
        
        with col_orig_info:
            mol = Chem.MolFromSmiles(smiles_input)
            if mol:
                mw = Descriptors.MolWt(mol)
                heavy_atoms = count_heavy_atoms(mol)
                st.markdown(f"""
                <div class="info-box">
                    <p><strong>Canonical SMILES:</strong></p>
                    <code class="molecule-box">{Chem.MolToSmiles(mol)}</code>
                    <p><strong>MW:</strong> {mw:.2f} | <strong>Heavy Atoms:</strong> {heavy_atoms}</p>
                </div>
                """, unsafe_allow_html=True)
        
        st.markdown("---")
        st.markdown("### 🧬 Bioisostere Replacements")
        
        # Show results in batches
        batch_size = 6
        for i in range(0, min(len(results), 30), batch_size):
            batch = results[i:i+batch_size]
            cols = st.columns(3)
            
            for j, result in enumerate(batch):
                with cols[j % 3]:
                    prod_smiles = result['Product_SMILES']
                    
                    # Get structure image
                    prod_img = get_pubchem_image(prod_smiles, 250)
                    
                    with st.container():
                        st.markdown('<div class="result-card">', unsafe_allow_html=True)
                        
                        if prod_img:
                            st.image(prod_img, caption=f"Product {i+j+1}", width=200)
                        
                        # Product info
                        prod_mol = Chem.MolFromSmiles(prod_smiles)
                        if prod_mol:
                            prod_mw = Descriptors.MolWt(prod_mol)
                            prod_heavy = count_heavy_atoms(prod_mol)
                            
                            st.markdown(f"""
                            <div style="background:#0d1b2a;padding:0.5rem;border-radius:5px;margin-top:0.5rem;">
                                <p style="margin:0;"><small>MW: {prod_mw:.2f} | Heavy Atoms: {prod_heavy}</small></p>
                            </div>
                            """, unsafe_allow_html=True)
                        
                        # SMILES
                        with st.expander("View SMILES"):
                            st.code(prod_smiles, language=None)
                        
                        st.markdown('</div>', unsafe_allow_html=True)
        
        if len(results) > 30:
            st.info(f"Showing 30 of {len(results)} results. For full results, download below.")
        
        # Download CSV
        df_results = pd.DataFrame(results)
        csv = df_results.to_csv(index=False)
        st.download_button(
            "📥 Download All Results (CSV)",
            csv,
            "nebula_bioisosteres.csv",
            "text/csv"
        )
    else:
        st.warning("No bioisostere replacements generated.")

else:
    st.info("👈 Enter a SMILES string and click Generate to start!")
    
    # Show example
    st.markdown("### 💡 Example")
    st.code("CC(=O)Oc1ccccc1C(=O)O  # Aspirin", language=None)

st.markdown("---")
st.markdown("""
<div style="text-align:center; color:#88a4c4; padding:1rem;">
    <p>NeBULA - Bioisostere Replacement Tool | For GM sir</p>
    <p style="font-size:0.8rem;">Based on: Huang et al. (2025) Medicine in Drug Discovery</p>
</div>
""", unsafe_allow_html=True)
