import streamlit as st
import re
from typing import List, Dict, Any, Tuple
import pandas as pd
import numpy as np
from chemstor import (
    get_compound_safety_data,
    get_name_and_smiles,
    classify_acid_base,
    get_mp_bp,
    compound_state,
    prioritize_pictograms,
    initialize_storage_groups,
    chemsort_multiple_order_3
)

# Set page configuration
st.set_page_config(
    page_title="ChemStor - Chemical Storage Manager",
    page_icon="🧪",
    layout="wide"
)

# Initialize session state variables
if 'stored_compounds' not in st.session_state:
    st.session_state.stored_compounds = []
if 'compound_input' not in st.session_state:
    st.session_state.compound_input = ''
if 'processed_compounds' not in st.session_state:
    st.session_state.processed_compounds = []
if 'storage_groups' not in st.session_state:
    st.session_state.storage_groups = initialize_storage_groups()
if 'selected_compound' not in st.session_state:
    st.session_state.selected_compound = None

# --- CSS Styling ---
st.markdown("""
<style>
    .storage-grid {
        display: grid;
        grid-template-columns: repeat(3, 1fr);
        gap: 10px;
    }
    .compound-button {
        margin: 2px 0;
        width: 100%;
    }
    .state-section {
        margin-top: 8px;
    }
    .main-header {
        color: #FFFFFF;
        text-align: center;
    }
    .storage-group-header {
        background-color: #000000;
        padding: 5px;
        border-radius: 5px;
    }
    .compound-list {
        max-height: 200px;
        overflow-y: auto;
    }
</style>
""", unsafe_allow_html=True)

# --- Functions ---
def add_compound():
    """Add a compound to the list of stored compounds."""
    compound = st.session_state.compound_input.strip()
    if compound and compound not in st.session_state.stored_compounds:
        st.session_state.stored_compounds.append(compound)
        st.session_state.compound_input = ''
        st.success(f"Added compound: {compound}")
    elif not compound:
        st.warning("Please enter a compound name.")
    else:
        st.info(f"Compound '{compound}' is already in the list.")

def remove_compound(compound):
    """Remove a compound from the list of stored compounds."""
    if compound in st.session_state.stored_compounds:
        st.session_state.stored_compounds.remove(compound)
        st.success(f"Removed compound: {compound}")

def display_pictogram(picto_name):
    """Return emoji representation of pictogram."""
    picto_map = {
        "Explosive": "💥",
        "Flammable": "🔥",
        "Oxidizer": "🔆",
        "Compressed Gas": "🧪",
        "Corrosive": "🧫",
        "Acute Toxic": "☠️",
        "Health Hazard": "⚠️",
        "Irritant": "❗",
        "Environmental Hazard": "🌍"
    }
    return picto_map.get(picto_name, "❓")

def process_compounds(compounds: List[str]):
    """Process the compounds and generate storage recommendations."""
    if not compounds:
        st.warning("No compounds to process.")
        return
    
    results = []
    with st.spinner("Processing compounds..."):
        for compound in compounds:
            # Display progress
            st.write(f"Processing {compound}...")
            
            # Get compound data using functions from chemstor
            cid, pictos, hazards = get_compound_safety_data(compound)
            
            if not cid:
                st.error(f"Could not find data for {compound}. Skipping.")
                continue
                
            name, iupac, smiles = get_name_and_smiles(cid)
            acid_base_class = classify_acid_base(name, iupac, smiles, hazards)
            mp_c, bp_c, mp_f, bp_f = get_mp_bp(compound)
            state = compound_state(mp_c, bp_c, mp_f, bp_f)
            sorted_picto = prioritize_pictograms(pictos)
            
            # Store all relevant information
            compound_info = {
                "name": compound,
                "iupac": iupac,
                "smiles": smiles,
                "sorted_pictograms": sorted_picto,
                "hazard_statements": hazards,
                "acid_base_class": acid_base_class,
                "state_room_temp": state,
                "melting_point_c": mp_c,
                "boiling_point_c": bp_c
            }
            results.append(compound_info)
    
    # Apply the storage grouping algorithm
    if results:
        storage_groups = initialize_storage_groups()
        sorted_storage = chemsort_multiple_order_3(results, storage_groups)
        st.session_state.processed_compounds = results
        st.session_state.storage_groups = sorted_storage
        return True
    else:
        st.error("Failed to process compounds. Please check your input.")
        return False

def display_compound_details(compound):
    """Display detailed information about a selected compound."""
    st.subheader(f"Details for {compound['name']}")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("### Basic Information")
        st.write(f"**Name:** {compound['name']}")
        st.write(f"**IUPAC:** {compound['iupac']}")
        st.write(f"**SMILES:** {compound['smiles']}")
        st.write(f"**State at Room Temperature:** {compound['state_room_temp']}")
        
        if compound['melting_point_c'] is not None:
            st.write(f"**Melting Point:** {compound['melting_point_c']:.2f}°C")
        else:
            st.write("**Melting Point:** Not available")
            
        if compound['boiling_point_c'] is not None:
            st.write(f"**Boiling Point:** {compound['boiling_point_c']:.2f}°C")
        else:
            st.write("**Boiling Point:** Not available")
            
    with col2:
        st.markdown("### Safety Information")
        
        # Display acid/base classification
        if isinstance(compound['acid_base_class'], tuple):
            st.write(f"**Acid/Base Classification:** {', '.join(compound['acid_base_class'])}")
        else:
            st.write(f"**Acid/Base Classification:** {compound['acid_base_class']}")
        
        # Display pictograms
        if compound['sorted_pictograms']:
            st.markdown("**Pictograms:**")
            for picto in compound['sorted_pictograms']:
                st.write(f"{display_pictogram(picto)} {picto}")
        else:
            st.write("**Pictograms:** None")
        
        # Display hazard statements
        if compound['hazard_statements']:
            st.markdown("**Hazard Statements:**")
            for hazard in compound['hazard_statements']:
                st.write(f"- {hazard}")
        else:
            st.write("**Hazard Statements:** None")

# --- UI Layout ---
st.markdown("<h1 class='main-header'>ChemStor - Chemical Storage Manager</h1>", unsafe_allow_html=True)

# Sidebar for input compounds
with st.sidebar:
    st.header("Input Compounds")
    st.text_input("Enter Compound Name", key='compound_input', on_change=add_compound)
    st.button("Add Compound", on_click=add_compound)
    
    # List of entered compounds
    st.markdown("### Stored Compounds")
    if st.session_state.stored_compounds:
        for compound in st.session_state.stored_compounds:
            col1, col2 = st.columns([0.8, 0.2])
            with col1:
                st.write(compound)
            with col2:
                if st.button("❌", key=f"remove_{compound}"):
                    remove_compound(compound)
                    st.rerun()
    else:
        st.write("No compounds added yet.")
    
    # Process button
    if st.button("Process Compounds"):
        if process_compounds(st.session_state.stored_compounds):
            st.success("Compounds processed successfully!")
        
    # Clear all button
    if st.button("Clear All"):
        st.session_state.stored_compounds = []
        st.session_state.processed_compounds = []
        st.session_state.storage_groups = initialize_storage_groups()
        st.session_state.selected_compound = None
        st.success("All data cleared!")

# Main content area - Show results if compounds are processed
if st.session_state.processed_compounds:
    # Display storage groups in a grid
    st.header("Storage Groups")
    
    # Filter non-empty storage groups
    non_empty_groups = {}
    for group_name, states in st.session_state.storage_groups.items():
        has_compounds = False
        for state_key, compounds in states.items():
            if compounds:
                has_compounds = True
                break
        if has_compounds:
            non_empty_groups[group_name] = states
    
    # If no compounds were sorted, display a message
    if not non_empty_groups:
        st.info("No compounds have been sorted into storage groups yet.")
    else:
        # Create a grid for storage groups
        cols_per_row = 3
        rows = np.ceil(len(non_empty_groups) / cols_per_row).astype(int)
        
        for row in range(rows):
            cols = st.columns(cols_per_row)
            for col_idx in range(cols_per_row):
                group_idx = row * cols_per_row + col_idx
                if group_idx < len(non_empty_groups):
                    group_name = list(non_empty_groups.keys())[group_idx]
                    states = non_empty_groups[group_name]
                    
                    with cols[col_idx]:
                        st.markdown(f"<div class='storage-group-header'><h3>{group_name.replace('_', ' ').title()}</h3></div>", unsafe_allow_html=True)
                        
                        # Display each state (solid, liquid, gas) if it has compounds
                        for state_key, compounds in states.items():
                            if compounds:
                                with st.expander(f"{state_key.capitalize()} ({len(compounds)})"):
                                    for compound in compounds:
                                        # Create a button for each compound
                                        picto_emoji = "".join([display_pictogram(p) for p in compound["sorted_pictograms"][:2]])
                                        button_label = f"{picto_emoji} {compound['name']}"
                                        if st.button(button_label, key=f"btn_{group_name}_{state_key}_{compound['name']}"):
                                            st.session_state.selected_compound = compound
                                            st.rerun()

    # Display selected compound details
    if st.session_state.selected_compound:
        st.markdown("---")
        display_compound_details(st.session_state.selected_compound)
else:
    # Show instructions if no compounds are processed yet
    st.info("Please add compounds in the sidebar and click 'Process Compounds' to see storage recommendations.")
    
    # Example section
    with st.expander("Examples of compounds you can try"):
        st.markdown("""
        - Acetone
        - Hydrogen peroxide
        - Sodium hydroxide
        - Hydrochloric acid
        - Methanol
        - Ethanol
        - Nitric acid
        - Benzene
        - Potassium permanganate
        - Toluene
        """)