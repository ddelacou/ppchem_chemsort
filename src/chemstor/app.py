import streamlit as st
import numpy as np
from typing import Dict, Any, Tuple, List, Optional, Union
import time

from chemstor.functions import (
    get_compound_safety_data,
    get_name_and_smiles,
    classify_acid_base,
    get_mp_bp,
    compound_state,
    prioritize_pictograms,
    initialize_storage_groups,
    chemsort_multiple_order_3
)

# Page config
st.set_page_config(
    page_title="ChemStorM - Chemical Storage Manager",
    page_icon="🧪",
    layout="wide"
)

# Initialize session state
if 'cache' not in st.session_state:
    st.session_state.cache = {}
if 'stored_compounds' not in st.session_state:
    st.session_state.stored_compounds = []
if 'processed_compounds' not in st.session_state:
    st.session_state.processed_compounds = []
if 'storage_groups' not in st.session_state:
    st.session_state.storage_groups = initialize_storage_groups()
if 'displayed_compounds' not in st.session_state:
    st.session_state.displayed_compounds = []

# Define a type for compound data
CompoundData = Dict[str, Any]

@st.cache_data(ttl=3600, show_spinner=False)
def cached_get_compound_safety_data(compound_name: str) -> Tuple[str, List[str], List[str]]:
    """Cache the PubChem safety data API call"""
    return get_compound_safety_data(compound_name)

@st.cache_data(ttl=3600, show_spinner=False)
def cached_get_name_and_smiles(cid: str) -> Tuple[str, str, str]:
    """Cache the PubChem name and SMILES data API call"""
    return get_name_and_smiles(cid)

@st.cache_data(ttl=3600, show_spinner=False)
def cached_get_mp_bp(compound_name: str) -> Tuple[Optional[float], Optional[float], Optional[float], Optional[float]]:
    """Cache the melting/boiling point API call"""
    return get_mp_bp(compound_name)

def create_error_compound_data(compound_name: str, error_message: str) -> CompoundData:
    """Create a compound data structure for error cases"""
    return {
        "name": compound_name,
        "iupac": "Error",
        "smiles": "Error",
        "sorted_pictograms": [],
        "hazard_statements": [error_message],
        "acid_base_class": "unknown",
        "state_room_temp": "unknown",
        "melting_point_c": None,
        "boiling_point_c": None,
        "has_error": True  # Flag to indicate this is an error entry
    }

def process_compound(compound_name: str) -> CompoundData:
    """Process a compound with progress tracking and caching"""
    progress_bar = None
    status_text = None

    try:
        # Check session cache first
        if compound_name in st.session_state.cache:
            return st.session_state.cache[compound_name]

        # Create progress indicators
        progress_bar = st.progress(0)
        status_text = st.empty()
        
        # Step 1: Safety data (30%)
        status_text.text("Fetching safety data...")
        cid, pictograms, hazards = cached_get_compound_safety_data(compound_name)
        if not cid:
            return create_error_compound_data(
                compound_name,
                f"Could not find safety data for {compound_name}"
            )
        progress_bar.progress(30)
        
        # Step 2: Chemical information (60%)
        status_text.text("Retrieving chemical information...")
        name, iupac, smiles = cached_get_name_and_smiles(cid)
        progress_bar.progress(60)
        
        # Step 3: Physical properties (90%)
        status_text.text("Processing physical properties...")
        mp_c, bp_c, mp_f, bp_f = cached_get_mp_bp(compound_name)
        progress_bar.progress(90)
        
        # Final processing
        acid_base = classify_acid_base(name, iupac, smiles, hazards)
        state = compound_state(mp_c, bp_c, mp_f, bp_f)
        sorted_picto = prioritize_pictograms(pictograms)
        
        # Clean up progress indicators
        progress_bar.progress(100)
        time.sleep(0.1)  # Short pause to show completion
        
        # Create compound data
        compound_data: CompoundData = {
            "name": compound_name,
            "iupac": iupac,
            "smiles": smiles,
            "sorted_pictograms": sorted_picto,
            "hazard_statements": hazards,
            "acid_base_class": acid_base,
            "state_room_temp": state,
            "melting_point_c": mp_c,
            "boiling_point_c": bp_c,
            "has_error": False
        }
        
        # Save to session cache
        st.session_state.cache[compound_name] = compound_data
        return compound_data
        
    except Exception as e:
        error_msg = f"Error processing {compound_name}: {str(e)}"
        return create_error_compound_data(compound_name, error_msg)
    
    finally:
        # Clean up progress indicators
        if status_text:
            status_text.empty()
        if progress_bar:
            progress_bar.empty()

def display_pictogram(picto_name: str) -> str:
    """Convert pictogram name to emoji"""
    return {
        "Explosive": "💥",
        "Flammable": "🔥",
        "Oxidizer": "⭐",
        "Compressed Gas": "💨",
        "Corrosive": "⚗️",
        "Acute Toxic": "☠️",
        "Health Hazard": "⚕️",
        "Irritant": "⚠️",
        "Environmental Hazard": "🌍"
    }.get(picto_name, "❓")

# Main layout
st.title("ChemStorM - Chemical Storage Manager")

# Sidebar
with st.sidebar:
    st.header("Add Compounds")
    
    # Compound input form
    with st.form("compound_form"):
        compound_input = st.text_input("Enter compound name:")
        submitted = st.form_submit_button("Add Compound")
        
    if submitted and compound_input:
        compound_name = compound_input.strip()
        if compound_name in st.session_state.stored_compounds:
            st.warning(f"'{compound_name}' is already in the list.")
        else:
            compound_data = process_compound(compound_name)
            if not compound_data.get("has_error", False):
                st.session_state.stored_compounds.append(compound_name)
                st.session_state.processed_compounds.append(compound_data)
                st.session_state.storage_groups = chemsort_multiple_order_3(
                    st.session_state.processed_compounds,
                    st.session_state.storage_groups
                )
                st.success(f"Added {compound_name}")
            else:
                st.error(compound_data["hazard_statements"][0])
            st.rerun()

    # Display stored compounds
    if st.session_state.stored_compounds:
        st.write("### Stored Compounds")
        for compound in st.session_state.stored_compounds:
            col1, col2 = st.columns([0.8, 0.2])
            with col1:
                st.write(f"• {compound}")
            with col2:
                if st.button("🗑️", key=f"del_{compound}"):
                    st.session_state.stored_compounds.remove(compound)
                    st.session_state.processed_compounds = [
                        c for c in st.session_state.processed_compounds 
                        if c['name'] != compound
                    ]
                    if compound in st.session_state.cache:
                        del st.session_state.cache[compound]
                    st.rerun()

    if st.button("Clear All"):
        st.session_state.stored_compounds = []
        st.session_state.processed_compounds = []
        st.session_state.storage_groups = initialize_storage_groups()
        st.session_state.displayed_compounds = []
        st.session_state.cache = {}
        st.rerun()

# Main content - Storage Groups
if st.session_state.processed_compounds:
    for group_name, states in st.session_state.storage_groups.items():
        if not any(compounds for compounds in states.values()):
            continue
            
        with st.expander(f"{group_name.replace('_', ' ').title()}", expanded=True):
            for state, compounds in states.items():
                if compounds:
                    st.subheader(f"{state.title()} State")
                    for compound in compounds:
                        col1, col2 = st.columns([0.7, 0.3])
                        with col1:
                            pictos = " ".join(
                                display_pictogram(p) 
                                for p in compound["sorted_pictograms"]
                            )
                            st.write(f"**{compound['name']}** {pictos}")
                        with col2:
                            if st.button("Details", key=f"details_{compound['name']}"):
                                if compound not in st.session_state.displayed_compounds:
                                    st.session_state.displayed_compounds.append(compound)

# Compound Details
if st.session_state.displayed_compounds:
    st.header("Compound Details")
    cols = st.columns(len(st.session_state.displayed_compounds))
    for idx, compound in enumerate(st.session_state.displayed_compounds):
        with cols[idx]:
            st.markdown(f"### {compound['name']}")
            if not compound.get("has_error", False):
                st.markdown(f"**IUPAC:** {compound['iupac']}")
                st.markdown(f"**SMILES:** {compound['smiles']}")
                st.markdown(f"**State:** {compound['state_room_temp']}")
                st.markdown(f"**Acid/Base:** {compound['acid_base_class']}")
                if compound['sorted_pictograms']:
                    st.markdown("**Hazards:**")
                    st.markdown(" ".join(
                        display_pictogram(p) 
                        for p in compound['sorted_pictograms']
                    ))
            else:
                st.error(compound['hazard_statements'][0])
            if st.button("Close", key=f"close_{compound['name']}"):
                st.session_state.displayed_compounds.remove(compound)
                st.rerun()