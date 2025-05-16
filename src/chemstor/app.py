import streamlit as st
from typing import List, Dict, Any
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



# Initialize session state
if 'stored_compounds' not in st.session_state:
    st.session_state.stored_compounds = []
if 'compound_input' not in st.session_state:
    st.session_state.compound_input = ''
if 'processing_result' not in st.session_state:
    st.session_state.processing_result = []
if 'storage_groups' not in st.session_state:
    st.session_state.storage_groups = initialize_storage_groups()

# --- Compound Functions ---
def add_compound():
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
    if compound in st.session_state.stored_compounds:
        st.session_state.stored_compounds.remove(compound)
        st.success(f"Removed compound: {compound}")

def process_compounds(compounds: List[str]):
    storage_groups = st.session_state.storage_groups
    results = []

    for compound in compounds:
        cid, pictos, hazards = get_compound_safety_data(compound)
        name, iupac, smiles = get_name_and_smiles(cid)
        acid_base_class = classify_acid_base(name, iupac, smiles, hazards)
        mp_c, bp_c, mp_f, bp_f = get_mp_bp(compound)
        state = compound_state(mp_c, bp_c, mp_f, bp_f)
        sorted_picto = prioritize_pictograms(pictos)

        compound_info = {
            "name": compound,
            "sorted_pictograms": sorted_picto,
            "hazard_statements": hazards,
            "acid_base_class": acid_base_class,
            "state_room_temp": state
        }
        results.append(compound_info)

    sorted_storage = chemsort_multiple_order_3(results, storage_groups)
    st.session_state.storage_groups = sorted_storage
    return sorted_storage

# --- UI Layout ---
st.set_page_config(layout="wide")
st.title("Chemical Compound Sorting Interface")

# Sidebar Input
with st.sidebar:
    st.header("Input Compounds")
    st.text_input("Enter Compound Name", key='compound_input', on_change=add_compound)
    st.button("Add Compound", on_click=add_compound)

    st.markdown("### Stored Compounds")
    if st.session_state.stored_compounds:
        for compound in st.session_state.stored_compounds:
            col1, col2 = st.columns([0.85, 0.15])
            with col1:
                st.write(compound)
            with col2:
                if st.button("❌", key=f"remove_{compound}"):
                    remove_compound(compound)
                    st.rerun()
    else:
        st.write("No compounds added yet.")

    if st.button("Apply"):
        if st.session_state.stored_compounds:
            sorted_storage = process_compounds(st.session_state.stored_compounds)
            st.session_state.processing_result = sorted_storage
            st.rerun()
        else:
            st.warning("Please add at least one compound before applying.")

# Results
if st.session_state.processing_result:
    st.header("Processed Compound Data")
    for info in st.session_state.processing_result:
        st.subheader(info['compound'])
        st.write(f"**Name:** {info['name']}")
        st.write(f"**SMILES:** {info['smiles']}")
        st.write(f"**Acid/Base Classification:** {info['acid_base_class']}")
        st.write(f"**State at Room Temperature:** {info['state_at_room_temp']}")
        st.write(f"**Pictograms:** {', '.join(info['pictograms']) if info['pictograms'] else 'None'}")
        st.write(f"**Safety Data:** {info['safety_data']}")

    st.markdown("---")
    st.subheader("Storage Groups")
    for group_name, compounds in st.session_state.storage_groups.items():
        if compounds:
            with st.expander(group_name):
                for c in compounds:
                    st.markdown(f"- {c}")
