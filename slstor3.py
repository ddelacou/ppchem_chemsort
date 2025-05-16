import streamlit as st
from typing import List, Dict, Any

# Placeholder for the chemsort_multiple_order function
def chemsort_multiple_order(compounds: List[str]) -> Dict[str, Any]:
    # Simulated processing logic
    return {compound: f"Processed data for {compound}" for compound in compounds}

# Initialize session state variables
if 'stored_compounds' not in st.session_state:
    st.session_state.stored_compounds = []
if 'compound_input' not in st.session_state:
    st.session_state.compound_input = ''

# Function to add a compound to the stored list
def add_compound():
    compound = st.session_state.compound_input.strip()
    if compound and compound not in st.session_state.stored_compounds:
        st.session_state.stored_compounds.append(compound)
        st.session_state.compound_input = ''  # Clear input field
        st.success(f"Added compound: {compound}")
    elif not compound:
        st.warning("Please enter a compound name.")
    else:
        st.info(f"Compound '{compound}' is already in the list.")

# Function to remove a compound from the stored list
def remove_compound(compound):
    if compound in st.session_state.stored_compounds:
        st.session_state.stored_compounds.remove(compound)
        st.success(f"Removed compound: {compound}")

st.title("Chemical Compound Input Interface")

# Sidebar for compound inputs
with st.sidebar:
    st.header("Input Compounds")

    # Text input with on_change callback to handle Enter key press
    st.text_input(
        label="Enter Compound Name",
        key='compound_input',
        on_change=add_compound
    )

    # Button to add compound (alternative to pressing Enter)
    st.button("Add Compound", on_click=add_compound)

    # Display stored compounds with remove buttons
    st.markdown("### Stored Compounds")
    if st.session_state.stored_compounds:
        for compound in st.session_state.stored_compounds:
            col1, col2 = st.columns([0.85, 0.15])
            with col1:
                st.write(compound)
            with col2:
                remove_button = st.button("❌", key=f"remove_{compound}")
                if remove_button:
                    remove_compound(compound)
                    st.rerun()
    else:
        st.write("No compounds added yet.")

    # Apply button to process compounds
    if st.button("Apply"):
        if st.session_state.stored_compounds:
            result = chemsort_multiple_order(st.session_state.stored_compounds)
            st.session_state['processing_result'] = result
            st.rerun()
        else:
            st.warning("Please add at least one compound before applying.")

# Main area to display processing results
if 'processing_result' in st.session_state:
    st.header("Processing Results")
    for compound, data in st.session_state['processing_result'].items():
        st.write(f"**{compound}**: {data}")
