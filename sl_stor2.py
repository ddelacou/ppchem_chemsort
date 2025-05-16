import streamlit as st
from typing import Dict, Any, List

chemical_date = {} #result of the chemsort_multiple_order_3 function

# Initialize session state
if 'selected_compounds' not in st.session_state:
    st.session_state.selected_compounds = {}

# Add compound to selection
def select_compound(compound):
    if compound['name'] not in st.session_state.selected_compounds:
        st.session_state.selected_compounds[compound['name']] = compound

# Remove compound from selection
def remove_compound(name):
    if name in st.session_state.selected_compounds:
        del st.session_state.selected_compounds[name]

# UI layout
st.set_page_config(layout="wide")
st.title("Chemical Category Explorer")

# Layout: 3 categories per row
categories = list(chemical_data.keys())
cols_per_row = 3

for i in range(0, len(categories), cols_per_row):
    row = categories[i:i + cols_per_row]
    row_cols = st.columns(cols_per_row)

    for col, category in zip(row_cols, row):
        with col:
            with st.expander(f"{category.replace('_', ' ').title()}", expanded=False):
                st.markdown("<h5 style='margin-bottom: 0.5rem;'>Gas</h5>", unsafe_allow_html=True)
                for compound in chemical_data[category].get('gas', []):
                    if st.button(compound['name'], key=f"{category}_gas_{compound['name']}"):
                        select_compound(compound)

                st.markdown("<h5 style='margin-bottom: 0.5rem;'>Liquid</h5>", unsafe_allow_html=True)
                for compound in chemical_data[category].get('liquid', []):
                    if st.button(compound['name'], key=f"{category}_liquid_{compound['name']}"):
                        select_compound(compound)

                st.markdown("<h5 style='margin-bottom: 0.5rem;'>Solid</h5>", unsafe_allow_html=True)
                for compound in chemical_data[category].get('solid', []):
                    if st.button(compound['name'], key=f"{category}_solid_{compound['name']}"):
                        select_compound(compound)

# Display selected compound info
if st.session_state.selected_compounds:
    st.markdown("---")
    st.subheader("Selected Compound Details")
    detail_cols = st.columns(len(st.session_state.selected_compounds))

    for i, (name, compound) in enumerate(st.session_state.selected_compounds.items()):
        with detail_cols[i]:
            st.markdown(f"### {compound['name']}")
            st.markdown(f"**State:** {compound['state_room_temp'].capitalize()}")
            st.markdown(f"**Acid/Base Class:** {compound['acid_base_class'].capitalize()}")
            st.markdown(f"**Pictograms:** {', '.join(compound['pictograms']) if compound['pictograms'] else 'None'}")
            if st.button("❌ Remove", key=f"remove_{compound['name']}"):
                remove_compound(compound['name'])
                st.rerun()  # Updated here
