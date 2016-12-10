# Chemistry Utility
This is what the chemistry database application building has been for - this application serves to automate much of the possibly tedious lookups and calculations involved in general chemistry.

## Features of Finished Application
- stoichiometry
- balancing of formulae
- molarity calculations
- information lookup
- temperature conversions
- mass and volume conversions
- free energy
- gas laws
- reactivity
- entropy
- possibly some templates for common procedures like titration or gas generation measurement

## Current API:
get_file_contents(filepath)  
load_json_database(filepath)

break_element_object(element_object, coefficient='1')  
break_component_section(component_section)  
get_substance_components(substance)  
get_reaction_components(formula)  
reaction_components_to_string(reaction_components)

locate_periodic_table_row(search_term)  
get_periodic_table_details_list()  
locate_periodic_table_column(search_term)  
periodic_table_retrieve(term, element)  
periodic_table_closest(term, max_results, column_search_term=None)

celsius_to_fahrenheit(temperature)  
fahrenheit_to_celsius(temperature)  
celsius_to_kelvin(temperature)  
kelvin_to_celsius(temperature)  
fahrenheit_to_kelvin(temperature)  
kelvin_to_fahrenheit(temperature)

moles_to_atoms(moles)  
atoms_to_moles(atoms)

molar_mass(substance)  
moles_substance(grams, substance)  
mass_substance(moles, substance)  
molarity(moles, liters)  
moles_from_molarity(molarity, liters)
molarity_from_mass(grams, substance, liters)  
moles_ideal_gas(liters)  
moles_to_liters_gas(moles)  
get_mass_percentages(substance)

mmhg_to_atm(mmhg)  
torr_to_atm(torr)  
psi_to_atm(psi)  
kpa_to_atm(kpa)  
inh2o_to_atm(inh2o)  
mm_to_in(mm)  

ideal_gas_initial_final_state(atm1, liters1, moles1, kelvin1, atm2, liters2, moles2, kelvin2)  
ideal_gas(atm, liters, moles, kelvin)  
gas_density(grams_per_liter, kelvin, grams_per_mole, atm)

max_reactivity(metal_one, metal_two)

balance_chemical_formula(reaction_components)  
stoichiometry(formula, reference_component, moles_limiting_reagent, target_component)

osmotic_pressure(molarity, kelvin, vanthoff_factor)

process_raw_wikipedia_table(detail, strip_extraneous_characters)
get_data_from_wikipedia(term)
lookup_term_wikipedia(term)
print_lookup_term_wikipedia(term)

print_moles_substance(grams, substance)  
print_mass_substance(moles, substance)  
print_balance_chemical_formula(formula)  
print_stoichiometry(formula, reference_component, moles_reference_component, target_component)  
print_ideal_gas(atm1, liters1, moles1, kelvin1, atm2=None, liters2=None, moles2=None, kelvin2=None)  
print_mass_percentages(substance)  
print_periodic_table_retrieve(search_term, element)  
print_periodic_table_closest(search_term, max_results, column_search_term=None)
