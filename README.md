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
print_lookup(term)

print_moles_substance(grams, substance)  
print_mass_substance(moles, substance)  
print_balance_chemical_formula(formula)  
print_stoichiometry(formula, reference_component, moles_reference_component, target_component)  
print_ideal_gas(atm1, liters1, moles1, kelvin1, atm2=None, liters2=None, moles2=None, kelvin2=None)  
print_mass_percentages(substance)  
print_periodic_table_retrieve(search_term, element)  
print_periodic_table_closest(search_term, max_results, column_search_term=None)

## Usage Examples

Returns data available from the infoboxes on Wikipedia:  
print_lookup('Al2(SO4)3')  
print_lookup('CuSO4')  
print_lookup('H2SO4')  
print_lookup('Hydrogen')

Returns the amount of moles for a given mass of a substance
print_moles_substance(36.4, 'CuSO4')  
print_moles_substance(100, 'H2O')  
print_moles_substance(molar_mass('H2O'), 'H2O') This will always return 1 as long as the input substances match  

Returns the mass of a given amount of moles of a substance  
print_mass_substance(2, 'MgSO4(H2O)7') mass of 2 moles of magnesium sulfate heptahydrate  

Any of the following formats are fine:  
print_balance_chemical_formula('C2H6 + O2 --> CO2 + H2O')  
print_balance_chemical_formula('C4H10 + O2 -> CO2 + H2O')  
print_balance_chemical_formula('FeCl3 + NH4OH > Fe(OH)3 + NH4Cl')  
print_balance_chemical_formula('P2I4 + P4 + H2O --> PH4I + H3PO4')  
print_balance_chemical_formula('Cu+HNO3-->Cu(NO3)2+NO2+H2O')  
print_balance_chemical_formula('S+HNO3-> H2SO4 + NO2 + H2O')

Solves stoichiometry calculations (automates balancing, moles conversions, and substance ratios):
print_stoichiometry('C6H12O6 + O2 --> CO2 + H2O', 'O2', 1.03*(10**-2), 'C6H12O6') Grams of glucose the human body consumes each minute if 1.03 x 10^-2 mol O2 is used for cellular respiration each minute  
print_stoichiometry('CuSO4 + Al --> Al2(SO4)3 + Cu', 'CuSO4', moles_substance(5, 'CuSO4(H2O)5'), 'Al') How much aluminum is needed to react with 5 grams of copper(II) sulfate pentahydrate?  
print(moles_to_liters_gas(stoichiometry('C2H6 + O2 --> CO2 + H2O', 'C2H6', moles_substance(20, 'C2H6'), 'O2')['moles']), 'liters of O2 gas at STP') How many liters of O2 is needed to react with 20 g of ethane?

