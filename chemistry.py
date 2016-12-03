import sys, ast
from sympy import *

def get_file_contents(filepath):
    try:
        with open(filepath, 'r', encoding="utf8") as fp:
            return fp.read()
    except:
        return ''

def load_json_database(filepath):
    contents, json_data = get_file_contents(filepath), []
    if len(contents) < 1:
        print('Failed to load ' + filepath)
        sys.exit()
    contents = contents.split('\n')
    for content in contents:
        try:
            json_row = ast.literal_eval(content)
            json_data.append(json_row)
        except:
            pass
    return json_data

periodic_table = load_json_database('periodic_table.json')
reactivity_series = load_json_database('reactivity_series.json')
electromotive_potentials = load_json_database('electromotive_potentials.json')

absolute_zero = -273.15

def celcius_to_farenheit(temperature):
    return ((temperature*9)/5)+32

def farenheit_to_celcius(temperature):
    return ((temperature-32)*5)/9

def celcius_to_kelvin(temperature):
    return temperature-absolute_zero

def kelvin_to_celcius(temperature):
    return temperature+absolute_zero

def farenheit_to_kelvin(temperature):
    return farenheit_to_celcius(celcius_to_kelvin(temperature))

def kelvin_to_farenheit(temperature):
    return celcius_to_farenheit(kelvin_to_celcius(temperature))

def moles_to_atoms(moles):
    moles * (6.0221 * (10**23))

def atoms_to_moles(atoms):
    atoms / (6.0221 * (10**23))

def break_element_object(element_object, coefficient='1'):
    element, element_start, element_length = {'moles': int(coefficient), 'element': '', 'subscript': 1}, 0, 0
    for c in range(len(element_object)):
        if element_object[c].isupper():
            element_start, element_length = c, 1
        elif element_object[c].islower():
            element_length += 1
    element['element'] = element_object[element_start:element_start + element_length]
    if element_start != 0:
        element['moles'] *= int(element_object[:element_start])
    if len(element_object) >= len(str(element['moles'])) + len(element['element']):
        element['subscript'] = int(element_object[element_start + element_length:])
    return element

def break_component_section(component_section):
    close_parenthesis = component_section.index(')')
    section, components, coefficient = component_section[1:close_parenthesis], [], '1'
    if component_section[close_parenthesis + 1].isdigit():
        i, coefficient = 2, component_section[close_parenthesis + 1]
        try:
            while component_section[close_parenthesis + i].isdigit():
                coefficient += component_section[close_parenthesis + i]
                i += 1
        except:
            pass
    for i in range(len(section)):
        c = section[i]
        upper, lower, digit, = c.isupper(), c.islower(), c.isdigit()
        if not upper and not lower and not digit:
            return False
        if upper:
            components.append([c])
        elif lower or digit:
            components[-1].append(c)
    components_, components = [''.join(component) for component in components], []
    for component in components_:
        appended = False
        for c in range(len(component)):
            if component[c].isdigit():
                components.append(component[:c] + str(int(coefficient) * int(component[c:])))
                appended = True
                break
        if not appended:
            components.append(component + coefficient)
    return components

def get_substance_components(substance):
    if ' ' in substance or '+' in substance or '-' in substance or '>' in substance or len(substance) < 1:
        return []
    components, coefficient, parenthesis_opened, check_parenthesis_coefficient = [], '1', False, False
    for i in range(len(substance)):
        c = substance[i]
        upper, lower, digit, open_parenthesis, close_parenthesis = c.isupper(), c.islower(), c.isdigit(), c is '(', c is ')'
        if not upper and not lower and not digit and not open_parenthesis and not close_parenthesis:
            return False
        if parenthesis_opened:
            if close_parenthesis:
                parenthesis_opened = False
                check_parenthesis_coefficient = True
            continue
        if check_parenthesis_coefficient and digit:
            continue
        check_parenthesis_coefficient = False
        if i == 0 and digit:
            coefficient = c
        elif upper:
            components.append([c])
        elif lower or digit:
            components[-1].append(c)
        elif open_parenthesis:
            parenthesis_opened = True
            components.extend(break_component_section(substance[i:]))
    components = [break_element_object(''.join(component), coefficient) for component in components]
    return components

def get_reaction_components(formula):
    f, sides, left, right, right_components, left_components = formula.replace(' ', ''), [], [], [], [], [] 
    if '-->' in f:
        sides = f.split('-->')
    elif '->' in f:
        sides = f.split('->')
    elif '>' in f:
        sides = f.split('>')
    else:
        sides = [f]
    if len(sides) > 0:
        left = sides[0]
    if len(sides) > 1:
        right = sides[1]
    if len(left) > 0:
        left_components = left.split('+')
    if len(right) > 0:
        right_components = right.split('+')
    result = {'reactants': [], 'products': []}
    for c in range(len(left_components)):
        result['reactants'].append(get_substance_components(left_components[c]))
    for c in range(len(right_components)):
        result['products'].append(get_substance_components(right_components[c]))
    return result

def max_reactivity(metal_one, metal_two):
    if metal_one == metal_two:
        return metal_one
    metal_one_row, metal_two_row = -1, -1
    for row in range(len(reactivity_series)):
        if metal_one == reactivity_series[row]['symbol'] or metal_one == reactivity_series[row]['name']:
            metal_one_row = row
        if metal_two == reactivity_series[row]['symbol'] or metal_two == reactivity_series[row]['name']:
            metal_two_row = row
    if metal_one_row < 0 or metal_two_row < 0:
        return ''
    if metal_one_row == metal_two_row:
        return metal_one
    if metal_one_row < metal_two_row:
        return metal_one
    else:
        return metal_two

def locate_periodic_table_row(search_term, key=None):
    if key != None:
        for row in range(len(periodic_table)):
            if str(search_term) == str(periodic_table[row][key]):
                return row
        for row in range(len(periodic_table)):
            if str(search_term) in str(periodic_table[row][key]):
                return row
        return -1
    for row in range(len(periodic_table)):
        if type(search_term) is int:
            if search_term == periodic_table[row]['number'] or search_term == periodic_table[row]['atomic_number']:
                return row
        elif type(search_term) is str:
            if search_term == periodic_table[row]['name'] or search_term == periodic_table[row]['symbol']:
                return row
    for row in range(len(periodic_table)):
        for k in periodic_table[row]:
            if str(search_term) in str(periodic_table[row][k]):
                return row
    return -1

def get_periodic_table_details_list():
    return [key for key in periodic_table[0]]

def molar_mass(substance):
    components, grams_per_mole = get_substance_components(substance), 0
    for component in components:
        row = locate_periodic_table_row(component['element'])
        grams_per_mole += periodic_table[row]['atomic_weight'] * component['subscript']
    return grams_per_mole

def moles_substance(grams, substance):
    return grams / molar_mass(substance)

def mass_substance(moles, substance):
    return moles * molar_mass(substance)

def molarity(moles, liters):
    return moles / liters

def molarity_from_mass(grams, substance, liters):
    return moles_substance(grams, substance) / liters
                       
def print_moles_substance(grams, substance):
    print(str(grams) + ' g of ' + substance + ' = ' + str(moles_substance(grams, substance)) + ' moles')

def mmhg_to_atm(mmhg):
    return mmhg / 760

def torr_to_atm(torr):
    return mmhg_to_atm(torr)

def psi_to_atm(psi):
    return psi / 14.6959

def kpa_to_atm(kpa):
    return kpa / 101.325

def inh2o_to_atm(inh2o):
    return inh2o / 407.189

def mm_to_in(mm):
    return mm * 0.0393701

def moles_ideal_gas(liters):
    return liters / 22.4

def ideal_gas_initial_final_state(atm1, liters1, moles1, kelvin1, atm2, liters2, moles2, kelvin2):
    R = 0.08205746
    P1, V1, n1, T1, P2, V2 ,n2, T2 = [-1 if a == None else a for a in [atm1, liters1, moles1, kelvin1, atm2, liters2, moles2, kelvin2]]
    numerator_left, denominator_left = P1 * V1, P2 * V2
    numerator_right, denominator_right = n1 * R * T1, n2 * R * T2
    if numerator_left < 0:
        return ((numerator_right / denominator_right) * denominator_left) / (numerator_left * -1)
    elif numerator_right < 0:
        return ((numerator_left / denominator_left) * denominator_right) / (numerator_right * -1)
    elif denominator_left < 0:
        return (numerator_left / (numerator_right / denominator_right)) / (denominator_left * -1)
    elif denominator_right < 0:
        return (numerator_right / (numerator_left / denominator_left)) / (denominator_right * -1)

def ideal_gas(atm, liters, moles, kelvin):
    R = 0.08205746
    P, V, n, T = [-1 if a == None else a for a in [atm, liters, moles, kelvin]]
    left, right = P * V, n * R * T
    if left < 0:
        return right / (left * -1)
    elif right < 0:
        return left / (right * -1)

def gas_density(grams_per_liter, kelvin, grams_per_mole, atm):
    R = 0.08205746
    d, T, M, P = [-1 if a == None else a for a in [grams_per_liter, kelvin, grams_per_mole, atm]]
    left, right = R * d * T, M * P
    if left < 0:
        return right / (left * -1)
    elif right < 0:
        return left / (right * -1)

def balance_chemical_equation(reaction_components):
    equation_components, equations, i = {'reactants': [], 'products': []}, [], 0
    for side in reaction_components:
        for component in reaction_components[side]:
            equation_components[side].append({'variable': chr(ord('a')+i), 'component': component})
            i += 1
    for component in equation_components['reactants']:
        for element in component['component']:
            equations.append({'reactants': str(element['subscript']) + '*' + component['variable'] + '+', 'products': ''})
            for side in equation_components:
                for c in equation_components[side]:
                    if c['variable'] == component['variable']: # i will likely have to change this
                        continue
                    for e in c['component']:
                        if e['element'] == element['element']:
                            equations[-1][side] += str(e['subscript']) + '*' + c['variable'] + '+'
            for side in equations[-1]:
                equations[-1][side] = sympify(equations[-1][side][:-1])
    for e in equations:
        if '+' not in str(e['reactants']) and '+' not in str(e['products']):
            e['reactants'] = e['reactants'].subs(symbols(str(e['reactants'])[-1]), 1)
            break
    for e in equations:
        print(e)
    return reaction_components
    
if __name__ == '__main__':
    reaction_components = get_reaction_components('CH4 + O2 --> CO2 + H2O')
    balanced_components = balance_chemical_equation(reaction_components)
##    for component in balanced_components:
##        for element in component:
##            print(element)
##        print()
    sys.exit()

##    2 Moles of a gas at a pressure of 2.00 atm occupies a volume of 22.4 L.
##    The temperature is 293 K. What will the pressure be if the volume is held
##    constant and the temperature is raised to 348K?
    print(ideal_gas_initial_final_state(2, 22.4, 2, 293, None, 22.4, 2, 348))

##    How many moles of Cl2 are in a container (usually called a vessel) if the
##    pressure is 2.5 atm, the temperature is 27 C, and the volume is 50 L? R is,
##    of course, 0.0821 L atm mole-1 K- 1.
    print(ideal_gas(2.5, 50, None, celcius_to_kelvin(27)))

##    The density of a gas is measured at 1.853 g / L at 745.5 mmHg and 23.8 °C. What is its molar mass?
    print(gas_density(1.853, celcius_to_kelvin(23.8), None, mmhg_to_atm(745.5)))

##    Anhydrous aluminum chloride sublimes at high temperatures.
##    What density will the vapor have at 225 degrees Celsius and 0.939 atm of pressure?
    print(gas_density(None, celcius_to_kelvin(225), molar_mass('AlCl3'), 0.939))
    sys.exit()
    
    #print_moles_substance(5, 'Al2(SO4)3')
    
##    left_components, right_components = get_reaction_components('2Al + 3CuSO4 -> Al2(SO4)3 + 3Cu')

##    left_components, right_components = get_reaction_components('2Al2(Ca22SO4C3)33(SO4)2')
##    if len(left_components):
##        for component in left_components:
##            print(component)
##    if len(right_components):
##        print()
##        for component in right_components:
##            print(component)
    
##    print()
##    for element in [[element['name'], element['electron_configuration'], element['boiling_point']] for element in periodic_table]:
##        print(''.join(['{}, '.format(detail) for detail in element])[:-2])
##    print()
##    for row in reactivity_series:
##        print(row)
##    print()
##    for row in electromotive_potentials:
##        print(row)
	
