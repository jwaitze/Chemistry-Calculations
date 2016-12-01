import sys, ast

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
    element, element_start, element_length = {'front': int(coefficient), 'element': '', 'back': 1}, 0, 0
    for c in range(len(element_object)):
        if element_object[c].isupper():
            element_start, element_length = c, 1
        elif element_object[c].islower():
            element_length += 1
    element['element'] = element_object[element_start:element_start + element_length]
    if element_start != 0:
        element['front'] *= int(element_object[:element_start])
    if len(element_object) > len(str(element['front'])) + len(element['element']):
        element['back'] = element_object[element_start + element_length:]
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
    components = [coefficient + ''.join(component) for component in components]
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
    result = [[], []]
    for c in range(len(left_components)):
        result[0].append(get_substance_components(left_components[c]))
    for c in range(len(right_components)):
        result[1].append(get_substance_components(right_components[c]))
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
    components, molar_mass = get_substance_components(substance), 0
    for component in components:
        row = locate_periodic_table_row(component)
        molar_mass += periodic_table[row]['atomic_weight']
    return molar_mass

def moles(mass, substance):
    return mass / molar_mass(substance)

def molarity(mass, substance, liters):
    return moles(mass, substance) / liters
    
if __name__ == '__main__':
##    left_components, right_components = get_reaction_components('2Al + 3CuSO4 -> Al2(SO4)3 + 3Cu')

    left_components, right_components = get_reaction_components('2Al2(Ca22SO4C3)33(SO4)2')
    if len(left_components):
        for component in left_components:
            print(component)
    if len(right_components):
        print()
        for component in right_components:
            print(component)
    
##    print()
##    for element in [[element['name'], element['electron_configuration'], element['boiling_point']] for element in periodic_table]:
##        print(''.join(['{}, '.format(detail) for detail in element])[:-2])
##    print()
##    for row in reactivity_series:
##        print(row)
##    print()
##    for row in electromotive_potentials:
##        print(row)
	
