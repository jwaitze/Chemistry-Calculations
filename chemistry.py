import sys, ast

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

def get_substance_components(substance):
    if ' ' in substance or '+' in substance or '-' in substance or '>' in substance or len(substance) < 1:
        return []
    components, coefficient = [], '1'
    for i in range(len(substance)):
        c = substance[i]
        upper, lower, digit, parenthesis = c.isupper(), c.islower(), c.isdigit(), c is '(' or c is ')'
        if not upper and not lower and not digit and not parenthesis:
            return False
        if i == 0 and digit:
            coefficient = c
        elif upper:
            components.append([c])
        elif lower or digit:
            components[-1].append(c)
    components = [''.join(component) for component in components]
    if coefficient != '1':
        components = [coefficient + component for component in components]
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

def get_file_contents(filepath):
    try:
        with open(filepath, 'r', encoding="utf8") as fp:
            return fp.read()
    except:
        return ''

def load_json_database(filepath):
    contents, json_data = get_file_contents(filepath), []
    if len(contents) < 1:
        print('test')
        return []
    contents = contents.split('\n')
    for content in contents:
        try:
            json_row = ast.literal_eval(content)
            json_data.append(json_row)
        except:
            pass
    return json_data

def is_valid_reaction(formula):
    if '->' not in formula:
        return False
    
if __name__ == '__main__':
    periodic_table = load_json_database('periodic_table.json')
    reactivity_series = load_json_database('reactivity_series.json')
    electromotive_potentials = load_json_database('electromotive_potentials.json')
    print(get_reaction_components('2Al + 3CuSO4 -> Al2(SO4)3 + 3Cu'))
    print()
    for element in [[element['name'], element['electron_configuration'], element['boiling_point']] for element in periodic_table]:
        print(''.join(['{}, '.format(detail) for detail in element])[:-2])
    print()
    for row in reactivity_series:
        print(row)
    print()
    for row in electromotive_potentials:
        print(row)
	
