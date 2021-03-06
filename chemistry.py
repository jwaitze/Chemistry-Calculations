import sys, ast, math, requests
from sympy import *
from bs4 import BeautifulSoup

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

def celsius_to_fahrenheit(temperature):
    return ((temperature*9)/5)+32

def fahrenheit_to_celsius(temperature):
    return ((temperature-32)*5)/9

def celsius_to_kelvin(temperature):
    return temperature-absolute_zero

def kelvin_to_celsius(temperature):
    return temperature+absolute_zero

def fahrenheit_to_kelvin(temperature):
    return fahrenheit_to_celsius(celsius_to_kelvin(temperature))

def kelvin_to_fahrenheit(temperature):
    return celsius_to_fahrenheit(kelvin_to_celsius(temperature))

def moles_to_atoms(moles):
    return moles * (6.0221 * (10**23))

def atoms_to_moles(atoms):
    return atoms / (6.0221 * (10**23))

def break_element_object(element_object, coefficient='1'):
    element, element_start, element_length = {'moles': int(coefficient), 'element': '', 'subscript': 1}, 0, 0
    for i, c in enumerate(element_object):
        if c.isupper():
            element_start, element_length = i, 1
        elif c.islower():
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
    elements = []
    for component in components:
        if component['element'] not in elements:
            elements.append(component['element'])
    for e in range(len(elements)):
        elements[e] = {'element': elements[e], 'subscript': 0}
        for component in components:
            if component['element'] == elements[e]['element']:
                elements[e]['moles'] = component['moles']
                elements[e]['subscript'] += component['subscript']
    return elements

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

def locate_periodic_table_row(search_term):
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

def locate_periodic_table_column(search_term):
    keys, results = get_periodic_table_details_list(), []
    for i, k in enumerate(keys):
        if search_term in k:
            results.append({'column': i, 'header': k})
    return results

def periodic_table_retrieve(term, element):
    results = []
    columns = locate_periodic_table_column(term)
    row = locate_periodic_table_row(element)
    for column in columns:
        results.append({'value': periodic_table[row][column['header']], 'header': column['header']})
    return results

def periodic_table_closest(value_term, max_results, column_search_term=None):
    results, keys = [], get_periodic_table_details_list()
    for result in range(max_results):
        closest = {'value': 10**6}
        for key in keys:
            if column_search_term is not None:
                if column_search_term not in key:
                    continue
            for element in periodic_table:
                if type(element[key]) is not int and type(element[key]) is not float:
                    continue
                next_closest = {'value': element[key], 'element': element['name'], 'symbol': element['symbol'], 'detail': key}
                if math.fabs(value_term - element[key]) <= math.fabs(value_term - closest['value']) and next_closest not in results:
                    closest = next_closest.copy()
        results.append(closest)
    return results

def molar_mass(substance):
    components, grams_per_mole = get_substance_components(substance), 0
    for component in components:
        row = locate_periodic_table_row(component['element'])
        grams_per_mole += periodic_table[row]['atomic_weight'] * component['subscript']
    return grams_per_mole

def print_molar_mass(substance):
    print(molar_mass(substance))

def moles_substance(grams, substance):
    return grams / molar_mass(substance)

def mass_substance(moles, substance):
    return moles * molar_mass(substance)

def molarity(moles, liters):
    return moles / liters

def moles_from_molarity(molarity, liters):
    return molarity * liters

def molarity_from_mass(grams, substance, liters):
    return moles_substance(grams, substance) / liters
                       
def print_moles_substance(grams, substance):
    moles = moles_substance(grams, substance)
    print(str(grams) + ' grams of ' + substance + ' = ' + str(moles) + ' moles')

def print_mass_substance(moles, substance):
    grams = mass_substance(moles, substance)
    print(str(moles) + ' moles of ' + substance + ' = ' + str(grams) + ' grams')

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

def moles_to_liters_gas(moles):
    return 22.4 * moles

def get_mass_percentages(substance):
    components, result = get_substance_components(substance), []
    for component in components:
        partial_substance = component['element'] + str(component['subscript'])
        percentage = molar_mass(partial_substance) / molar_mass(substance)
        result.append({'substance': partial_substance, 'percentage': percentage})
    return result

def print_get_mass_percentages(substance):
    percentages = get_mass_percentages(substance)
    for p in percentages:
        print(p['substance'] + ':', p['percentage'])

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

def balance_chemical_formula(reaction_components):
    equation_components, equations, i, coefficients = {'reactants': [], 'products': []}, [], 0, {}
    for side in reaction_components:
        for component in reaction_components[side]:
            equation_components[side].append({'variable': chr(ord('a')+i), 'component': component})
            i += 1
    for component in equation_components['reactants']:
        for element in component['component']:
            equations.append({'reactants': str(element['subscript']) + '*' + component['variable'] + '+', 'products': ''})
            for side in equation_components:
                for c in equation_components[side]:
                    if c['variable'] == component['variable']:
                        continue
                    for e in c['component']:
                        if e['element'] == element['element']:
                            equations[-1][side] += str(e['subscript']) + '*' + c['variable'] + '+'
            for side in equations[-1]:
                equations[-1][side] = sympify(equations[-1][side][:-1])
    new_equations = [Eq(e['reactants'], e['products']) for e in equations]
    solved = solve(new_equations)
    if type(solved) is list:
        solved = solved[0]
    equations = [{'reactants': k, 'products': solved[k]} for k in solved]
    equations_str, best_variable = ''.join([str(e['reactants']) + str(e['products']) for e in equations]), '-'
    for e in equations:
        if '+' not in str(e['reactants']) and '+' not in str(e['products']):
            for side in e:
                v = str(e[side])[-1]
                if not v.isupper() and not v.islower():
                    continue
                if equations_str.count(v) >= equations_str.count(best_variable):
                    best_variable = v
    symbol = symbols(best_variable)
    for e in equations:
        for side in e:
            if str(symbol) in str(e[side]):
                e[side] = e[side].subs(symbol, 1)
                coefficients[str(symbol)] = 1
    while len(coefficients) != i:
        for e in equations:
            if Eq(e['reactants'], e['products']) == True:
                continue
            for j in range(i):
                v = chr(ord('a')+j)
                answer = solveset(Eq(e['reactants'], e['products']), symbols(v))
                for s in answer:
                    if ' ' in str(s) or '*' in str(s):
                        continue
                    if not str(s).isdigit() and '/' not in str(s) and '.' not in str(s):
                        continue
                    if not str(s)[0].isdigit() or not str(s)[-1].isdigit():
                        continue
                    if v not in coefficients:
                        if '/' in str(s):
                            coefficients[v] = str(s)#eval(str(s))
                        elif '.' in str(s):
                            coefficients[v] = float(s)
                        else:
                            coefficients[v] = int(s)
                    break
        for coefficient in coefficients:
            for e in equations:
                for side in e:
                    e[side] = e[side].subs(symbols(coefficient), coefficients[coefficient])
    result, i = {}, 0
    for side in reaction_components:
        result[side] = []
        for component in reaction_components[side]:
            result[side].append([])
            symbol = chr(ord('a')+i)
            for element in component:
                result[side][-1].append(element.copy())
                result[side][-1][-1]['moles'] = coefficients[symbol]
            i += 1
    return result

def reaction_components_to_string(reaction_components):
    formula = ''
    for component in reaction_components['reactants']:
        if component[0]['moles'] != 1:
            formula += str(component[0]['moles'])
        for element in component:
            formula += element['element']
            if element['subscript'] != 1:
                formula += str(element['subscript'])
        formula += ' + '
    formula = formula[:-3] + ' --> '
    for component in reaction_components['products']:
        if component[0]['moles'] != 1:
            formula += str(component[0]['moles'])
        for element in component:
            formula += element['element']
            if element['subscript'] != 1:
                formula += str(element['subscript'])
        formula += ' + '
    formula = formula[:-3]
    return formula

def print_balance_chemical_formula(formula):
    reaction_components = get_reaction_components(formula)
    balanced_components = balance_chemical_formula(reaction_components)
    balanced_formula = reaction_components_to_string(balanced_components)
    print(balanced_formula)

def stoichiometry(formula, reference_component, moles_limiting_reagent, target_component):
    reaction_components = get_reaction_components(formula)
    balanced_components = balance_chemical_formula(reaction_components)
    balanced_formula = reaction_components_to_string(balanced_components)
    ratio = {}
    for side in balanced_components:
        for component in balanced_components[side]:
            substance = ''
            for element in component:
                substance += element['element']
                if element['subscript'] > 1:
                    substance += str(element['subscript'])
            if substance == reference_component:
                ratio['reference_component'] = component[0]['moles']
            elif substance == target_component:
                ratio['target_component'] = component[0]['moles']
    for c in ratio:
        if type(ratio[c]) is str:
            ratio[c] = eval(ratio[c])
    moles_target_component = (moles_limiting_reagent / ratio['reference_component']) * ratio['target_component']
    mass_target_component = moles_target_component * molar_mass(target_component)
    return {'balanced_formula': balanced_formula, 'moles': moles_target_component, 'grams': mass_target_component}

def print_stoichiometry(formula, reference_component, moles_reference_component, target_component):
    result = stoichiometry(formula, reference_component, moles_reference_component, target_component)
    print(result['balanced_formula'])
    print(reference_component, '->', moles_reference_component, 'moles or',
          mass_substance(moles_reference_component, reference_component), 'grams')
    print(target_component, '->', result['moles'], 'moles or', result['grams'], 'grams')

def print_ideal_gas(atm1, liters1, moles1, kelvin1, atm2=None, liters2=None, moles2=None, kelvin2=None):
    inputs, result = [atm1, liters1, moles1, kelvin1, atm2, liters2, moles2, kelvin2], 0
    if set([None]) == set(inputs[4:]):
        result = ideal_gas(inputs[0], inputs[1], inputs[2], inputs[3])
    else:
        result = ideal_gas_initial_final_state(inputs[0], inputs[1], inputs[2], inputs[3],
                                            inputs[4], inputs[5], inputs[6], inputs7)
    print(result)

def print_mass_percentages(substance):
    mass_percentages = get_mass_percentages(substance)
    for percentage in mass_percentages:
        print(percentage['substance'], '=', percentage['percentage'])

def print_periodic_table_retrieve(search_term, element):
    search_results = periodic_table_retrieve(search_term, element)
    for row in search_results:
        print(row)

def print_periodic_table_closest(search_term, max_results, column_search_term=None):
    results = periodic_table_closest(search_term, max_results, column_search_term)
    for i, row in enumerate(results):
        print(str(i) + ':', row['detail'], 'of', row['element'] + '(' + row['symbol'] + ')', 'is', row['value'])

def gibbs_free_energy():
    return

def specific_heat():
    return

def entropy():
    return

def enthalpy():
    return

def osmotic_pressure(molarity, kelvin, vanthoff_factor):
    R = 0.08206 # atm / mol K
    return vanthoff_factor * molarity * R * kelvin

def process_raw_wikipedia_table(detail, strip_extraneous_characters):
    body = detail[1].replace('\xa0', ' ').replace('\u200b', '').strip()
    title = detail[0].lower().replace(', ', '_').replace('\xa0', ' ')
    title = title.replace(' ', '_').replace('.', '').strip()
    title = title.replace('\'', '').replace('\n', '')
    for i in range(1, 25):
        title = title.replace('[' + str(i) + ']', '')
        body = body.replace('[' + str(i) + ']', '')
    if body.replace('-', '').replace('.', '').isdigit() and '-' not in body[1:] and body.count('.') < 2:
        if '.' in body:
            body = float(body)
        else:
            body = int(body)
    extraneous_delimeters = ']) '
    if strip_extraneous_characters and type(body) is str and sum([1 for e in extraneous_delimeters if e in body]) != 0:
        stripped = body
        for e in extraneous_delimeters:
            if e in stripped:
                stripped = stripped[:stripped.index(e)]
        strip_chars = '[(,' + extraneous_delimeters
        for c in strip_chars:
            stripped = stripped.replace(c, '')
        if stripped.replace('.', '').replace('-', '').isdigit() and '-' not in stripped[1:] and stripped.count('.') < 2:
            if '.' in stripped:
                body = float(stripped)
            else:
                body = int(stripped)
    if type(body) is str:
        split_delimeters = [', ', '\n']
        for s in split_delimeters:
            if s in body:
                body = body.split(s)
                body = [d for d in body if len(d) != 0 and d != ' ']
    if type(body) is list:
        body = ', '.join(body)
    return title, body

def get_data_from_wikipedia(term):
    data = []
    url = 'https://en.wikipedia.org/wiki/%s' % term.replace(' ', '_')
    r = requests.get(url)
    soup = BeautifulSoup(r.content, 'html.parser')
    infoboxes = soup.find_all('table', 'infobox')
    trs = infoboxes[0].find_all('tr')
    for tr in trs:
        ths = tr.find_all('th')
        tds = tr.find_all('td')
        try:
            if len(ths[0].text) > 0 and len(tds[0].text) > 0:
                data.append([ths[0].text, tds[0].text])
        except:
            try:
                if len(tds[0].text) > 0 and len(tds[1].text) > 0:
                    data.append([tds[0].text, tds[1].text])
            except:
                pass
    for detail in data:
        detail[0], detail[1] = process_raw_wikipedia_table(detail, False)
    return data

def lookup_term_wikipedia(term):
    data = get_data_from_wikipedia(term)
    return [row for row in data if len(row[0]) > 0 and len(str(row[1])) > 0]

def print_lookup_term_wikipedia(term):
    data = [[row[0] + ':', row[1]] if row[0][-1] != ':' else row for row in lookup_term_wikipedia(term)]
    for row in data:
        print(row[0], row[1])

def print_lookup(term):
    print_lookup_term_wikipedia(term)

def c(*x):
    xl = len(x)
    if xl == 0:
        print(get_file_contents('print_help.txt'))
    elif xl == 1:
        if type(x[0]) is str:
            if '+' in x[0] and '>' in x[0]:
                print_balance_chemical_formula(x[0])
            elif x[0][:4].lower() == 'mm: ':
                print_molar_mass(x[0][4:])
            elif x[0][:4] == 'm%: ':
                print_get_mass_percentages(x[0][4:])
            else:
                print_lookup(term)
        else:
            print_periodic_table_closest(x[0], 10)
    elif xl == 2:
        print_moles_substance(x[0], x[1])
    elif xl == 4:
        print_stoichiometry(x[0], x[1], x[2], x[3])
    
if __name__ == '__main__':

    sys.exit()

    print_lookup_term_wikipedia('Al2(SO4)3')

    sys.exit()
    
    print_balance_chemical_formula('C2H6 + O2 --> CO2 + H2O')
    print_balance_chemical_formula('C4H10 + O2 --> CO2 + H2O')
    print_balance_chemical_formula('S8 + O3 --> SO2')
    print_balance_chemical_formula('KFe3AlSi3O10(OH)2 + Cu + O2 + H2S --> KAlSi3O8 + CuFeS2 + H2O')
    print_balance_chemical_formula('Ca(OH)2 + H3PO4 --> Ca3(PO4)2 + H2O')
    print_balance_chemical_formula('KClO3 --> KCl + O2')
    print_balance_chemical_formula('Fe + O2 --> Fe2O3')
    print_balance_chemical_formula('C2H6 + O2 --> CO2 + H2O')
    print_balance_chemical_formula('Zn + HCl --> ZnCl2 + H2')
    print_balance_chemical_formula('H2 + Cl2 --> HCl')
    print_balance_chemical_formula('Al2(CO3)3 + H3PO4 --> AlPO4 + CO2 + H2O')
    print_balance_chemical_formula('C5H12 + O2 --> CO2 + H2O')
    print_balance_chemical_formula('CuSO4 + Al --> Al2(SO4)3 + Cu')
    print_balance_chemical_formula('H2 + O2 --> H2O')
    print_balance_chemical_formula('H2 + O3 --> H2O')
    print_balance_chemical_formula('S8 + F2 --> SF6')
    print_balance_chemical_formula('P4 + O2 --> P2O5')
    print_balance_chemical_formula('FeCl3 + NH4OH --> Fe(OH)3 + NH4Cl')
    print_balance_chemical_formula('P2I4 + P4 + H2O --> PH4I + H3PO4')
    print_balance_chemical_formula('Cu + HNO3 --> Cu(NO3)2 + NO2 + H2O')
    print_balance_chemical_formula('S + HNO3 --> H2SO4 + NO2 + H2O')
    
    sys.exit()

    print(osmotic_pressure(molarity(moles_substance(5, 'C6H12O6'), 0.1), fahrenheit_to_kelvin(70.1), 1))
    print(osmotic_pressure(molarity(moles_substance(210, 'C12H22O11'), 0.1), fahrenheit_to_kelvin(70.1), 1))
    print(osmotic_pressure(molarity(moles_substance(74.5, 'CaCl2'), 0.1), fahrenheit_to_kelvin(70.1), 3))
    print(osmotic_pressure(molarity(moles_substance(35.9, 'NaCl'), 0.1), fahrenheit_to_kelvin(70.1), 2))

    sys.exit()

## 1. The human body needs at least 1.03 x 10-2 mol O2 every minute. If all of this oxygen
## is used for the cellular respiration reaction that breaks down glucose, how many grams of
## glucose does the human body consume each minute?
    print_stoichiometry('C6H12O6 + O2 --> CO2 + H2O', 'O2', 1.03*(10**-2), 'C6H12O6')

    sys.exit()

    result = print_stoichiometry('KClO3 --> KCl + O2', 'KClO3', moles_substance(39.498-38.388, 'KClO3'), 'O2')
    print('or', ideal_gas(1, None, result['moles'], fahrenheit_to_kelvin(70.7)), 'liters of O2 gas at 70.7 °F')

    sys.exit()

    print_stoichiometry('KI + Pb(NO3)2 --> PbI2 + KNO3', 'KI', moles_substance(1.187, 'KI'), 'PbI2')

    sys.exit()

    print_stoichiometry('CuSO4 + Al --> Al2(SO4)3 + Cu', 'CuSO4', moles_substance(5, 'CuSO4(H2O)5'), 'Al')

    sys.exit()

    print(moles_to_liters_gas(stoichiometry('C2H6 + O2 --> CO2 + H2O', 'C2H6', moles_substance(20, 'C2H6'), 'O2')['moles']), 'liters of O2 gas at STP')

    sys.exit()
    
    formula = 'S + HNO3 --> H2SO4 + NO2 + H2O'
    print(formula)
    reaction_components = get_reaction_components(formula)
    balanced_components = balance_chemical_formula(reaction_components)
##    for side in reaction_components:
##        for component in reaction_components[side]:
##            print(component)
##        print()
##    print('\n')
##    for side in balanced_components:
##        for component in balanced_components[side]:
##            print(component)
##        print()
##    print()
    balanced_formula = reaction_components_to_string(balanced_components)
    print(balanced_formula)
    sys.exit()

##    2 Moles of a gas at a pressure of 2.00 atm occupies a volume of 22.4 L.
##    The temperature is 293 K. What will the pressure be if the volume is held
##    constant and the temperature is raised to 348K?
    print(ideal_gas_initial_final_state(2, 22.4, 2, 293, None, 22.4, 2, 348))

##    How many moles of Cl2 are in a container (usually called a vessel) if the
##    pressure is 2.5 atm, the temperature is 27 C, and the volume is 50 L? R is,
##    of course, 0.0821 L atm mole-1 K- 1.
    print(ideal_gas(2.5, 50, None, celsius_to_kelvin(27)))

##    The density of a gas is measured at 1.853 g / L at 745.5 mmHg and 23.8 °C. What is its molar mass?
    print(gas_density(1.853, celsius_to_kelvin(23.8), None, mmhg_to_atm(745.5)))

##    Anhydrous aluminum chloride sublimes at high temperatures.
##    What density will the vapor have at 225 degrees Celsius and 0.939 atm of pressure?
    print(gas_density(None, celsius_to_kelvin(225), molar_mass('AlCl3'), 0.939))
    sys.exit()
    
    print_moles_substance(5, 'Al2(SO4)3')
    
##    print()
##    for element in [[element['name'], element['electron_configuration'], element['boiling_point']] for element in periodic_table]:
##        print(''.join(['{}, '.format(detail) for detail in element])[:-2])
##    print()
##    for row in reactivity_series:
##        print(row)
##    print()
##    for row in electromotive_potentials:
##        print(row)
	
