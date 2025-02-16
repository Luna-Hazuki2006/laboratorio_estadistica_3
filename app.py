import csv
import pandas
from scipy.stats import f, studentized_range
from pprint import pprint
import math

info = {'∑xt': 0, '∑xt^2': 0, '(∑xt)^2': 0, 'nt': 0, '(∑xt)^2/nt': 0, 'variables': []}

with open('data.csv', newline='', encoding='utf-8') as data: 
    datos = csv.DictReader(data)
    nombres = list(datos.fieldnames)
    for variable in nombres: 
        info['variables'].append({'nombre': variable, 'lista': []})
    [[[info['variables'][i]['lista'].append(float(x[nombres[i]]))] for i in range(len(nombres))] for x in datos]

for variable in info['variables']: 
    variable['lista^2'] = [x ** 2 for x in variable['lista']]
    variable['∑x'] = sum(variable['lista'])
    variable['∑x^2'] = sum(variable['lista^2'])
    variable['(∑x)^2'] = variable['∑x'] ** 2
    variable['n'] = len(variable['lista'])
    variable['((∑x)^2)/n'] = variable['(∑x)^2'] / variable['n']
    variable['x'] = variable['∑x'] / variable['n']
    info['∑xt'] += variable['∑x']
    info['∑xt^2'] += variable['∑x^2']
    info['(∑xt)^2'] += variable['(∑x)^2']
    info['nt'] += variable['n']
    info['(∑xt)^2/nt'] += variable['((∑x)^2)/n']

C = (info['∑xt'] ** 2) / info['nt']
SCT = info['∑xt^2'] - C
SCTR = info['(∑xt)^2/nt'] - C
SCE = info['∑xt^2'] - info['(∑xt)^2/nt']

glt = len(info['variables']) - 1
gle = info['nt'] - len(info['variables'])

if (glt + gle) != (info['nt'] - 1): print('ALGO ESTÁ MUY MALLLLLLLLLLLLLLLLLLL')

ftab = f.ppf(q=0.95, dfn=glt, dfd=gle)

MCTR = SCTR / glt
MCE = SCE / gle

nj = len(info['variables'][0]['lista'])

fcal = MCTR / MCE

DHS = (studentized_range.ppf(q=0.95, k=glt + 1, df=gle)) * (math.sqrt(MCE / nj))

pprint(info)
print(f'C {C}')
print(f'SCT {SCT}')
print(f'SCTR {SCTR}')
print(f'SCE {SCE}')
print(SCE + SCTR)
print(glt)
print(gle)
print(ftab)
print(f'MCTR {MCTR}')
print(f'MCE {MCE}')
print(f'fcal {fcal}')
print(nj)
print(glt + 1)
print(gle)
print(DHS)

independientes = []

for i in range(len(info['variables'])): 
    for j in range(i + 1, len(info['variables'])): 
        resta = info['variables'][i]['x'] - info['variables'][j]['x']
        if resta > DHS: 
            if info['variables'][i]['nombre'] not in independientes: 
                independientes.append(info['variables'][i]['nombre'])
            if info['variables'][j]['nombre'] not in independientes: 
                independientes.append(info['variables'][j]['nombre'])
            print(resta)

pprint(independientes)

coeficientes_ecuaciones = {}
coeficientes_ecuaciones['legenda'] = [
    {'variable': 'y', 'nombre': x} if independientes.index(x) == 0 else 
    {'variable': f'x{independientes.index(x)}', 'nombre': x} for x in independientes]

coeficientes_ecuaciones['tabla'] = {}
for x in coeficientes_ecuaciones['legenda']: 
    print(x['nombre'])
    coeficientes_ecuaciones['tabla'][x['variable']] = list(filter(lambda y: y['nombre'] == x['nombre'], info['variables']))[0]['lista']
    coeficientes_ecuaciones[f"∑{x['variable']}"] = sum(coeficientes_ecuaciones['tabla'][x['variable']])
    if x['variable'] == 'y': continue
    coeficientes_ecuaciones['tabla'][f"{x['variable']}^2"] = list(filter(lambda y: y['nombre'] == x['nombre'], info['variables']))[0]['lista^2']

for i in range(len(coeficientes_ecuaciones['legenda'])): 
    for j in range(i + 1, len(coeficientes_ecuaciones['legenda'])):
        primer_nombre = coeficientes_ecuaciones['legenda'][i]['variable']
        segundo_nombre = coeficientes_ecuaciones['legenda'][j]['variable']
        primero = coeficientes_ecuaciones['tabla'][primer_nombre]
        segundo = coeficientes_ecuaciones['tabla'][segundo_nombre]
        resultado = [x * y for x, y in zip(primero, segundo)]
        resultado_nombre = f'{primer_nombre}*{segundo_nombre}'
        coeficientes_ecuaciones['tabla'][resultado_nombre] = resultado
        coeficientes_ecuaciones[f'∑{resultado_nombre}'] = sum(coeficientes_ecuaciones['tabla'][resultado_nombre])

pprint(coeficientes_ecuaciones)

x_solas = list(filter(lambda x: ('*' not in x) and ('^' not in x) and ('y' not in x) and ('∑' in x), coeficientes_ecuaciones.keys()))

x_solas_arreglados = [f'∑x{y}' for y in sorted([int(x[2:]) for x in x_solas])]
x_solas_arreglados.insert(0, '∑y')

matriz = []
lista = []

print(x_solas_arreglados)