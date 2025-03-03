import csv
import pandas as pd
from scipy.stats import f, studentized_range
from pprint import pprint
from fractions import Fraction
import math
import seaborn as sns
import matplotlib.pyplot as plt
from fastapi import FastAPI, Request
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

app = FastAPI()

app.mount("/static", StaticFiles(directory="./static"), name="static")

templates = Jinja2Templates(directory="./templates")

info = {'∑xt': 0, '∑xt^2': 0, '(∑xt)^2': 0, 'nt': 0, '(∑xt)^2/nt': 0, 'variables': []}

independientes = []

pedazos_independientes = []
tukey = []
pares = []
impares = []
correlaciones = []
coeficientes_ecuaciones = {}
matriz_resultado = []
resultados_finales = []
matriz = []
resultados = []
x_solas = []
x_solas_arreglados = []
textos = {}

with open('data.csv', newline='', encoding='utf-8') as data: 
    datos = csv.DictReader(data)
    nombres = list(datos.fieldnames)
    for variable in nombres: 
        info['variables'].append({'nombre': variable, 'lista': []})
    [[[info['variables'][i]['lista'].append(float(x[nombres[i]]))] for i in range(len(nombres))] for x in datos]

def info_correlaciones(r): 
    if r == 1: return 'La correlación es perfecta directa'
    elif 1 > r > 0.75: return 'La correlación es fuerte directa'
    elif 0.75 > r > 0.25: return 'La correlación es intermedia directa'
    elif 0.25 > r > 0: return 'La correlación es débil directa'
    elif r == 0: return 'No hay relación'
    elif 0 > r > -0.25: return 'La correlación es débil indirecta'
    elif -0.25 > r > -0.75: return 'La correlación es intermedia indirecta'
    elif -0.75 > r > -1: return 'La correlación es fuerte indirecta'
    elif r == -1: return 'La correlación es perfecta indirecta'
    else: return 'HUBO UN ERROR GRAVE EN EL PROCESO DE LAS FÓRMULAS DE LA CORRELACIÓN'

def mostrar_info(tabla : dict[list[dict] | float]): 
    textos['contingencia'] = tabla
    tamaño = len(tabla['variables'][0]['lista'])
    nombres = [variable["nombre"] for variable in tabla['variables']]
    tamaños = []
    for lista in tabla['variables']: 
        totaluno = max([len(f'{round(x, 4)}') for x in lista['lista']])
        totaldos = max([len(f'{round(x, 4)}') for x in lista['lista^2']])
        tamaños.append((totaluno, totaldos))
    titulo = "| "
    subtitulo = "| "
    for nombre, pedazos in zip(nombres, tamaños): 
        titulo += f'{nombre} ' + (" " * ((pedazos[0] + pedazos[1]) - len(nombre))) + ' | '
        subtitulo += f"x " + (" " * (pedazos[0] - len('x'))) + " | "
        subtitulo += f"x^2 " + (" " * (pedazos[1] - len('x^2'))) + " | "
    print('-' * len(titulo))
    print(titulo)
    print('-' * len(subtitulo))
    print(subtitulo)
    print('-' * len(subtitulo))
    for i in range(tamaño):
        fila = '| '
        for lista, pedazos in zip(tabla['variables'], tamaños): 
            primero = f"{round(lista['lista'][i], 4)}"
            fila += f"{primero} " + (" " * (pedazos[0] - len(primero))) + " | "
            segundo = f"{round(lista['lista^2'][i], 4)}"
            fila += f"{segundo} " + (" " * (pedazos[1] - len(segundo))) + ' | '
        print(fila)
        print('-' * len(fila))
    sumatorias = list(filter(lambda x: '∑' in x or x == 'nt', tabla.keys()))
    for sumas in sumatorias: print(f'{sumas} = {round(tabla[sumas], 4)}')

def mostrar_ecuaciones(tabla: dict[list | float]): 
    print('\nCUADRO DE LOS CÁLCULOS PARA LA RECTA DE REGRESIÓN LINEAL\n')
    for nombre in tabla['legenda']: 
        print(f'{nombre["nombre"]} = {nombre["variable"]}')
    titulos = list(tabla['tabla'].keys())
    tamaño = len(tabla['tabla'][titulos[0]])
    tamaños = []
    for nombre in titulos: 
        total = max([len(f'{round(x, 4)}') for x in tabla['tabla'][nombre]])
        tamaños.append((total))
    titulo = '| '
    for nombre, pedazo in zip(titulos, tamaños): 
        titulo += f'{nombre} ' + (" " * (pedazo - len(nombre))) + " | "
    print("-" * len(titulo))
    print(titulo)
    print("-" * len(titulo))
    for i in range(tamaño): 
        fila = "| "
        for nombre, pedazo in zip(titulos, tamaños): 
            dato = f"{round(tabla['tabla'][nombre][i], 4)}"
            fila += f"{dato} " + (" " * (pedazo - len(dato))) + " | "
        print(fila)
        print('-' * len(fila))
    sumas = list(filter(lambda x: '∑' in x, tabla.keys()))
    for suma in sumas: print(f'{suma} = {tabla[suma]}')
    print('\n')
    

def mostrar_tukey(tabla: list[dict]): 
    textos['tukey'] = tabla
    print('\nTABLA DE TUKEY\n')
    tamaño = len(info['variables'])
    variables = ['y' if i == 0 else f'x{i}' for i in range(tamaño)]
    tamaños = []
    for variable in variables: 
        try: 
            tamaños.append(max([len(str(round(dato['resultado'], 4))) for dato in list(filter(lambda x: x['y'] == variable, tabla))]))
        except: tamaños.append(1)
    maximo = max([len(x) for x in variables])
    titulo = '| ' + ('/' * maximo) + ' | ' 
    for variable, pedazo in zip(variables, tamaños): 
        titulo += f'{variable} ' + (' ' * (pedazo - len(variable))) + '| '
    print('-' * len(titulo))
    print(titulo)
    print('-' * len(titulo))
    for i in range(tamaño):
        fila = f'| {variables[i]} ' + (' ' * (maximo - len(variables[i]))) + '| '
        for variable, pedazo in zip(variables, tamaños): 
            try: 
                numero = [dato['resultado'] for dato in list(filter(lambda x: x['y'] == variable, tabla))][i]
                fila += f'{round(numero, 4)} ' + (' ' * (pedazo - len(str(round(numero, 4))))) + '| '
            except: 
                fila += '/' + ('/' * (pedazo - len(variable))) + '/ | '
        print(fila)
        print('-' * len(fila))
        
def mostrar_correlacion(tabla: dict): 
    print(f'\nCorrelaciones entre {tabla["variables"][0]} y {tabla["variables"][1]}')
    print(f'{tabla["variables"][0]} = x\n{tabla["variables"][1]} = y')
    tamaño = len(tabla['x'])
    lista = ['x', 'y', 'x^2', 'y^2', 'x*y']
    tamaños = [max([len(str(round(y, 4))) for y in tabla[x]]) for x in lista]
    titulo = '| '
    for variable, pedazo in zip(lista, tamaños): 
        titulo += f'{variable} ' + (' ' * (pedazo - len(variable))) + ' | '
    print('-' * len(titulo))
    print(titulo)
    print('-' * len(titulo))
    for i in range(tamaño): 
        fila = '| '
        for variable, pedazo in zip(lista, tamaños): 
            fila += f'{round(tabla[variable][i], 4)} ' + (' ' * (pedazo - len(str(round(tabla[variable][i], 4))))) + ' | '
        print(fila)
        print('-' * len(fila))

    print(f'r = {round(tabla["r"], 4)}')
    tabla['correlacion'] = info_correlaciones(tabla['r'])
    print(tabla['correlacion'])
    sns.set_theme(style="darkgrid")
    datos = pd.DataFrame({tabla['variables'][0]: tabla['x'], tabla['variables'][1]: tabla['y']})
    sns.jointplot(x=tabla['variables'][0], y=tabla['variables'][1], data=datos,
                    kind="reg", truncate=False,
                    color="m")
    lugar = f'./static/{tabla["variables"][0]}-{tabla["variables"][1]}.png'
    plt.savefig(lugar)
    tabla['lugar'] = lugar
    # plt.show()

def mostrar_matriz(tabla: list[list], lista: list): 
    print('\n')
    tamaños = []
    for i in range(len(tabla)):
        tamaños.append(max([len(str(round(x[i], 4))) for x in tabla]))
    tamaño_resultado = max([len(str(round(x, 4))) for x in lista])
    for i in range(len(tabla)): 
        fila = '| '
        for j in range(len(tabla)): 
            numero = str(round(tabla[i][j], 4))
            fila += f'{numero} ' + (' ' * (tamaños[j] - len(numero))) + ' | '
        fila += f'{round(lista[i], 4)}' + (' ' * (tamaño_resultado - len(str(round(lista[i], 4))))) + ' | '
        print('-' * len(fila))
        print(fila)
    print('-' * len(fila) + '\n')

def crear_contingencia(): 
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
    textos['C'] = C
    textos['SCT'] = SCT
    textos['SCTR'] = SCTR
    textos['SCE'] = SCE

    glt = len(info['variables']) - 1
    gle = info['nt'] - len(info['variables'])
    textos['glt'] = glt
    textos['gle'] = gle

    if (glt + gle) != (info['nt'] - 1): print('ALGO ESTÁ MUY MALLLLLLLLLLLLLLLLLLL')

    ftab = f.ppf(q=0.95, dfn=glt, dfd=gle)
    textos['ftab'] = ftab

    MCTR = SCTR / glt
    MCE = SCE / gle
    textos['MCTR'] = MCTR
    textos['MCE'] = MCE

    nj = len(info['variables'][0]['lista'])
    textos['nj'] = nj

    fcal = MCTR / MCE
    textos['fcal'] = fcal

    DHS = (studentized_range.ppf(q=0.95, k=glt + 1, df=gle)) * (math.sqrt(MCE / nj))
    textos['DHS'] = DHS

    # pprint(info)
    mostrar_info(info)
    print(f'C = {round(C, 4)}')
    print(f'SCT = {round(SCT, 4)}')
    print(f'SCTR = {round(SCTR, 4)}')
    print(f'SCE = {round(SCE, 4)}')
    print(f'SCE + SCTR = {round(SCE + SCTR, 4)}')
    print(f'glt = {round(glt, 4)}')
    print(f'gle = {round(gle, 4)}')
    print(f'ftab = {round(ftab, 4)}')
    print(f'MCTR = {round(MCTR, 4)}')
    print(f'MCE = {round(MCE, 4)}')
    print(f'fcal = {round(fcal, 4)}')
    print(f'nj = {round(nj, 4)}')
    print(f'glt + 1 = {glt + 1}')
    print(f'DHS = {round(DHS, 4)}')
    return DHS, fcal, ftab, nj

def crear_tukey(DHS : float, fcal : float, ftab : float): 
    for i in range(len(info['variables'])): 
        for j in range(i + 1, len(info['variables'])): 
            resta = info['variables'][i]['x'] - info['variables'][j]['x']
            if resta > DHS: 
                if info['variables'][i]['nombre'] not in pares: 
                    pares.append(info['variables'][i]['nombre'])
                if info['variables'][j]['nombre'] not in impares: 
                    impares.append(info['variables'][j]['nombre'])
                # print("info['variables'][i]['nombre']round(resta, 4))
                pedazos_independientes.append({
                    'x': info['variables'][i]['nombre'], 
                    'y': info['variables'][j]['nombre']})
            tukey.append({
                'x': 'y' if i == 0 else f'x{i}', 
                'y': 'y' if j == 0 else f'x{j}', 
                'resultado': resta})
    # independientes = ['primero', 'segundo', 'tercero']
    independientes.extend(filter(lambda x: x not in impares, pares))
    independientes.extend(filter(lambda x: x not in pares, impares))

    hipotesis = ''
    if fcal > ftab:
        hipotesis = f'''
    Se rechaza la hipótesis nula (H0), ya que fcal ({round(fcal, 4)}) > ftab ({round(ftab, 4)}), 
    lo cual quiere decir que se acepta la hipótesis alternativa (Ha), 
    lo cual implica que si existe diferencia significativa entre las variables 
    ''' 
    elif ftab > fcal: 
        hipotesis = f'''
    Se acepta la hipótesis nula (H0), ya que fcal ({round(fcal, 4)}) < ftab ({round(ftab, 4)}), 
    lo cual implica que no existen diferencias significativas entre las variables 
    '''
    for i in range(len(independientes)): 
        if i == 0: hipotesis += f'{independientes[i]}'
        elif i == len(independientes) - 1: hipotesis += f' y {independientes[i]} '
        else: hipotesis += f', {independientes[i]}'
    hipotesis += 'sobre las variables '

    nombres = [x['nombre'] for x in info['variables']]

    dependientes = list(filter(lambda x: x not in independientes, nombres))

    for i in range(len(dependientes)): 
        if i == 0: hipotesis += f'{dependientes[i]}'
        elif i == len(dependientes) - 1: hipotesis += f' y {dependientes[i]} '
        else: hipotesis += f', {dependientes[i]}'

    textos['hipotesis'] = hipotesis
    print(hipotesis)

    mostrar_tukey(tukey) 

    pprint(independientes)
    pprint(pedazos_independientes)
    if len(independientes) == 0: return False
    return True

def crear_correlaciones(): 
    for independiente in pedazos_independientes: 
        lista = {'variables': [independiente['x'], independiente['y']]}
        for variable in list(independiente.keys()): 
            data = list(filter(lambda x: x['nombre'] == independiente[variable], info['variables']))[0]
            lista[variable] = data['lista']
            lista[f'∑{variable}'] = sum(lista[variable])
            lista[f'{variable}^2'] = data['lista^2']
            lista[f'∑{variable}^2'] = sum(lista[f'{variable}^2'])
        lista['x*y'] = [x * y for x, y in zip(lista['x'], lista['y'])]
        lista['∑x*y'] = sum(lista['x*y'])
        n = len(lista['y'])
        lista['r'] = (lista['∑x*y'] - ((lista['∑x'] * lista['∑y']) / n)) / (
            math.sqrt((lista['∑x^2'] - ((lista['∑x'] ** 2) / n)) * (lista['∑y^2'] - ((lista['∑y'] ** 2) / n))))
        lista['b'] = (lista['∑x'] * lista['∑y'] - n * lista['∑x*y']) / ((lista['∑x'] ** 2) - n * lista['∑x^2'])
        lista['a'] = (lista['∑y'] - lista['∑x'] * lista['b']) / n
        print(f'ÿ = {round(lista["a"], 4)} + {round(lista["b"], 4)} * X')
        lista['ÿ'] = lista['a'] / (-lista['b'])
        mostrar_correlacion(lista)
        correlaciones.append(lista)

# pprint(correlaciones)

def crear_tablas_ecuaciones(): 
    coeficientes_ecuaciones['legenda'] = [
        {'variable': 'y', 'nombre': x} if independientes.index(x) == 0 else 
        {'variable': f'x{independientes.index(x)}', 'nombre': x} for x in independientes]

    coeficientes_ecuaciones['tabla'] = {}
    for x in coeficientes_ecuaciones['legenda']: 
        # print(f"{x['nombre']} = {x['variable']}")
        coeficientes_ecuaciones['tabla'][x['variable']] = list(filter(lambda y: y['nombre'] == x['nombre'], info['variables']))[0]['lista']
        coeficientes_ecuaciones[f"∑{x['variable']}"] = sum(coeficientes_ecuaciones['tabla'][x['variable']])
        if x['variable'] == 'y': continue
        coeficientes_ecuaciones['tabla'][f"{x['variable']}^2"] = list(filter(lambda y: y['nombre'] == x['nombre'], info['variables']))[0]['lista^2']
        coeficientes_ecuaciones[f"∑{x['variable']}^2"] = sum(coeficientes_ecuaciones['tabla'][f"{x['variable']}^2"])

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

    # pprint(coeficientes_ecuaciones)
    mostrar_ecuaciones(coeficientes_ecuaciones)

def crear_matrices(nj): 
    x_solas = list(filter(lambda x: ('*' not in x) and ('^' not in x) and ('y' not in x) and ('∑' in x), coeficientes_ecuaciones.keys()))

    x_solas_arreglados = [f'∑x{y}' for y in sorted([int(x[2:]) for x in x_solas])]
    x_solas_arreglados.insert(0, '∑y')

    # matriz = []

    lista = [nj if x == '∑y' else coeficientes_ecuaciones[x] for x in x_solas_arreglados]
    matriz.append(lista)

    # resultados = []
    resultados.append(coeficientes_ecuaciones['∑y'])
    for i in range(1, len(x_solas_arreglados)): 
        lista = []
        lista.append(coeficientes_ecuaciones[x_solas_arreglados[i]])
        for j in range(1, len(x_solas_arreglados)): 
            if j == i: 
                lista.append(coeficientes_ecuaciones[f"{x_solas_arreglados[i]}^2"])
                continue
            llave = list(filter(lambda x: ('*' in x) and (x_solas_arreglados[i][1:] in x) and (x_solas_arreglados[j][1:] in x), coeficientes_ecuaciones.keys()))[0]
            lista.append(coeficientes_ecuaciones[llave])
        matriz.append(lista)
        llave = list(filter(lambda x: ('*' in x) and (x_solas_arreglados[i][1:] in x) and ('y' in x), coeficientes_ecuaciones.keys()))[0]
        resultados.append(coeficientes_ecuaciones[llave])

    global matriz_resultado
    matriz_resultado = matriz.copy()
    global resultados_finales
    resultados_finales = resultados.copy()

    for i in range(len(matriz_resultado)): 
        invertido = 1 / matriz_resultado[i][i]
        matriz_resultado[i] = [invertido * x for x in matriz_resultado[i]]
        resultados_finales[i] *= invertido
        for j in range(len(matriz_resultado)): 
            if j == i: continue
            negativo = -matriz_resultado[j][i]
            matriz_resultado[j] = [y + negativo * x for x, y in zip(matriz_resultado[i], matriz_resultado[j])]
            resultados_finales[j] = resultados_finales[j] + negativo * resultados_finales[i]

    pprint(x_solas_arreglados)
    mostrar_matriz(matriz, resultados)
    mostrar_matriz(matriz_resultado, resultados_finales)
    print('RESULTADOS FINALES')
    textos['fracciones_finales'] = []
    for numero in resultados_finales: 
        textos['fracciones_finales'].append(Fraction(numero).limit_denominator())
        print(f'B{resultados_finales.index(numero)} = {Fraction(numero).limit_denominator()}')
        print(f'B{resultados_finales.index(numero)} = {round(numero, 4)}')
    print(sum(resultados_finales))
    # if sum([y if matriz[0].index(x) == 0 else x * y for x, y in zip(matriz[0], resultados_finales)]) == resultados[0]: 
    #     print('todo está bien')
    # else: print('TODO ESTÁ MUY MAAAAAAAL')
    textos['sumatorias_finales'] = []
    for x in matriz: 
        print('//////////////////////////')
        data = [numero * constante for numero, constante in zip(x, resultados_finales)]
        y = 'ÿ = '
        for i in range(len(x)): 
            if i == 0: 
                y += f'{round(x[i], 4)} * B{i}'
                continue
            y += f' + {round(x[i], 4)} * B{i}'
        print(y)
        # print(data)
        resultado = sum(data)
        textos['sumatorias_finales'].append(resultado)
        print(f'ÿ = {round(resultado, 4)}')
        indice = matriz.index(x)
        print(f'∑{'y' if indice == 0 else f'y*x{indice}'} = {round(resultados[indice], 4)}')

def main(): 
    try: 
        dhs, fcal, ftab, nj = crear_contingencia()
        textos['verdad'] = crear_tukey(dhs, fcal, ftab)
        if textos['verdad']: 
            crear_correlaciones()
            crear_tablas_ecuaciones()
            crear_matrices(nj)
        else: print('NO SE PUDO APLICAR TUKEY')
    except: print('NO SE PUDO APLICAR TUKEY')

main()

@app.get('/')
def mostrar(request : Request):
    return templates.TemplateResponse('index.html', {
        'request': request, 'hipotesis': textos['hipotesis'], 'DHS': textos['DHS'], 'MCE': textos['MCE'], 
        'C': textos['C'], 'nj': textos['nj'], 'fcal': textos['fcal'], 'ftab': textos['ftab'], 
        'SCT': textos['SCT'], 'SCTR': textos['SCTR'], 'SCE': textos['SCE'], 'MCTR': textos['MCTR'], 
        'glt': textos['glt'], 'gle': textos['gle'],  'round': round, 'verdad': textos['verdad'],
        'contingencia': textos['contingencia'], 'tukey': textos['tukey'], 'matriz': matriz, 
        'resultados': resultados, 'matriz_resultado': matriz_resultado, 'sumatorias_finales': textos['sumatorias_finales'], 
        'resultados_finales': resultados_finales, 'fracciones_finales': textos['fracciones_finales'], 
        'correlaciones': correlaciones, 'regresiones': coeficientes_ecuaciones})


# if __name__ == '__main__': 
#     main()