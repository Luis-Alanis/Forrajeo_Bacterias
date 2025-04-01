from copy import copy
from multiprocessing import Manager, Pool
import time
from bacteria import bacteria
import numpy
import copy

from fastaReader import fastaReader

import time
import pandas as pd
import matplotlib.pyplot as plt


if __name__ == "__main__":
    numeroDeBacterias = 2  # Número moderado de bacterias para evitar un uso excesivo de memoria.
    numRandomBacteria = 1   # Algunas bacterias aleatorias para diversificar la búsqueda.
    iteraciones = 5        # Número de iteraciones razonable para obtener buenos resultados sin sobrecargar el sistema.
    tumbo = 2               # Número de gaps a insertar, ajustado para un balance entre complejidad y rendimiento.
    nado = 2                # Mantén este valor bajo para evitar cálculos excesivos.
    secuencias = list()
    
    secuencias = fastaReader().seqs
    names = fastaReader().names
    
        
    
  
         
    
    
    #hace todas las secuencias listas de caracteres
    for i in range(len(secuencias)):
        #elimina saltos de linea
        secuencias[i] = list(secuencias[i])
        

    

    globalNFE = 0                            #numero de evaluaciones de la funcion objetivo
    
    

    dAttr = 0.2  # Distancia de atracción.
    wAttr = 0.1  # Peso de atracción, ajustado para un comportamiento estable.
    dRepel = 0.1 # Distancia de repulsión.
    wRep = 0.05  # Peso de repulsión, ajustado para evitar estancamiento.
    
   

  
    
    manager = Manager()
    numSec = len(secuencias)
    print("numSec: ", numSec)
    
    poblacion = manager.list(range(numeroDeBacterias))
    names = manager.list(names)
    NFE = manager.list(range(numeroDeBacterias))
    
    
    # print(secuencias)



    def poblacionInicial():    #lineal
        #crece la poblacion al numero de bacterias
        for i in range(numeroDeBacterias):
            bacterium = []
            for j in range(numSec):
                bacterium.append(secuencias[j])
            poblacion[i] = list(bacterium)
           
   




    def printPoblacion():
        for i in range(numeroDeBacterias):
            print(poblacion[i])
            
    

    #---------------------------------------------------------------------------------------------------------
    operadorBacterial = bacteria(numeroDeBacterias)    
    veryBest = [None, None, None] #indice, fitness, secuencias
    
    #registra el tiempo de inicio
    start_time = time.time()
    
    print("poblacion inicial ...")

    resultados = []
    for corrida in range(30):
        start_time = time.time()
        operadorBacterial = bacteria(numeroDeBacterias)
        veryBest = [None, None, None]  # índice, fitness, secuencias
        globalNFE = 0

        poblacionInicial() 
        
        for it in range(iteraciones):
            print("poblacion inicial creada - Tumbo ...")
            operadorBacterial.tumbo(numSec, poblacion, tumbo)
            print("Tumbo Realizado - Cuadrando ...")
            operadorBacterial.cuadra(numSec, poblacion)
            print("poblacion inicial cuadrada - Creando granLista de Pares...")
            operadorBacterial.creaGranListaPares(poblacion)
            print("granList: creada - Evaluando Blosum Parallel")
            operadorBacterial.evaluaBlosum()  #paralelo
            print("blosum evaluado - creando Tablas Atract Parallel...")

            operadorBacterial.creaTablasAtractRepel(poblacion, dAttr, wAttr, dRepel, wRep)


            operadorBacterial.creaTablaInteraction()
            print("tabla Interaction creada - creando tabla Fitness")
            operadorBacterial.creaTablaFitness()
            print("tabla Fitness creada ")
            globalNFE += operadorBacterial.getNFE()
            bestIdx, bestFitness = operadorBacterial.obtieneBest(globalNFE)
            if (veryBest[0] == None) or (bestFitness > veryBest[1]): #Remplaza el mejor 
                veryBest[0] = bestIdx
                veryBest[1] = bestFitness
                veryBest[2] = copy.deepcopy(poblacion[bestIdx])
            operadorBacterial.replaceWorst(poblacion, veryBest[0])
            operadorBacterial.resetListas(numeroDeBacterias)
    
    execution_time = time.time() - start_time
    resultados.append({
        "corrida": corrida + 1,
        "fitness": veryBest[1],
        "tiempo_ejecucion": execution_time,
        "interaccion": operadorBacterial.tablaInteraction[veryBest[0]],
        "calificacion_blosum": operadorBacterial.blosumScore[veryBest[0]]
    })
    df_resultados = pd.DataFrame(resultados)
    df_resultados.to_csv("resultados.csv", index=False)

    print("Very Best: ", veryBest)
    #imprime el tiempo de ejecucion
    print("--- %s seconds ---" % (time.time() - start_time))

    # Leer los resultados
    df = pd.read_csv("resultados.csv")

    # Mostrar tabla resumen
    print(df.describe())

    # Graficar fitness vs corrida
    plt.figure(figsize=(10, 6))
    plt.plot(df["corrida"], df["fitness"], marker="o", label="Fitness")
    plt.xlabel("Corrida")
    plt.ylabel("Fitness")
    plt.title("Fitness por Corrida")
    plt.legend()
    plt.grid()
    plt.savefig("fitness_por_corrida.png")
    plt.show()

    # Graficar tiempo de ejecución vs corrida
    plt.figure(figsize=(10, 6))
    plt.plot(df["corrida"], df["tiempo_ejecucion"], marker="o", label="Tiempo de Ejecución")
    plt.xlabel("Corrida")
    plt.ylabel("Tiempo de Ejecución (s)")
    plt.title("Tiempo de Ejecución por Corrida")
    plt.legend()
    plt.grid()
    plt.savefig("tiempo_por_corrida.png")
    plt.show()