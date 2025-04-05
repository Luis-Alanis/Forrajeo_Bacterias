from copy import copy
from multiprocessing import Manager, Pool
import time
from bacteria import bacteria
import numpy as np  # Cambiar import numpy a import numpy as np
import copy

from fastaReader import fastaReader

import time
import pandas as pd
import matplotlib.pyplot as plt
import psutil
import os
import gc


def get_memory_parameters():
    """
    Determina parámetros óptimos según la memoria disponible del sistema.
    Ajusta tamaños de población y otros parámetros para evitar errores de memoria.
    """
    try:
        # Obtener memoria disponible en GB
        mem = psutil.virtual_memory()
        available_memory_gb = mem.available / (1024 ** 3)
        
        # Establecer parámetros basados en memoria disponible
        if available_memory_gb > 8:
            # Máquina con mucha memoria
            return {
                "numeroDeBacterias": 4,
                "numRandomBacteria": 2,
                "iteraciones": 8, 
                "tumbo": 3,
                "nado": 2,
                "chunk_size_factor": 6,  # Divisor para chunks
                "max_processes": min(8, os.cpu_count()),
                "memory_status": "high"
            }
        elif available_memory_gb > 4:
            # Memoria media
            return {
                "numeroDeBacterias": 3,
                "numRandomBacteria": 1,
                "iteraciones": 6,
                "tumbo": 2,
                "nado": 2,
                "chunk_size_factor": 8,
                "max_processes": min(4, os.cpu_count()),
                "memory_status": "medium"
            }
        else:
            # Poca memoria disponible
            return {
                "numeroDeBacterias": 2,
                "numRandomBacteria": 1,
                "iteraciones": 5,
                "tumbo": 2,
                "nado": 1,
                "chunk_size_factor": 10,
                "max_processes": min(2, os.cpu_count()),
                "memory_status": "low"
            }
    except:
        # Valores por defecto si no se puede acceder a la información de memoria
        return {
            "numeroDeBacterias": 2,
            "numRandomBacteria": 1,
            "iteraciones": 5,
            "tumbo": 2,
            "nado": 1,
            "chunk_size_factor": 10,
            "max_processes": 2,
            "memory_status": "unknown"
        }


if __name__ == "__main__":
    # Obtener parámetros optimizados para la memoria disponible
    params = get_memory_parameters()
    
    # Usar los parámetros ajustados
    numeroDeBacterias = params["numeroDeBacterias"]
    numRandomBacteria = params["numRandomBacteria"]
    iteraciones = params["iteraciones"]
    tumbo = params["tumbo"]
    nado = params["nado"]
    
    print(f"Memoria disponible: {params['memory_status']}")
    print(f"Configuración: Bacterias={numeroDeBacterias}, Iteraciones={iteraciones}")
    
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
                # Clonar la secuencia para evitar referencias cruzadas
                seq_clone = list(secuencias[j])
                bacterium.append(seq_clone)
            poblacion[i] = bacterium  # Almacenar directamente como lista, no tupla
           
   




    def printPoblacion():
        for i in range(numeroDeBacterias):
            print(poblacion[i])
            
    

    #---------------------------------------------------------------------------------------------------------
    operadorBacterial = bacteria(numeroDeBacterias, memory_params={
        "chunk_size_factor": params["chunk_size_factor"],
        "max_processes": params["max_processes"]
    })    
    veryBest = [None, None, None] #indice, fitness, secuencias
    
    #registra el tiempo de inicio
    start_time = time.time()
    
    print("poblacion inicial ...")

    # Inicializar el archivo con encabezados
    with open("resultadosMejorados.csv", "w") as f:
        f.write("corrida,fitness,tiempo_ejecucion,interaccion,calificacion_blosum,memoria_usado_pct\n")

    resultados = []
    for corrida in range(30):
        print(f"Corrida numero: {corrida + 1}")
        start_time = time.time()
        
        # Monitoreo de memoria
        mem_start = psutil.virtual_memory().percent

        poblacionInicial() 
        
        for it in range(iteraciones):
            print(f"Iteración {it+1}/{iteraciones}")
            
            # Verificar memoria disponible en cada iteración y ajustar estrategia
            mem_available = psutil.virtual_memory().available / (1024 ** 3)
            current_memory_pct = psutil.virtual_memory().percent
            
            print(f"Memoria disponible: {mem_available:.2f} GB ({current_memory_pct:.1f}% usado)")
            
            # Si estamos por debajo del umbral crítico, forzar liberación de memoria
            if mem_available < 0.75:  # Menos de 750MB disponibles
                print("¡Memoria baja detectada! Optimizando uso de recursos...")
                gc.collect()
                # Reducir temporalmente parámetros de procesamiento
                operadorBacterial.memory_params["chunk_size_factor"] += 5
                operadorBacterial.memory_params["max_processes"] = 1
            else:
                # Restaurar parámetros originales si hay suficiente memoria
                operadorBacterial.memory_params["chunk_size_factor"] = params["chunk_size_factor"]
                operadorBacterial.memory_params["max_processes"] = params["max_processes"]
            
            print("poblacion inicial creada - Tumbo ...")
            operadorBacterial.tumbo(numSec, poblacion, tumbo)
            print("Tumbo Realizado - Cuadrando ...")
            operadorBacterial.cuadra(numSec, poblacion)
            print("poblacion inicial cuadrada - Creando granLista de Pares...")
            
            # Liberar memoria antes de operaciones intensivas
            gc.collect()
            
            operadorBacterial.creaGranListaPares(poblacion)
            print("granList: creada - Evaluando Blosum Parallel")
            # Usar la versión optimizada para reducir uso de memoria
            try:
                operadorBacterial.optimizedEvaluaBlosum()
                print("blosum evaluado - creando Tablas Atract Parallel...")
            except Exception as e:
                print(f"Error en evaluación BLOSUM: {e}")
                # Intentar con valores más conservadores
                operadorBacterial.memory_params["chunk_size_factor"] = 20
                operadorBacterial.memory_params["max_processes"] = 1
                print("Reintentando con parámetros más conservadores...")
                operadorBacterial.optimizedEvaluaBlosum()
            
            # Liberar memoria después de operaciones intensivas
            gc.collect()
            
            try:
                operadorBacterial.creaTablasAtractRepel(poblacion, dAttr, wAttr, dRepel, wRep)
                operadorBacterial.creaTablaInteraction()
                print("tabla Interaction creada - creando tabla Fitness")
                # Usar la versión mejorada de fitness
                operadorBacterial.creaTablaFitnessEnhanced()
            except Exception as e:
                print(f"Error en cálculo de fitness: {e}")
                # Plan de contingencia: usar la versión simple si la mejorada falla
                operadorBacterial.creaTablaFitness()
            
            print("tabla Fitness creada ")
            globalNFE += operadorBacterial.getNFE()
            
            try:
                bestIdx, bestFitness = operadorBacterial.obtieneBest(globalNFE)
                
                # Aplicar diferentes estrategias de optimización en iteraciones alternadas
                if it % 3 == 0 and it > 0:
                    print("Aplicando optimización inteligente de gaps...")
                    operadorBacterial.smartGapOptimization(poblacion, bestIdx)
                    # Re-calcular mejor después de optimización
                    bestIdx, bestFitness = operadorBacterial.obtieneBest(globalNFE)
                elif it % 3 == 1 and it > 0:
                    print("Aplicando alineamiento mejorado...")
                    operadorBacterial.enhancedAlignment(poblacion, bestIdx)
                    # Re-calcular mejor después de optimización
                    bestIdx, bestFitness = operadorBacterial.obtieneBest(globalNFE)
                
                if (veryBest[0] is None) or (bestFitness > veryBest[1]):
                    veryBest[0] = bestIdx
                    veryBest[1] = bestFitness
                    veryBest[2] = copy.deepcopy(poblacion[bestIdx])
                
                operadorBacterial.replaceWorst(poblacion, veryBest[0])
                
                # Aplicar diversificación para mantener diversidad
                print("Aplicando diversificación de población...")
                poblacion = operadorBacterial.diversification(poblacion, bestIdx, it, iteraciones)
            except Exception as e:
                print(f"Error en optimización: {e}")
                # Mantener el mejor encontrado hasta ahora
                if veryBest[0] is not None:
                    operadorBacterial.replaceWorst(poblacion, veryBest[0])
            
            # Liberar memoria explícitamente
            if it % 2 == 1 or mem_available < 1.0:
                gc.collect()
            
            print(f"Best Fitness: {bestFitness:.2f}, Global NFE: {globalNFE}")
    
        # Registrar uso de memoria final
        mem_end = psutil.virtual_memory().percent
        execution_time = time.time() - start_time
        
        # Crear un DataFrame solo con el resultado de esta corrida
        df_current = pd.DataFrame({
            "corrida": [corrida + 1],
            "fitness": [veryBest[1]],
            "tiempo_ejecucion": [execution_time],
            "interaccion": [operadorBacterial.tablaInteraction[veryBest[0]]],
            "calificacion_blosum": [operadorBacterial.blosumScore[veryBest[0]]],
            "memoria_usado_pct": [mem_end]
        })
        
        # Agregar al archivo CSV sin encabezados
        df_current.to_csv("resultadosMejorados.csv", mode="a", header=False, index=False)
        
        print(f"Corrida {corrida + 1} completada: Fitness={veryBest[1]:.2f}, Tiempo={execution_time:.2f}s")

    print("Very Best: ", veryBest)
    #imprime el tiempo de ejecucion
    print("--- %s seconds ---" % (time.time() - start_time))

    # Leer los resultados completos
    df = pd.read_csv("resultadosMejorados.csv")

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

    plt.figure(figsize=(10, 6))
    plt.plot(df["corrida"], df["interaccion"], marker="o", label="Interacción")
    plt.xlabel("Corrida")
    plt.ylabel("Interacción")
    plt.title("Interacción por Corrida")
    plt.legend()
    plt.grid()
    plt.savefig("interaccion_por_corrida.png")
    plt.show()

    # Graficar calificación BLOSUM vs corrida
    plt.figure(figsize=(10, 6))
    plt.plot(df["corrida"], df["calificacion_blosum"], marker="o", label="Calificación BLOSUM")
    plt.xlabel("Corrida")
    plt.ylabel("Calificación BLOSUM")
    plt.title("Calificación BLOSUM por Corrida")
    plt.legend()
    plt.grid()
    plt.savefig("calificacion_blosum_por_corrida.png")
    plt.show()