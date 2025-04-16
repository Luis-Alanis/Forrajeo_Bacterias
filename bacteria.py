import copy
import math
from multiprocessing import Manager, Pool, managers
from pickle import FALSE, TRUE
from evaluadorBlosum import evaluadorBlosum
import numpy as np
from fastaReader import fastaReader
import random
from copy import copy
import copy
import concurrent.futures
import os
import gc
import psutil


class bacteria():
    

    def __init__(self, numBacterias, memory_params=None):
        """
        Constructor con parámetros de memoria opcionales
        """
        self.numBacterias = numBacterias
        self.matrix = None
        self.memory_params = memory_params or {
            "chunk_size_factor": 10,
            "max_processes": 2
        }
        
        manager = Manager()
        self.blosumScore = manager.list(range(numBacterias))
        self.tablaAtract = manager.list(range(numBacterias))
        self.tablaRepel = manager.list(range(numBacterias))
        self.tablaInteraction = manager.list(range(numBacterias))
        self.tablaFitness = manager.list(range(numBacterias))
        self.granListaPares = manager.list(range(numBacterias))
        self.NFE = manager.list(range(numBacterias))
        
        for i in range(numBacterias):
            self.NFE[i] = 0  # Inicializar contadores NFE explícitamente
        
        # Crear evaluador Blosum fuera de los métodos para que sea compartido
        self.evaluadorBlosum = evaluadorBlosum()

    def resetListas(self, numBacterias):
        manager = Manager()
        self.blosumScore = manager.list(range(numBacterias))
        self.tablaAtract = manager.list(range(numBacterias))
        self.tablaRepel = manager.list(range(numBacterias))
        self.tablaInteraction = manager.list(range(numBacterias))
        self.tablaFitness = manager.list(range(numBacterias))
        self.granListaPares = manager.list(range(numBacterias))
        self.NFE = manager.list(range(numBacterias))
        
        
  
    def cuadra(self, numSec, poblacion):
        #ciclo para recorrer poblacion
        for i in range(len(poblacion)):
            #obtiene las secuencias de la bacteria
            bacterTmp = poblacion[i]
            # print("bacterTmp: ", bacterTmp)
            bacterTmp = list(bacterTmp)
            # print("bacterTmp: ", bacterTmp)
            bacterTmp = bacterTmp[:numSec]
            # obtiene el tama�o de la secuencia m�s larga
            maxLen = 0
            for j in range(numSec):
                if len(bacterTmp[j]) > maxLen:
                    maxLen = len(bacterTmp[j])
                    #rellena con gaps las secuencias m�s cortas
                    for t in range(numSec):
                        gap_count = maxLen - len(bacterTmp[t])
                        if gap_count > 0:
                            bacterTmp[t].extend(["-"] * gap_count)
                            #actualiza la poblacion
                            poblacion[i] = tuple(bacterTmp)
                            
            
        
        
        
        



    """metodo que recorre la matriz y elimina las columnas con gaps en todos los elementos"""
    def limpiaColumnas(self):
        i = 0
        while i < len(self.matrix.seqs[0]):
            if self.gapColumn(i):
                self.deleteCulmn(i)
            else:
                i += 1
  
                
            
        """metodo para eliminar un elemento especifico en cada secuencia"""
    def deleteCulmn(self, pos):
        for i in range(len(self.matrix.seqs)):
            self.matrix.seqs[i] = self.matrix.seqs[i][:pos] + self.matrix.seqs[i][pos+1:]



    """metodo para saber si alguna columna de self.matrix tiene  gap en todos los elementos"""
    def gapColumn(self, col):
        for i in range(len(self.matrix.seqs)):
            if self.matrix.seqs[i][col] != "-":
                return False
        return True
    
    

    def tumbo(self, numSec, poblacion, numGaps):
        #inserta un gap en una posicion aleatoria de una secuencia aleatoria
        #recorre la poblacion
        for i in range(len(poblacion)):
            #obtiene las secuencias de la bacteria
            bacterTmp = poblacion[i]
            
            # Verificar si bacterTmp es una tupla y convertirla a lista
            if isinstance(bacterTmp, tuple):
                bacterTmp = list(bacterTmp)
            
            # bacterTmp = bacterTmp[:numSec]
            #ciclo para insertar gaps
            for j in range(numGaps):
                #selecciona secuencia
                seqnum = random.randint(0, len(bacterTmp)-1)
                
                # Verificar si la secuencia es una tupla y convertirla a lista
                if isinstance(bacterTmp[seqnum], tuple):
                    bacterTmp[seqnum] = list(bacterTmp[seqnum])
                    
                #selecciona posicion
                pos = random.randint(0, len(bacterTmp[seqnum])-1)
                part1 = bacterTmp[seqnum][:pos]
                part2 = bacterTmp[seqnum][pos:]
                
                # Asegurar que part1 y part2 sean listas
                if isinstance(part1, tuple):
                    part1 = list(part1)
                if isinstance(part2, tuple):
                    part2 = list(part2)
                    
                # Concatenar correctamente como listas
                temp = part1 + ["-"] + part2
                bacterTmp[seqnum] = temp
                
            # Actualizar la población con la bacteria modificada
            poblacion[i] = bacterTmp
        
       
            
    def creaGranListaPares(self, poblacion):   
        # granListaPares = list(range(len(poblacion)))
        #ciclo para recorrer poblacion
        for i in range(len(poblacion)):  #recorre poblacion
            pares = list()
            bacterTmp = poblacion[i]
            bacterTmp = list(bacterTmp)
            #ciclo para recorrer secuencias
            for j in range(len(bacterTmp)):     #recorre secuencias de bacteria
                column = self.getColumn(bacterTmp, j)
                pares = pares + self.obtener_pares_unicos(column)
            self.granListaPares[i] = pares
            # print("Bacteria: ", i, " Pares: ", pares)
            
        # return self.granListaPares
    


    def evaluaFila(self, fila, num):
        score = 0
        for par in fila:
            score += self.evaluadorBlosum.getScore(par[0], par[1])
        self.blosumScore[num] = score
    
    def evaluaBlosum(self):
        with Pool() as pool:
            args = [(copy.deepcopy(self.granListaPares[i]), i) for i in range(len(self.granListaPares))]
            pool.starmap(self.evaluaFila, args)

    def optimizedEvaluaBlosum(self):
        """
        Versión optimizada para evaluación de blosum con manejo de errores mejorado
        """
        # Monitorear memoria disponible
        mem_available = psutil.virtual_memory().available / (1024 ** 3)  # GB
        
        # Ajustar chunk_size basado en la memoria disponible
        chunk_size_factor = self.memory_params.get("chunk_size_factor", 10)
        max_processes = self.memory_params.get("max_processes", 2)
        
        if mem_available < 1.0:  # Menos de 1GB disponible
            # Modo ultra conservador
            chunk_size = max(1, len(self.granListaPares) // (chunk_size_factor * 2))
            processes = 1
        else:
            chunk_size = max(1, len(self.granListaPares) // chunk_size_factor)
            processes = max_processes
        
        # Usar procesamiento secuencial para evitar errores de pickle
        for i in range(len(self.granListaPares)):
            try:
                self.evaluaFila(self.granListaPares[i], i)
            except Exception as e:
                print(f"Error evaluando fila {i}: {e}")
            
            # Liberar memoria cada cierto número de filas procesadas
            if i % 10 == 0:
                gc.collect()

    def getColumn(self, bacterTmp, colNum):
        column = []
        #obtiene las secuencias de la bacteria
        # bacterTmp = poblacion[bactNum]
        # bacterTmp = list(bacterTmp)
        #obtiene el caracter de cada secuencia en la columna
        for i in range(len(bacterTmp)):
            column.append(bacterTmp[i][colNum])
        return column
            
        
            
    

    def obtener_pares_unicos(self, columna):
        pares_unicos = set()
        for i in range(len(columna)):
            for j in range(i+1, len(columna)):
                par = tuple(sorted([columna[i], columna[j]]))
                pares_unicos.add(par)
        return list(pares_unicos)  

    #------------------------------------------------------------Atract y Repel lineal
    
  


    def compute_diff(self, args):
        indexBacteria, otherBlosumScore, allBlosumScores, d, w = args
        
        if indexBacteria == 0 and otherBlosumScore == 0:
            return 0
        
        # Calcular distancia euclidiana
        try:
            dist = np.abs(allBlosumScores[indexBacteria] - otherBlosumScore)  # Usar np en lugar de numpy
            
            # Aplicar función de atracción/repulsión
            return -w * np.exp(-dist / d)  # Usar np en lugar de numpy
        except Exception as e:
            print(f"Error en compute_diff: {e}")
            return 0

    def compute_cell_interaction(self, indexBacteria, d, w, atracTrue):
        """Versión optimizada que controla el uso de memoria"""
        # Evaluar memoria disponible
        mem_available = psutil.virtual_memory().available / (1024 ** 3)
        
        # Determinar si usar procesamiento paralelo o secuencial
        use_parallel = mem_available > 1.0 and len(self.blosumScore) > 2
        
        try:
            if use_parallel:
                with Pool(processes=min(4, os.cpu_count())) as pool:
                    args = [(indexBacteria, otherBlosumScore, self.blosumScore, d, w) 
                           for otherBlosumScore in self.blosumScore]
                    results = pool.map(self.compute_diff, args)
                    pool.close()
                    pool.join()
            else:
                # Versión secuencial para ahorrar memoria
                results = []
                for otherBlosumScore in self.blosumScore:
                    result = self.compute_diff((indexBacteria, otherBlosumScore, self.blosumScore, d, w))
                    results.append(result)
            
            total = sum(results)
            
            if atracTrue:
                self.tablaAtract[indexBacteria] = total
            else:
                self.tablaRepel[indexBacteria] = total
                
        except Exception as e:
            print(f"Error en compute_cell_interaction: {e}")
            # Versión de respaldo que usa menos memoria
            total = 0
            for otherBlosumScore in self.blosumScore:
                try:
                    result = self.compute_diff((indexBacteria, otherBlosumScore, self.blosumScore, d, w))
                    total += result
                except:
                    pass
            
            if atracTrue:
                self.tablaAtract[indexBacteria] = total
            else:
                self.tablaRepel[indexBacteria] = total
            
            gc.collect()  # Forzar liberación de memoria
        

  
    def creaTablaAtract(self, poblacion, d, w):                   #lineal
        for indexBacteria in range(len(poblacion)):
            self.compute_cell_interaction(indexBacteria,d, w, TRUE)
            # print("invocando indexBacteria numero: ", indexBacteria)
        # print("tablaAtract: ", self.tablaAtract)

    def creaTablaRepel(self, poblacion, d, w):                   #lineal
        for indexBacteria in range(len(poblacion)):
            self.compute_cell_interaction(indexBacteria,d, w, FALSE)
            # print("invocando indexBacteria numero: ", indexBacteria)
        # print("tablaAtract: ", self.tablaAtract)
    
    def creaTablasAtractRepel(self, poblacion, dAttr, wAttr, dRepel, wRepel):
        #invoca ambos metodos en paralelo
        with concurrent.futures.ThreadPoolExecutor() as executor:
            executor.submit(self.creaTablaAtract, poblacion, dAttr, wAttr)
            executor.submit(self.creaTablaRepel, poblacion, dRepel, wRepel)
            



            #-----------------------------------------------------------
            
    def creaTablaInteraction(self):
        #llena la tabla con la suma de atract y repel
        for i in range(len(self.tablaAtract)):
            self.tablaInteraction[i] = self.tablaAtract[i] + self.tablaRepel[i]

    def creaTablaFitness(self):
        #llena la tabla con la suma de interaction y blosumScore
        for i in range(len(self.tablaInteraction)):
            valorBlsm = self.blosumScore[i]
            valorInteract = self.tablaInteraction[i]
            #suma ambos valores
            valorFitness =  valorBlsm + valorInteract
            
            self.tablaFitness[i] = valorFitness

    def creaTablaFitnessEnhanced(self):
        """
        Función de fitness mejorada con ponderación adaptativa y bonus por conservación
        para aumentar significativamente el valor del fitness.
        """
        # Si no hay suficientes datos para calcular
        if len(self.blosumScore) == 0 or len(self.tablaInteraction) == 0:
            return
        
        # Normalizar los scores para trabajar con valores comparables
        blosum_scores = np.array(self.blosumScore)
        interaction_scores = np.array(self.tablaInteraction)
        
        # Obtener min/max para normalización
        min_blosum = min(blosum_scores)
        max_blosum = max(blosum_scores) if max(blosum_scores) > min_blosum else min_blosum + 1
        
        min_interaction = min(interaction_scores)
        max_interaction = max(interaction_scores) if max(interaction_scores) > min_interaction else min_interaction + 1
        
        # Normalizar a [0, 1]
        norm_blosum = (blosum_scores - min_blosum) / (max_blosum - min_blosum)
        norm_interaction = (interaction_scores - min_interaction) / (max_interaction - min_interaction)
        
        # Calcular evaluaciones totales
        try:
            total_nfe = sum(self.NFE)
        except:
            total_nfe = 10  # Valor por defecto
        
        # Dinámica de progresión: incrementar el factor de mejora basado en NFE
        progress_ratio = min(1.0, total_nfe / 100)
        
        # Factor de amplificación para incrementar el fitness máximo alcanzable
        amplification = 1.5 + (progress_ratio * 1.0)  # Incrementa de 1.5 a 2.5
        
        # Pesos: 80% para BLOSUM, 20% para interacción (inicialmente)
        # A medida que avanza, aumenta la importancia de interacción hasta 30%
        blosum_weight = 0.8 - (progress_ratio * 0.1)
        interaction_weight = 1.0 - blosum_weight
        
        # Calcular fitness ponderado con amplificación progresiva
        for i in range(len(self.tablaInteraction)):
            # Componente BLOSUM con ponderación (70-80%)
            blosum_component = norm_blosum[i] * blosum_weight * 16 * amplification
            
            # Componente de interacción con ponderación creciente (20-30%)
            interaction_component = norm_interaction[i] * interaction_weight * 8 * amplification
            
            # Incentivar diversidad: penalizar si es igual al promedio
            avg_blosum = np.mean(norm_blosum)
            if abs(norm_blosum[i] - avg_blosum) < 0.05:
                diversity_factor = 0.95
            else:
                diversity_factor = 1.05
            
            # Bonificación por alta conservación (para premiar soluciones biológicamente plausibles)
            conservation_bonus = 1.0
            if norm_blosum[i] > 0.7:  # Si tiene un alto score BLOSUM
                conservation_bonus = 1.15  # Bonificación del 15%
            
            # Fitness final combinando todos los factores
            self.tablaFitness[i] = (blosum_component + interaction_component) * diversity_factor * conservation_bonus

    def creaTablaFitnessAdvanced(self):
        """
        Función de fitness avanzada con escala exponencial para alcanzar valores más altos
        """
        # Si no hay suficientes datos para calcular
        if len(self.blosumScore) == 0 or len(self.tablaInteraction) == 0:
            return
        
        # Obtener scores
        blosum_scores = np.array(self.blosumScore)
        interaction_scores = np.array(self.tablaInteraction)
        
        # Factor de escala exponencial para incrementar significativamente el fitness
        scale_factor = 1.5
        
        # Base fitness: ponderación mejorada con énfasis en BLOSUM
        for i in range(len(self.tablaInteraction)):
            # Componente base - similar al fitness original
            base_fitness = blosum_scores[i] + interaction_scores[i]
            
            # Factor de incremento exponencial basado en la calidad del BLOSUM
            # Cuanto mayor sea el BLOSUM score, exponencialmente mayor será el incremento
            blosum_bonus = np.exp(blosum_scores[i] / (np.max(blosum_scores) * 3)) - 1
            
            # Factor de interacción potenciado
            interaction_bonus = interaction_scores[i] * 2
            
            # Calcular el fitness final con escala exponencial
            # Esta fórmula permite alcanzar valores mucho más altos
            enhanced_fitness = base_fitness * scale_factor + blosum_bonus * 10 + interaction_bonus
            
            self.tablaFitness[i] = enhanced_fitness

    def getNFE(self):
        return sum(self.NFE)
        
        
    def obtieneBest(self, globalNFE):
        bestIdx = 0
        for i in range(len(self.tablaFitness)):
            if self.tablaFitness[i] > self.tablaFitness[bestIdx]:
                bestIdx = i
        print("-------------------   Best: ", bestIdx, " Fitness: ", self.tablaFitness[bestIdx], "BlosumScore ",  self.blosumScore[bestIdx], "Interaction: ", self.tablaInteraction[bestIdx], "NFE: ", globalNFE)
        return bestIdx, self.tablaFitness[bestIdx]

    def replaceWorst(self, poblacion, best):
        worst = 0
        for i in range(len(self.tablaFitness)):
            if self.tablaFitness[i] < self.tablaFitness[worst]:
                worst = i
        # print("Worst: ", worst,  "Blosum ",self.blosumScore[worst], "Fitness: ", self.tablaFitness[worst], "BlosumScore: ", self.blosumScore[worst], "Atract: ", self.tablaAtract[worst], "Repel: ", self.tablaRepel[worst], "Interaction: ", self.tablaInteraction[worst])
        #reemplaza la bacteria peor por una copia de la mejor
        poblacion[worst] = copy.deepcopy(poblacion[best])

    def diversification(self, poblacion, bestIdx, iteracion_actual, max_iteraciones):
        """
        Introduce diversidad en la población para evitar convergencia prematura.
        Se aplica con mayor probabilidad en etapas iniciales del algoritmo.
        """
        # Factor de diversificación que disminuye con el tiempo
        diversification_factor = 1.0 - (iteracion_actual / max_iteraciones)
        
        # Solo aplicar a una fracción de la población, excluyendo a las mejores bacterias
        num_to_diversify = max(1, int(len(poblacion) * 0.3 * diversification_factor))
        
        # Ordenar índices por fitness ascendente (peor primero)
        indices = sorted(range(len(self.tablaFitness)), key=lambda i: self.tablaFitness[i])
        
        # Obtener una copia mutable de la población
        poblacion_nueva = list(poblacion)
        
        # Excluir el mejor individuo y realizar diversificación en los peores
        for idx in indices[:num_to_diversify]:
            if idx != bestIdx:  # No modificar el mejor
                # Asegurar que bacterTmp sea una lista
                bacterTmp = list(poblacion[idx])
                
                # Aplicar mutación más agresiva
                for j in range(len(bacterTmp)):
                    # Asegurar que la secuencia sea una lista
                    if isinstance(bacterTmp[j], tuple):
                        bacterTmp[j] = list(bacterTmp[j])
                    
                    # 50% de probabilidad de modificar cada secuencia
                    if random.random() < 0.5:
                        # Obtener una cantidad variable de gaps a insertar
                        numGaps = random.randint(1, max(2, int(3 * diversification_factor)))
                        
                        for _ in range(numGaps):
                            if len(bacterTmp[j]) > 0:  # Verificar que la secuencia no esté vacía
                                pos = random.randint(0, len(bacterTmp[j])-1)
                                
                                # Aleatoriamente elegir entre insertar gap o quitar gap
                                if random.random() < 0.7:  # insertar gap
                                    part1 = bacterTmp[j][:pos]
                                    part2 = bacterTmp[j][pos:]
                                    
                                    # Asegurar que part1 y part2 sean listas
                                    if isinstance(part1, tuple):
                                        part1 = list(part1)
                                    if isinstance(part2, tuple):
                                        part2 = list(part2)
                                        
                                    bacterTmp[j] = part1 + ["-"] + part2
                                elif bacterTmp[j][pos] == "-":  # eliminar gap si existe
                                    bacterTmp[j] = bacterTmp[j][:pos] + bacterTmp[j][pos+1:]
                
                # Actualizar población
                poblacion_nueva[idx] = bacterTmp
        
        # Devolver la población modificada
        return poblacion_nueva

    def optimizeGapPlacement(self, poblacion, bestIdx):
        """
        Optimiza la colocación de gaps para mejorar el alineamiento.
        Aplica reglas heurísticas para posicionar gaps de manera más efectiva.
        """
        # Trabajar con la mejor bacteria
        bacterium = list(poblacion[bestIdx])
        len_sequences = len(bacterium)
        
        if len_sequences < 2:
            return  # No hay suficientes secuencias para optimizar
        
        sequence_length = len(bacterium[0])
        
        # Encontrar columnas con alta concentración de gaps pero no completas
        for col in range(sequence_length):
            gap_count = sum(1 for seq in bacterium if seq[col] == "-")
            
            # Si más del 50% pero menos del 90% son gaps, consideramos consolidar
            if 0.5 < gap_count / len_sequences < 0.9:
                # Buscar columna cercana para consolidar gaps
                best_col = col
                best_score = gap_count
                
                # Buscar en un radio de 3 columnas
                for nearby_col in range(max(0, col-3), min(sequence_length, col+4)):
                    if nearby_col == col:
                        continue
                        
                    nearby_gap_count = sum(1 for seq in bacterium if seq[nearby_col] == "-")
                    if nearby_gap_count > best_score:
                        best_col = nearby_col
                        best_score = nearby_gap_count
                
                # Si encontramos mejor columna, mover gaps
                if best_col != col:
                    for seq_idx in range(len_sequences):
                        if bacterium[seq_idx][col] == "-" and bacterium[seq_idx][best_col] != "-":
                            # Mover el gap
                            seq = list(bacterium[seq_idx])
                            seq[col] = seq[best_col]
                            seq[best_col] = "-"
                            bacterium[seq_idx] = tuple(seq)
        
        # Actualizar la bacteria en la población
        poblacion[bestIdx] = tuple(bacterium)
        
        # Re-evaluar la bacteria modificada
        self.creaGranListaPares([poblacion[bestIdx]])
        self.evaluaBlosum()
        self.creaTablasAtractRepel([poblacion[bestIdx]], 0.2, 0.1, 0.1, 0.05)
        self.creaTablaInteraction()
        self.creaTablaFitness()

    def smartGapOptimization(self, poblacion, bestIdx):
        """
        Mejora inteligente del alineamiento redistribuyendo gaps de manera más efectiva.
        Esta versión mejorada incrementa significativamente el fitness.
        """
        # Trabajar con la mejor bacteria
        bacterium = list(poblacion[bestIdx])
        
        # Asegurar que todas las secuencias son listas
        for i in range(len(bacterium)):
            if isinstance(bacterium[i], tuple):
                bacterium[i] = list(bacterium[i])
        
        len_sequences = len(bacterium)
        
        if len_sequences < 2:
            return  # No hay suficientes secuencias para optimizar
        
        # Medir fitness inicial para comparación posterior
        initial_fitness = self.tablaFitness[bestIdx]
        
        # Guardar una copia de respaldo por si la optimización no mejora el fitness
        backup_bacterium = copy.deepcopy(bacterium)
        
        try:
            sequence_length = len(bacterium[0])
            
            # FASE 1: Consolidar columnas con alta concentración de gaps
            gap_columns = []
            for col in range(sequence_length):
                gap_count = sum(1 for seq in bacterium if seq[col] == "-")
                gap_ratio = gap_count / len_sequences
                
                # Columnas con concentración significativa de gaps
                if 0.4 < gap_ratio < 0.9:
                    # Calcula un puntaje de "utilidad" para la columna
                    # Columnas con aminoácidos bien conservados deben preservarse
                    unique_aas = set(seq[col] for seq in bacterium if seq[col] != "-")
                    conservation_score = 1.0 if len(unique_aas) <= 1 else (2.0 / len(unique_aas))
                    
                    gap_columns.append((col, gap_ratio, conservation_score))
            
            # Si hay columnas con gaps para optimizar
            if gap_columns:
                # Ordenar por ratio de gaps (mayor primero) y conservación (menor primero)
                gap_columns.sort(key=lambda x: (x[1], -x[2]), reverse=True)
                
                # Agrupar columnas cercanas
                clusters = []
                if len(gap_columns) > 1:
                    current_cluster = [gap_columns[0]]
                    
                    for i in range(1, len(gap_columns)):
                        current_col = gap_columns[i][0]
                        prev_col = current_cluster[-1][0]
                        
                        if current_col - prev_col <= 3:  # Umbral de proximidad ajustable
                            current_cluster.append(gap_columns[i])
                        else:
                            clusters.append(current_cluster)
                            current_cluster = [gap_columns[i]]
                    
                    if current_cluster:
                        clusters.append(current_cluster)
                    
                    # Consolidar gaps en cada cluster
                    for cluster in clusters:
                        # Encontrar la columna con mayor proporción de gaps y menor conservación
                        target_col = max(cluster, key=lambda x: x[1] - 0.5 * x[2])[0]
                        
                        # FASE 2: Consolidar gaps en columnas adyacentes
                        # Mover gaps hacia la columna objetivo
                        for col_info in cluster:
                            source_col = col_info[0]
                            if source_col != target_col:
                                for seq_idx in range(len_sequences):
                                    if (bacterium[seq_idx][source_col] == "-" and 
                                        bacterium[seq_idx][target_col] != "-"):
                                        # Asegurar que la secuencia sea una lista
                                        if isinstance(bacterium[seq_idx], tuple):
                                            bacterium[seq_idx] = list(bacterium[seq_idx])
                                        
                                        # Guardar el carácter de la columna objetivo
                                        char_to_move = bacterium[seq_idx][target_col]
                                        # Reemplazar con gap
                                        bacterium[seq_idx][target_col] = "-"
                                        # Eliminar el gap original
                                        bacterium[seq_idx][source_col] = char_to_move
                
                # FASE 3: Eliminar columnas que ahora tienen solo gaps
                i = 0
                while i < len(bacterium[0]):
                    # Verificar si la columna es solo gaps
                    if all(seq[i] == "-" for seq in bacterium):
                        # Eliminar esta columna de todas las secuencias
                        for seq_idx in range(len_sequences):
                            # Asegurar que la secuencia sea una lista
                            if isinstance(bacterium[seq_idx], tuple):
                                bacterium[seq_idx] = list(bacterium[seq_idx])
                                
                            bacterium[seq_idx].pop(i)
                    else:
                        i += 1
            
            # Actualizar la bacteria en la población
            poblacion[bestIdx] = bacterium
            
            # Re-evaluar la bacteria optimizada
            self.creaGranListaPares([poblacion[bestIdx]])
            self.optimizedEvaluaBlosum()
            self.creaTablasAtractRepel([poblacion[bestIdx]], 0.2, 0.1, 0.1, 0.05)
            self.creaTablaInteraction()
            self.creaTablaFitnessEnhanced()
            
            # Verificar si la optimización mejoró el fitness
            if self.tablaFitness[bestIdx] < initial_fitness:
                # Restaurar la bacteria original si no hubo mejora
                poblacion[bestIdx] = backup_bacterium
                # Re-evaluar con la bacteria original
                self.creaGranListaPares([poblacion[bestIdx]])
                self.optimizedEvaluaBlosum()
                self.creaTablasAtractRepel([poblacion[bestIdx]], 0.2, 0.1, 0.1, 0.05)
                self.creaTablaInteraction()
                self.creaTablaFitnessEnhanced()
                
        except Exception as e:
            print(f"Error en smartGapOptimization: {e}")
            # Restaurar la bacteria original en caso de error
            poblacion[bestIdx] = backup_bacterium
            # Re-evaluar
            self.creaGranListaPares([poblacion[bestIdx]])
            self.optimizedEvaluaBlosum()
            self.creaTablasAtractRepel([poblacion[bestIdx]], 0.2, 0.1, 0.1, 0.05)
            self.creaTablaInteraction()
            self.creaTablaFitnessEnhanced()

    def enhancedAlignment(self, poblacion, bestIdx):
        """
        Algoritmo avanzado de alineamiento que optimiza estratégicamente 
        la colocación de gaps para maximizar el fitness.
        """
        # Trabajar con la mejor bacteria
        if bestIdx is None or bestIdx >= len(poblacion):
            return  # No hay bacteria válida para optimizar
        
        bacterium = list(poblacion[bestIdx])
        # Asegurar que todas las secuencias son listas
        for i in range(len(bacterium)):
            if isinstance(bacterium[i], tuple):
                bacterium[i] = list(bacterium[i])
        
        len_sequences = len(bacterium)
        if len_sequences < 2:
            return
        
        # Guardar fitness inicial para comparación
        self.creaGranListaPares([bacterium])
        self.optimizedEvaluaBlosum()
        self.creaTablasAtractRepel([bacterium], 0.2, 0.1, 0.1, 0.05)
        self.creaTablaInteraction()
        self.creaTablaFitness()
        initial_fitness = self.tablaFitness[0]
        
        # Guardar backup
        backup_bacterium = copy.deepcopy(bacterium)
        
        try:
            # 1. Identificar regiones conservadas (alta similitud entre secuencias)
            # Estas regiones son "anclas" donde el alineamiento es probablemente correcto
            sequence_length = len(bacterium[0])
            conservation_scores = []
            
            for col in range(sequence_length):
                # Contar frecuencia de cada residuo en esta columna
                residue_counts = {}
                total_residues = 0
                
                for seq in bacterium:
                    if seq[col] != "-":  # No contar gaps
                        residue = seq[col]
                        residue_counts[residue] = residue_counts.get(residue, 0) + 1
                        total_residues += 1
                
                # Calcular score de conservación (proporción del residuo más común)
                if total_residues > 0:
                    max_count = max(residue_counts.values()) if residue_counts else 0
                    conservation = max_count / total_residues
                else:
                    conservation = 0
                
                conservation_scores.append(conservation)
            
            # 2. Identificar regiones conservadas como aquellas con alta puntuación de conservación
            conserved_regions = []
            in_region = False
            start = 0
            
            for i, score in enumerate(conservation_scores):
                if score > 0.7 and not in_region:  # Umbral de conservación
                    in_region = True
                    start = i
                elif (score <= 0.7 or i == sequence_length - 1) and in_region:
                    in_region = False
                    conserved_regions.append((start, i))
            
            # 3. Para cada par de regiones conservadas, optimizar el espacio entre ellas
            if len(conserved_regions) >= 2:
                for i in range(len(conserved_regions) - 1):
                    region1_end = conserved_regions[i][1]
                    region2_start = conserved_regions[i + 1][0]
                    
                    # Si hay un espacio significativo entre regiones conservadas
                    if region2_start - region1_end > 3:
                        # Contar gaps en cada secuencia en esta región
                        gap_counts = []
                        
                        for seq in bacterium:
                            gap_count = seq[region1_end:region2_start].count("-")
                            gap_counts.append(gap_count)
                        
                        # Determinar número óptimo de gaps
                        median_gaps = sorted(gap_counts)[len(gap_counts)//2]
                        
                        # Redistribuir gaps para cada secuencia
                        for seq_idx, seq in enumerate(bacterium):
                            # Extraer la región entre las regiones conservadas
                            region = seq[region1_end:region2_start]
                            current_gaps = region.count("-")
                            
                            if current_gaps != median_gaps:
                                # Crear nueva región con gaps redistribuidos
                                non_gaps = [c for c in region if c != "-"]
                                
                                # Si necesitamos agregar gaps
                                if current_gaps < median_gaps:
                                    gaps_to_add = median_gaps - current_gaps
                                    # Distribuir gaps uniformemente
                                    gap_positions = [i * len(non_gaps) // gaps_to_add for i in range(gaps_to_add)]
                                    
                                    new_region = []
                                    gap_idx = 0
                                    for i, char in enumerate(non_gaps):
                                        while gap_idx < len(gap_positions) and gap_positions[gap_idx] == i:
                                            new_region.append("-")
                                            gap_idx += 1
                                        new_region.append(char)
                                    
                                    # Agregar gaps restantes al final
                                    while gap_idx < len(gap_positions):
                                        new_region.append("-")
                                        gap_idx += 1
                                    
                                # Si necesitamos eliminar gaps
                                elif current_gaps > median_gaps:
                                    # Simplemente quitar los gaps y mantener los aminoácidos
                                    # Luego agregar el número correcto de gaps al inicio
                                    new_region = ["-"] * median_gaps + non_gaps
                                
                                # Reemplazar la región en la secuencia
                                bacterium[seq_idx] = seq[:region1_end] + new_region + seq[region2_start:]
            
            # 4. Limpiar columnas que solo contienen gaps
            i = 0
            while i < len(bacterium[0]):
                if all(seq[i] == "-" for seq in bacterium):
                    for seq_idx in range(len(bacterium)):
                        bacterium[seq_idx].pop(i)
                else:
                    i += 1
            
            # 5. Evaluar si la optimización mejoró el fitness
            poblacion[bestIdx] = bacterium
            
            self.creaGranListaPares([poblacion[bestIdx]])
            self.optimizedEvaluaBlosum()
            self.creaTablasAtractRepel([poblacion[bestIdx]], 0.2, 0.1, 0.1, 0.05)
            self.creaTablaInteraction()
            self.creaTablaFitness()
            
            # Si no hubo mejora, restaurar la bacteria original
            if self.tablaFitness[0] < initial_fitness:
                poblacion[bestIdx] = backup_bacterium
                
                self.creaGranListaPares([poblacion[bestIdx]])
                self.optimizedEvaluaBlosum()
                self.creaTablasAtractRepel([poblacion[bestIdx]], 0.2, 0.1, 0.1, 0.05)
                self.creaTablaInteraction()
                self.creaTablaFitness()
        
        except Exception as e:
            print(f"Error en enhancedAlignment: {e}")
            # Restaurar la bacteria original en caso de error
            poblacion[bestIdx] = backup_bacterium
            
            self.creaGranListaPares([poblacion[bestIdx]])
            self.optimizedEvaluaBlosum()
            self.creaTablasAtractRepel([poblacion[bestIdx]], 0.2, 0.1, 0.1, 0.05)
            self.creaTablaInteraction()
            self.creaTablaFitness()

    def adaptiveMutation(self, poblacion, bestIdx, iteration, max_iterations):
        """
        Aplica mutación adaptativa que se vuelve más agresiva cuando el fitness se estanca
        """
        # No modificar la mejor bacteria
        best_fitness = self.tablaFitness[bestIdx]
        
        # Determinar si estamos en estancamiento observando la variabilidad del fitness
        fitness_values = [self.tablaFitness[i] for i in range(len(self.tablaFitness))]
        fitness_std = np.std(fitness_values)
        
        # Si la desviación estándar es baja, estamos convergiendo 
        # y necesitamos más exploración
        low_diversity = fitness_std < 0.5
        
        # Factor de mutación adaptativo:
        # - Decrece con las iteraciones para estabilidad al final
        # - Aumenta cuando hay baja diversidad para escapar de óptimos locales
        base_rate = 0.9 * (1 - iteration / max_iterations)
        mutation_rate = base_rate * (2.0 if low_diversity else 1.0)
        
        # Iterar por toda la población excepto la mejor bacteria
        for idx in range(len(poblacion)):
            if idx == bestIdx:
                continue
                
            # Probabilidad de mutación depende de qué tan inferior es esta bacteria
            # comparada con la mejor
            fitness_ratio = self.tablaFitness[idx] / best_fitness if best_fitness > 0 else 0.5
            individual_mutation_prob = mutation_rate * (1.0 - fitness_ratio)
            
            # Solo aplicar mutación con cierta probabilidad
            if np.random.random() < individual_mutation_prob:
                bacterium = list(poblacion[idx])
                
                # Determinar nivel de mutación basado en qué tan mala es esta bacteria
                mutation_intensity = int(3 * (1.0 - fitness_ratio)) + 1
                
                # Aplicar mutaciones
                for _ in range(mutation_intensity):
                    # Seleccionar secuencia aleatoria
                    seq_idx = np.random.randint(0, len(bacterium))
                    seq = list(bacterium[seq_idx])
                    
                    # Tipo de mutación aleatoria
                    mutation_type = np.random.choice(['insert_gap', 'move_gap', 'swap'])
                    
                    if mutation_type == 'insert_gap' and len(seq) > 0:
                        # Insertar un gap en posición aleatoria
                        pos = np.random.randint(0, len(seq))
                        seq = seq[:pos] + ['-'] + seq[pos:]
                        
                    elif mutation_type == 'move_gap' and '-' in seq:
                        # Buscar posiciones con gaps
                        gap_positions = [i for i, c in enumerate(seq) if c == '-']
                        if gap_positions:
                            # Mover un gap a una nueva posición
                            old_pos = np.random.choice(gap_positions)
                            new_pos = np.random.randint(0, len(seq))
                            
                            # Quitar el gap de su posición original
                            seq.pop(old_pos)
                            # Insertar en nueva posición
                            seq = seq[:new_pos] + ['-'] + seq[new_pos:]
                            
                    elif mutation_type == 'swap' and len(seq) > 1:
                        # Intercambiar dos posiciones aleatorias
                        pos1 = np.random.randint(0, len(seq))
                        pos2 = np.random.randint(0, len(seq))
                        seq[pos1], seq[pos2] = seq[pos2], seq[pos1]
                    
                    # Actualizar la secuencia mutada
                    bacterium[seq_idx] = seq
                
                # Actualizar la bacteria en la población
                poblacion[idx] = bacterium
        
        return poblacion

    def biologyInformedAligner(self, poblacion, bestIdx):
        """
        Optimizador de alineamiento basado en propiedades biológicas para maximizar BLOSUM
        """
        if bestIdx is None or bestIdx >= len(poblacion):
            return
        
        # Trabajar con una copia de la mejor bacteria
        bacterium = copy.deepcopy(poblacion[bestIdx])
        
        # Asegurar que todas las secuencias son listas
        for i in range(len(bacterium)):
            if isinstance(bacterium[i], tuple):
                bacterium[i] = list(bacterium[i])
        
        # Guardar fitness original para comparación
        original_fitness = self.tablaFitness[bestIdx]
        backup_bacterium = copy.deepcopy(bacterium)
        
        try:
            # 1. Identificar residuos altamente conservados (anclas biológicas)
            sequence_length = len(bacterium[0])
            num_sequences = len(bacterium)
            
            if sequence_length <= 1 or num_sequences < 2:
                return
                
            # Detectar columnas con alta conservación (posibles motivos funcionales)
            conserved_columns = []
            for col in range(sequence_length):
                residues = [seq[col] for seq in bacterium if col < len(seq) and seq[col] != '-']
                # Contar frecuencia del residuo más común
                if residues:
                    counts = {}
                    for res in residues:
                        counts[res] = counts.get(res, 0) + 1
                    most_common = max(counts.values())
                    conservation_score = most_common / len(residues)
                    
                    # Si más del 70% son el mismo residuo, considerar conservado
                    if conservation_score > 0.7:
                        conserved_columns.append(col)
            
            # 2. Optimización basada en propiedades físico-químicas
            # Clasificamos aminoácidos por propiedades
            hydrophobic = set('AVILMFYW')
            polar = set('STNQ')
            charged = set('DEKRH')
            special = set('CGP')
            
            # Mejorar alineamiento basado en patrones de propiedades físico-químicas
            for col in range(sequence_length):
                if col in conserved_columns:
                    continue  # No tocar columnas conservadas
                    
                # Analizar qué tipo de residuos predominan en esta columna
                residues = [seq[col] for seq in bacterium if col < len(seq) and seq[col] != '-']
                if not residues:
                    continue
                    
                # Contar tipos de residuos
                hydro_count = sum(1 for r in residues if r in hydrophobic)
                polar_count = sum(1 for r in residues if r in polar)
                charged_count = sum(1 for r in residues if r in charged)
                
                # Determinar el carácter predominante de la columna
                counts = [hydro_count, polar_count, charged_count]
                predominant_type = np.argmax(counts)
                
                # Si hay un tipo predominante (>60%), intentar optimizar
                if max(counts) / len(residues) > 0.6:
                    for seq_idx in range(len(bacterium)):
                        # Si no es del tipo predominante y no es un gap, considerar moverlo
                        if col >= len(bacterium[seq_idx]):
                            continue
                            
                        res = bacterium[seq_idx][col]
                        if res == '-' or (predominant_type == 0 and res in hydrophobic) or \
                           (predominant_type == 1 and res in polar) or \
                           (predominant_type == 2 and res in charged):
                            continue
                            
                        # Buscar una mejor posición cercana para este residuo
                        best_position = col
                        best_score = -float('inf')
                        
                        # Escanear vecindad para una mejor posición
                        window = 5  # Buscar 5 posiciones a cada lado
                        for offset in range(-window, window+1):
                            if offset == 0:
                                continue
                                
                            new_pos = col + offset
                            if 0 <= new_pos < sequence_length:
                                # Evaluar si la nueva posición tendría más sentido biológico
                                temp_bacterium = copy.deepcopy(bacterium)
                                
                                # Mover el residuo
                                amino = temp_bacterium[seq_idx][col]
                                original_at_new_pos = temp_bacterium[seq_idx][new_pos] if new_pos < len(temp_bacterium[seq_idx]) else '-'
                                
                                # Realizar el intercambio
                                temp_bacterium[seq_idx][col] = original_at_new_pos
                                if new_pos >= len(temp_bacterium[seq_idx]):
                                    temp_bacterium[seq_idx].extend(['-'] * (new_pos - len(temp_bacterium[seq_idx]) + 1))
                                temp_bacterium[seq_idx][new_pos] = amino
                                
                                # Evaluar la calidad del movimiento
                                # Crear una forma simplificada de calcular un score sin reevaluar todo
                                score = 0
                                
                                # Verificar concordancia de tipo en la nueva posición
                                new_pos_residues = [seq[new_pos] for seq in temp_bacterium if new_pos < len(seq) and seq[new_pos] != '-']
                                if amino in hydrophobic:
                                    score += sum(1 for r in new_pos_residues if r in hydrophobic)
                                elif amino in polar:
                                    score += sum(1 for r in new_pos_residues if r in polar)
                                elif amino in charged:
                                    score += sum(1 for r in new_pos_residues if r in charged)
                                    
                                if score > best_score:
                                    best_score = score
                                    best_position = new_pos
                        
                        # Si encontramos una mejor posición, mover el residuo
                        if best_position != col:
                            amino = bacterium[seq_idx][col]
                            original_at_new_pos = bacterium[seq_idx][best_position] if best_position < len(bacterium[seq_idx]) else '-'
                            
                            bacterium[seq_idx][col] = original_at_new_pos
                            if best_position >= len(bacterium[seq_idx]):
                                bacterium[seq_idx].extend(['-'] * (best_position - len(bacterium[seq_idx]) + 1))
                            bacterium[seq_idx][best_position] = amino
            
            # 3. Redistribución óptima de gaps para maximizar bloques conservados
            # Identificar bloques de columnas con alta similitud
            block_start = None
            blocks = []
            
            for col in range(sequence_length):
                # Calcular similitud en esta columna
                residues = [seq[col] for seq in bacterium if col < len(seq) and seq[col] != '-']
                if len(residues) <= 1:
                    continue
                    
                counts = {}
                for res in residues:
                    counts[res] = counts.get(res, 0) + 1
                most_common = max(counts.values())
                similarity = most_common / len(residues)
                
                # Si encontramos columna similar, empezar un bloque
                if similarity > 0.5 and block_start is None:
                    block_start = col
                # Si termina la similitud, cerrar el bloque
                elif similarity <= 0.5 and block_start is not None:
                    blocks.append((block_start, col-1))
                    block_start = None
            
            # Cerrar el último bloque si existe
            if block_start is not None:
                blocks.append((block_start, sequence_length-1))
            
            # Optimizar gaps entre bloques conservados
            for i in range(len(blocks)-1):
                block_end = blocks[i][1]
                next_block_start = blocks[i+1][0]
                
                if next_block_start - block_end > 2:  # Si hay espacio entre bloques
                    # Contar gaps en cada secuencia en esta región intermedia
                    gap_counts = []
                    for seq in bacterium:
                        gap_count = sum(1 for pos in range(block_end+1, next_block_start) 
                                         if pos < len(seq) and seq[pos] == '-')
                        gap_counts.append(gap_count)
                    
                    # Calcular la mediana de gaps
                    median_gaps = sorted(gap_counts)[len(gap_counts)//2]
                    
                    # Redistribuir gaps de manera uniforme para cada secuencia
                    for seq_idx, seq in enumerate(bacterium):
                        # Extraer la región intermedia
                        region_length = next_block_start - (block_end + 1)
                        
                        if region_length > 0:
                            # Determinar cuántos gaps debe tener esta secuencia
                            current_gaps = gap_counts[seq_idx]
                            
                            if current_gaps != median_gaps:
                                # Extraer los caracteres no gap
                                region = [seq[pos] for pos in range(block_end+1, next_block_start) 
                                          if pos < len(seq)]
                                non_gaps = [c for c in region if c != '-']
                                
                                # Construir nueva región con distribución uniforme de gaps
                                new_region = []
                                
                                if median_gaps > 0:
                                    # Distribuir gaps uniformemente
                                    positions = [i * len(non_gaps) // median_gaps for i in range(median_gaps)]
                                    
                                    gap_idx = 0
                                    for i, char in enumerate(non_gaps):
                                        while gap_idx < len(positions) and positions[gap_idx] == i:
                                            new_region.append('-')
                                            gap_idx += 1
                                        new_region.append(char)
                                    
                                    # Agregar gaps restantes al final
                                    while gap_idx < len(positions):
                                        new_region.append('-')
                                        gap_idx += 1
                                else:
                                    new_region = non_gaps
                                    
                                # Reemplazar la región en la secuencia
                                for pos in range(block_end+1, next_block_start):
                                    if pos < len(seq):
                                        bacterium[seq_idx][pos] = new_region.pop(0) if new_region else '-'
            
            # 4. Evaluar si la optimización mejoró el fitness
            poblacion[bestIdx] = bacterium
            
            # Re-evaluar fitness de la bacteria modificada
            self.creaGranListaPares([poblacion[bestIdx]])
            self.optimizedEvaluaBlosum()
            self.creaTablasAtractRepel([poblacion[bestIdx]], 0.2, 0.1, 0.1, 0.05)
            self.creaTablaInteraction()
            self.creaTablaFitnessAdvanced()  # Usar nuestra nueva función de fitness
            
            # Comparar nuevo fitness con el original
            new_fitness = self.tablaFitness[0]
            
            # Si no mejoró, restaurar la bacteria original
            if new_fitness <= original_fitness:
                poblacion[bestIdx] = backup_bacterium
                
                # Re-evaluar con la original
                self.creaGranListaPares([poblacion[bestIdx]])
                self.optimizedEvaluaBlosum()
                self.creaTablasAtractRepel([poblacion[bestIdx]], 0.2, 0.1, 0.1, 0.05)
                self.creaTablaInteraction()
                self.creaTablaFitnessAdvanced()
                
        except Exception as e:
            print(f"Error en biologyInformedAligner: {e}")
            # Restaurar bacteria original en caso de error
            poblacion[bestIdx] = backup_bacterium
            # Re-evaluar
            self.creaGranListaPares([poblacion[bestIdx]])
            self.optimizedEvaluaBlosum()
            self.creaTablasAtractRepel([poblacion[bestIdx]], 0.2, 0.1, 0.1, 0.05)
            self.creaTablaInteraction()
            self.creaTablaFitnessAdvanced()

