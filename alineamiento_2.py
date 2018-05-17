
def leerGenoma(ficheroGenoma):
	'''obtiene el genoma leido desde un fichero
        parametro ficheroGenoma: path del fichero'''
	genoma= ''
	with open(ficheroGenoma,'r') as f:
		for linea in f:
			if not linea[0] == '>':
				genoma += linea.rstrip()
	return genoma

def leerReads(ficheroFastq):
	'''Crea una lista de reads que son leidos desde un fichero'''
	secuencias=[]
	with open(ficheroFastq,'r') as fr:
		while True:
			fr.readline() 
			read = fr.readline().rstrip()
			fr.readline()
			fr.readline()
			if len(read)==0:
				break
			secuencias.append(read)
	return secuencias
	
def crearIndice(genoma,size_k_mero):
	'''crea un diccionario con los k-meros del genoma
		clave: es el k-mero
		valor: una lista de posiciones del k-mero'''
	indice={}
	for i in range(len(genoma)-size_k_mero+1):
		k_mero = genoma[i:i+size_k_mero]
		#Si la clave ya existe en el diccionario, agrega un elemento al valor
		#caso contrario, agrega la clave y valor al diccionario
		if( k_mero in indice): 
			indice[k_mero].append(i)
		else:
			indice[k_mero] = [i]
	return indice


def alineamientoAprox(p,t):
	'''Realiza un alineamiento aproximado de P con respecto a TabError
	Retorna la cantidad minima de edicion y la posicion de inicio del alineamiento'''
#Inicializacion de la matriz-----------------------------------------------------------
	D=[]
	for i in range(len(p)+1):
		D.append([0]*(len(t)+1))
		
	for i in range(len(p)+1):
		D[i][0] = i
	
#Calculo de la distancia edicion --------------------------------------------------------	
	for i in range(1,len(p)+1):
		for j in range(1,len(t)+1):
			distHor = D[i][j-1]+1
			distVer=D[i-1][j]+1
			if p[i-1] == t[j-1]:
				distDiag = D[i-1][j-1]
			else:
				distDiag = D[i-1][j-1]+1
			D[i][j] = min(distHor, distVer, distDiag)
#Se elige el menor elemento de la ultima fila----------------------------------------------
	distEdicion = min(D[-1])
	c = D[-1].index(distEdicion)#la posicion columna del ultimo elemento
	f = len(D)-1# la ultima fila
	
	recorrido = []#Almacenara el recorrido 

	while(f > 0 and c > 0):
		diag = D[f-1][c-1]
		arriba = D[f-1][c]
		izquierda = D[f][c-1]
		elemento = D[f][c]
		
		recorrido.append((f,c))

		#Si hay sustitucion o match
		if(diag <= arriba and diag <= izquierda):

			f=f-1
			c=c-1
		#Si hay insercion
		if(arriba < diag and arriba < izquierda):
			
			f=f-1
		#Si hay borrado
		if(izquierda < diag and izquierda < arriba):
			
			c=c-1
			#la cantidad de edicion y donde empieza el alineamiento menos el caracter vacio
	return distEdicion,recorrido[-1][1]-1
	
#Busqueda en el indice------------------------------------------------------------------------
def obtenerHits(particion):
	'''Obtiene una lista con las ocurrencias de la particion'''
	hits = []
	if(particion in indice):
		hits = indice[particion]
	
	return hits

#Consulta las posiciones al alinear read con k ediciones
def consultarEdicion(read,genoma,size_k_mero,k):
	'''Retorna una lista de ocurrencias de read en genoma
	La lista contiene una tupla de (edicion, posicion dentro de genoma)'''
	
	#Crear particiones del read, donde cada particion tiene k_mero de tamaño 
	# para realizar la consulta al indice
	read_particion = []
	fin = len(read)
	inicio = 0
	while(fin >= inicio + size_k_mero):
		particion = read[inicio:inicio+size_k_mero]
		
		read_particion.append((particion,inicio)) #(particion, posicion dentro de read)
		
		inicio+=size_k_mero
		
	#Buscar donde ocurre los hit para cada particion del read dentro del genoma
	ocurrencias = []
	vista = set()
	
	for part, pinicio in read_particion:
		
		for hit in obtenerHits(part):
			
			#lado izq. de Genoma para incluir en la matriz
			izq = max(0,hit - pinicio - k)
			#lado der. de Genoma para incluir en la matriz
			der = min( len(genoma), hit - pinicio + len(read) + k )
			distEdicion,posicion = alineamientoAprox(read,genoma[izq:der])
			posicion+=izq
			
			if(distEdicion <= k and (distEdicion,posicion) not in vista):
				ocurrencias.append((distEdicion,posicion))
				vista.add((distEdicion,posicion))
				
			
	return ocurrencias



#Programa---------------------------------------------------------------------------------------


# path del fichero que contiene el genoma y los reads

ficheroGenoma= "phix.fa"
ficheroFastq = "reads_phix_1.fastq"
k=2 #nro de ediciones permitidas

print("Leyendo fichero del genoma...")
genoma = leerGenoma(ficheroGenoma)

print("Leyendo fichero del reads...")
lista_reads = leerReads(ficheroFastq)
 
#Calculamos la longitud del k-mero para el indice
size_k_mero = len(lista_reads[0])//(k+1)

print("Creando Indice con ",size_k_mero,"-mero...")
indice = crearIndice(genoma,size_k_mero)

print("Realizando alineamiento de reads...")

#Diccionario donde la clave esla posicion donde comienza la alineacion y el valor es una tupla 
#con el nucleotido y la distancia edicion

read_alineados = {}
for read in lista_reads:
	ocurrencias = consultarEdicion(read,genoma,size_k_mero,k)
	for distEdicion,posicion in ocurrencias:
		for i in range(len(read)):
			if( i+posicion in read_alineados): 
				nucleotido, edicion = read_alineados[i+posicion]
				
				if(distEdicion <= edicion):
				
					read_alineados[i+posicion]=(read[i], distEdicion)
			else:
				read_alineados[i+posicion] = (read[i],distEdicion)

genoma_ind = ""
for i in sorted(read_alineados.keys()):
	
	nucleotido,edicion = read_alineados[i]
	genoma_ind = genoma_ind+nucleotido
	#print(i, "->", nucleotido, " ",edicion)
			
			
print("Escribiendo fichero...")		
with open("genoma22.fasta","w") as f:
	print(">gi",file=f)
	inicio = 0
	fin= len(genoma_ind)
	while(inicio < fin):
		print(genoma_ind[inicio:inicio+70],file=f)
		inicio=inicio+70

#conocer las diferencias
c=0
for i in range(len(genoma)):
	if(genoma_ind[i:i+1] != genoma[i:i+1]):
		c=c+1
		
print(c)




	


