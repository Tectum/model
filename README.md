# funciones generales

def getCoords(i,width):
    ''' Transforma un indice en una coordenada'''
    x_coord = int(i / width)
    y_coord = int(i % width)
    return x_coord, y_coord

def toIndex((x,y),width):
    ''' Transforma una coordenada en un indice'''
    return x * width + y

def calc_dist(a,b,width):
    ''' Calcula la distancia entre dos indices que cooresponden a coordenadas'''
    a = getCoords(a,width)
    b = getCoords(b,width)
    return sqrt((a[0]-b[0])**2+(a[1]-b[1])**2)

def calc_dist_2(x,y): 
    ''' Calcula la distancia entre un indice y el centro del campo receptivo
    de una neurona del Ipc''' 
    y = (y[0]*3 + 7, y[1]*3 + 7)  # centro de los RF del ipc en el estimulo
    d0 = abs(x[0] - y[0])
    d1 = abs(x[1] - y[1])
    return sqrt(d0**2+d1**2)    
        
def stimulus(posx,posy,side,rate):
    ''' Genera una matriz que da cuenta de la posicion del estimulo en el espacio bidimensional.
    Esto es lo que recibe el PoissonGroup'''
    stimulus = np.zeros((size_VF, size_VF))  # siempre un cuadrado de campo visual
    stimulus[posx : posx + side*conv, posy : posy + side*conv] = 1.0
    stimulus = stimulus.ravel() * rate
    return stimulus

def retraso(x,y): 
    ''' Calcula el retraso de la respuesta de una neurona del Ipc frente a estimulos en distintas
    posiciones de su campo receptivo. Ambos argumentos se deben expresar en coordenadas'''
    distance = calc_dist_2(x, y)
    return  (3.2 * distance + 32) * ms  

# funciones que describen las conecciones entre grupos neuronales

def stim_to_ipc(from_,to_):
    ''' Pesos sinapticos en la conexion entre estimulo y neuronas del Ipc'''
    m = 15; sep = 3        
    for i in range(N_estimulo):
        post = []
        a,b = getCoords(from_,sqrt(N_estimulo))
        for k in pivotes_ipc:
            if (k[0] <= a <= k[0] + m) and (k[1] <= b <= k[1] + m):
                post.append(toIndex((k[0]/sep,k[1]/sep),sqrt(N_ipc)))   
        for j in range(N_ipc):
            if to_ in post:
                return 15 * mV
            else: 
                return 0 * mV          

def  ipc_to_bb(from_,to_):
    ''' Pesos sinapticos en la conexion entre neuronas del Ipc y bottlebrushes'''
    size_paintb = 2; m = 15        
    for i in range(N_ipc):
        x,y = pivotes_ipc[from_]        
        for j in range(N_bottlebrushes):
            a,b = getCoords(to_,sqrt(N_bottlebrushes))
            if (x+(m/2) <= a <= x+(m/2)+size_paintb) and (y+(m/2) <= b <= y+(m/2)+size_paintb):     # (c == a) and c == b):
                return 1 * mV
            else:
                return 0 * mV          

ipc_to_ipc = lambda i, j:((exp(-(calc_dist(i,j,19)) * 1.)*10)-9) * mV 

def stim_to_imc(from_,to_):
    ''' Pesos sinapticos en la conexion entre el estimulo y neuronas del Imc'''
    ancho = 15; largo = 70; sep_imc = 3    
    for i in range(N_estimulo):
        post = []
        a,b = getCoords(from_,sqrt(N_estimulo))
        for k in pivotes_imc:
            if (k[0] <= a <= k[0] + largo) and (k[1] <= b <= k[1] + ancho):
                post.append(toIndex((k[0]/sep_imc,k[1]/sep_imc),sqrt(N_imc)))   
        for j in range(N_imc):
            if to_ in post:
                return 15 * mV
            else: 
                return 0 * mV

def imc_to_ipc(from_,to_):
    ''' Pesos sinapticos en la conexion entre las neuronas del Imc y las neuronas del Ipc'''
    for i in range(N_imc):
        for j in range(N_ipc):
            (a,b) = getCoords(to_,sqrt(N_ipc))
            if abs(b - from_) <= 1:
                return 0 * mV
            else:
                return -3.0 * mV 

def imc_to_imc(from_,to_):
    ''' Pesos sinapticos en la conexion entre neuronas del Imc. Inhibicion lateral'''
    for i in range(N_imc):
        for j in range(N_imc):
            if from_ != to_:
                return -15 * mV
            else: 
                return 0 * mV

def ipc_to_imc(from_,to_):
    ''' Pesos sinapticos en la conexion entre las neuronas del Ipc y las neuronas del Imc. 
    Proyeccion homotopica excitatoria'''
    for i in range(N_ipc):
        for j in range(N_imc):
            (a,b) = getCoords(from_,sqrt(N_ipc))
            if abs(to_ - b) <= 1:
                return 15 * mV
            else:
                return 0 * mV                

def pivotes_ipc():
    ''' Genera una lista de coordenadas que marcan las posiciones de los campos receptivos
    de las neuronas del Ipc independientemente del tamaño de estos. Ademas se define el 
    tamaño''' 
    pivotes_ipc = []
    sep = 3  # separacion entre los campos receptivos, da una medida de superposicion
    for i in range(int(sqrt(N_ipc))):   # numero de pivotes por lado
        w = sep*i
        for j in range(int(sqrt(N_ipc))):
            z = sep*j
            pivotes_ipc.append((w,z))
    return pivotes_ipc
            
def pivotes_imc():    
    ''' Genera una lista de coordenadas que marcan las posiciones de los campos receptivos
    de las neuronas del Imc independientemente del tamaño de estos. Ademas se define el 
    tamaño (largo y ancho)'''
    pivotes_imc = []
    sep_imc = 3 # separacion entre los campos receptivos, da una medida de superposicion
    for i in range(1): # numero de pivotes por lado
        w = sep_imc*i
        for j in range(19):
            z = sep_imc*j
            pivotes_imc.append((w,z))
    return pivotes_imc              

# funciones relacionadas con la estimulacion

def change_position2(t):
    ''' Cambio en la posicion del estimulo durante la simulacion para dos estimulos'''
    posx1 = x0_stim1 + np.sin(angle1) * (vel_converted_1) * t 
    posy1 = y0_stim1 + np.cos(angle1) * (vel_converted_1) * t 
    stim1 = stimulus(posx1, posy1, side1, rate1)
    if t > t_initial_stim2:
        posx2 = x0_stim2 + np.sin(angle2) * (vel_converted_2) * (t-t_initial_stim2) 
        posy2 = y0_stim2 + np.cos(angle2) * (vel_converted_2) * (t-t_initial_stim2) 
        stim2 = stimulus(posx2, posy2, side2, rate2)
        return stim1 + stim2
    else:
        return stim1 
      
def freq_RGC(vel,tam):
    ''' Tuneo de las ganglionares al tamaño y a la velocidad de los estimulos'''
    tun_vel = 404.7* (1-(np.exp(-vel/22.6)))
    tun_tam = 1 * np.exp(-(tam-2)**2/2*1**2)
    return tun_vel * tun_tam * Hz
    
''' Dos estimulos. Se definen los parametros del campo visual y de los estimulos.
Se crea el PoissonGroup'''
## caracteriticas del estimulo 1
N_estimulo = 70*70    
init1 = (20,10)
angle1 = 0
x0_stim1 = init1[0]
y0_stim1 = init1[1]
velocity_1 = 8 * Hz # dada en grados por segundo
conv = 2 # factor de conversion de grados a bottlebrushes     
vel_converted_1 = velocity_1 * conv  # cantidad de bb que recorre por segundo
side1 = 2
rate1 = freq_RGC(velocity_1,side1)
## caracteristicas del estimulo 2 
init2 = (40,10)
angle2 = 0
x0_stim2 = init2[0]
y0_stim2 = init2[1]
velocity_2 = 12 * Hz   # velocidad maxima para caer dentro del locus es 20º/s
vel_converted_2 = velocity_2 * conv     
side2 = 2
rate2 = freq_RGC(velocity_2,side2)
t_initial_stim2 = 0 * msecond # tiempo en el cual aparece el segundo estimulo, dado en segundos
lado = 70  
size_VF = lado
stim = PoissonGroup(N_estimulo, change_position2)

# creacion de los grupos neuronales

''' Define el numero de neuronas por grupo neuronal, las ecuaciones que rigen su 
dinamica de voltaje y se crean los NeuronGroup. Son tres grupos neuronales: los
bottlebrushes en la capa 5 del tectum, las neuronas del Ipc y del Imc'''
N_bb = 70    
N_bottlebrushes = N_bb**2
N_ipc = 19*19
N_imc = 19*1

C = 200 * pF
gL = 10 * nS
taum = C / gL
EL = -58 * mV
VT = -50 * mV
DeltaT = 2 * mV
Vcut = 0 *mV

tauw = 3500*ms   #3500ms inicialmente
a = 2*C/(3500*ms)
b = 0.01*nA  # 100 queda bien 
Vr = -60* mV

eqs_i = """
dvm/dt=(-gL*(vm-EL)+gL*DeltaT*exp((vm-VT)/DeltaT)+I-w)/C : volt
dw/dt=(a*(vm-EL)-w)/tauw : amp
I : amp
"""
tau = 18 * ms  # 20 ms con adaptacion en parvo
tau_e = 2 * ms # AMPA synapse
eqs = '''
dv/dt=(I-v)/tau : volt
dI/dt=-I/tau_e : volt
'''
L5 = NeuronGroup(N_bottlebrushes, model=eqs, reset=0 * mV, threshold=10 * mV)
Ipc = NeuronGroup(N_ipc, model=eqs_i, threshold=Vcut, reset="vm=Vr;w+=b", freeze=True)
Imc = NeuronGroup(N_imc, model=eqs, reset=0 * mV, threshold=10 * mV,refractory=0*ms)
Ipc.vm = EL

# creacion de las conecciones entre grupos neuronales 

''' Sinapsis entre el estimulo y los bottlebrushes. Se basa en la dinamica probabilistica
descrita por Luksch'''
# parametros de la sinapsis 
p_max = 0.87
t_0 = 2025. * ms 
S = Synapses(stim, L5, model = '''w : 1
				  p : 1
				  prob : 1
				  lastupdate : second''',
		     pre = '''prob = p_max * (1 - exp(-((t - lastupdate)/t_0)))
                       p = (lastupdate==0)*p_max + (lastupdate!=0)*prob
                       v += w * p''')
S[:,:]='i==j' # conexion homotopica
S.w[:,:] = 11 * mV # peso sinaptico no suficiente para generar espigas en las postsinapticas

''' Se crean las listas necesarias para la construccion de las matrices de pesos sinapticos
para los campos receptivos del Ipc y Imc'''
pivotes_ipc = pivotes_ipc()
pivotes_imc = pivotes_imc()          

''' Se crea la conexion entre el estimulo y las neuronas del Ipc'''
RF_ipc = Connection(stim, Ipc, sparseness=0.11, delay=True, structure='dense', weight=stim_to_ipc)

''' Se crea la conexion entre las neuronas del Ipc y los bottlebrushes'''
paintbrush = Connection(Ipc, L5, sparseness=1., weight=ipc_to_bb)

''' Se crea la conexion entre el estimulo y las neuronas del Imc'''
RF_imc = Connection(stim, Imc, sparseness=0.8, delay=True, structure='dense', weight=stim_to_imc)
        
''' Se crea la conexion entre las neuronas del Imc y las neuronas del Ipc'''
heterotopic = Connection(Imc, Ipc, sparseness=1., weight=imc_to_ipc) 

''' Se crea la conexion lateral inhibitoria entre neuronas del Imc'''
inh_imc = Connection(Imc,Imc, sparseness=1., weight=imc_to_imc)

''' Se crea la conexion homotopica excitatoria del Ipc al Imc'''
ipc_imc = Connection(Ipc,Imc,sparseness=1.,weight=ipc_to_imc) 

# output de la simulacion 
''' Registrar los tiempos de espigas'''
M = SpikeMonitor(L5)
N = SpikeMonitor(Ipc)
O = SpikeMonitor(Imc)
L = SpikeMonitor(stim)

''' Registrar variables de estado'''
rec_bot = StateMonitor(L5,'v',record=True, timestep=100)
rec_ipc = StateMonitor(Ipc,'vm',record=True, timestep=100)
rec_stim = StateMonitor(stim,'rate', record= True, timestep=100)

''' Correr la simulacion una cantidad de segundos'''
run(2000*msecond)
