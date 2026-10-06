import sys
from Bio import SeqIO
from matplotlib import pyplot as plt

def buscar_origem(arquivo_fasta, dnaA_box, origem_real=None):
    """
    Busca a origem de replicação de um genoma bacteriano.

    Parâmetros:
    arquivo_fasta (str): Caminho para o arquivo fasta.
    dnaA_box (str): Sequencia de nucleotídeos
    origem_real (int, optional): Posição da origem de replicação real.

    """

    #Passos:
    #1. Ler o arquivo fasta
    #2. Calcular a assimetria GC acumulada
    #3. Encontrar o mínimo e máximo da assimetria GC acumulada (candidatos a oriC e ter respectivamente)
    #4. Refinar a busca da origem de replicação utilizando a abundância de A e T por janela e a abundância de dnaA box por janela

    genoma = SeqIO.read(arquivo_fasta, "fasta")
    sequencia = str(genoma.seq)
    
    assimetria_GC_acumulada = calcular_assimetria_GC_acumulada(sequencia)
    
    minimo = min(assimetria_GC_acumulada) #valor minimo, indicador de oriC
    maximo = max(assimetria_GC_acumulada)

    #Caso haja mais de um mínimo ou máximo, retorna uma lista com as posições
    if assimetria_GC_acumulada.count(minimo) == 1:
        minimo_assimetria_GC = [assimetria_GC_acumulada.index(minimo)]
    elif assimetria_GC_acumulada.count(minimo) > 1:
        minimo_assimetria_GC = [i for i, x in enumerate(assimetria_GC_acumulada) if x == minimo]

    #REFINAMENTO
    margem_para_refinamento = 10000
    tamanho_janela_de_leitura = 1000
    passo_de_leitura = 250

    #Posição inicial e final da região candidata a oriC
    regiao_candidata = [minimo_assimetria_GC[0] - margem_para_refinamento,
                         minimo_assimetria_GC[-1] + margem_para_refinamento]

    origem_estimada = refinamento(sequencia, regiao_candidata, dnaA_box, tamanho_janela_de_leitura, passo_de_leitura)

    #plot assimetria_GC_acumulada
    plt.plot(range(len(assimetria_GC_acumulada)), assimetria_GC_acumulada)
    plt.xlabel("Posição")
    plt.ticklabel_format(style='plain')
    plt.ylabel("Assimetria GC Acumulada")
    plt.title(genoma.id)
    #plt.show()
    plt.savefig(genoma.id + "_desvioGC.png")

    #plotar outro grafico apenas com a região candidata a oriC
    plt.clf()
    plt.plot(range(regiao_candidata[0], regiao_candidata[1]), assimetria_GC_acumulada[regiao_candidata[0]:regiao_candidata[1]])
    if origem_real is not None:
        plt.axvline(x=origem_real, color='g', linestyle='-', label='Origem de replicação real')
    plt.axvline(x=origem_estimada, color='r', linestyle='--', label='Origem de replicação estimada')
    plt.xlabel("Posição")
    plt.ticklabel_format(style='plain')
    plt.ylabel("Assimetria GC Acumulada")
    plt.title(genoma.id)
    plt.savefig(genoma.id + "_erro.png")
    return None 

def calcular_assimetria_GC_acumulada(sequencia):
    """
    Calcula a assimetria GC acumulada de um sequencia bacteriano.

    Parâmetros:
    sequencia (str): Seqência de nucleotídeos.

    Retorna:
    assimetria_GC_acumulada [list]: Lista com os valores de assimetria GC acumulada.
    """
 
    assimetria_GC_acumulada = []
    
    acumulo = 0
    for base in sequencia:
        if base == 'G':
            acumulo = (assimetria_GC_acumulada[-1] + 1 if assimetria_GC_acumulada else 1)
        elif base == 'C':
            acumulo = (assimetria_GC_acumulada[-1] - 1 if assimetria_GC_acumulada else -1)
        else:
            acumulo = (assimetria_GC_acumulada[-1] if assimetria_GC_acumulada else 0)

        assimetria_GC_acumulada.append(acumulo)
        
    return assimetria_GC_acumulada

def calcular_abundancia_AT_por_janela(sequencia, tamanho_janela_de_leitura = 1000, passo_de_leitura = 250):
    """
    Calcula a abundância de A e T em uma sequência de nucleotídeos por janela.

    Parâmetros:
    sequencia (str): Seqência de nucleotídeos.
    tamanho_janela_de_leitura (int): Tamanho da janela.
    passo_de_leitura (int): Passo para mover a janela.

    Retorna:
    abundancia_AT [list]: Lista com os valores de abundância de A e T por janela.
    """
    abundancia_AT = []
    
    for i in range(0, len(sequencia) - tamanho_janela_de_leitura + 1, passo_de_leitura):
        janela = sequencia[i:i+tamanho_janela_de_leitura]

        quantidade_A = janela.count('A')
        quantidade_T = janela.count('T')
        abundancia_AT.append((quantidade_A + quantidade_T)/tamanho_janela_de_leitura)

    return abundancia_AT

def calcular_abundancia_motif_por_janela(sequencia, motif, tamanho_janela_de_leitura = 1000, passo_de_leitura = 250):
    """
    Calcula a abundância de um motif em uma sequência de nucleotídeos por janela.

    Parâmetros:
    sequencia (str): Seqência de nucleotídeos.
    motif (str): Motif a ser buscado.
    tamanho_janela_de_leitura (int): Tamanho da janela.
    passo_de_leitura (int): Passo para mover a janela.

    Retorna:
    abundancia_motif [list]: Lista com os valores de abundância do motif por janela.
    """
    abundancia_motif = []
    
    for i in range(0, len(sequencia) - tamanho_janela_de_leitura + 1, passo_de_leitura):
        janela = sequencia[i:i+tamanho_janela_de_leitura]
        abundancia_motif.append(janela.count(motif))
        
    return abundancia_motif

def refinamento(sequencia, regiao_candidata, dnaA_box, tamanho_janela_de_leitura = 1000, passo_de_leitura = 250):
    """
    Refinamento da busca da origem de replicação.

    Parâmetros:
    sequencia (str): Sequência completa.
    regiao_candidata (tuple): [0]posição inicial [1]posição final da região.
    dnaA_box (str): Sequência do motif DnaA box.
    tamanho_janela_de_leitura (int): Tamanho da janela.
    passo_de_leitura (int): Passo para mover a janela.

    Retorna:
    origem_estimada (int): posição da origem de replicação.
    """
    
    regiao_candidata_sequencia = sequencia[regiao_candidata[0]:regiao_candidata[1]]

    #Refinar a região candidata a oriC
    abundancia_AT = calcular_abundancia_AT_por_janela(regiao_candidata_sequencia, tamanho_janela_de_leitura, passo_de_leitura)
    abundancia_dnaA_box = calcular_abundancia_motif_por_janela(regiao_candidata_sequencia, dnaA_box, tamanho_janela_de_leitura, passo_de_leitura)

    pontuacao = []
    for i in range(len(abundancia_AT)):
        pontuacao.append(2*abundancia_AT[i] + 5*abundancia_dnaA_box[i])

    origem_estimada = regiao_candidata[0] + (pontuacao.index(max(pontuacao)) * passo_de_leitura)
    print(f"Origem de replicação estimada na posição: {origem_estimada}, "+
           f"AT%: {abundancia_AT[pontuacao.index(max(pontuacao))]}, "+
           f"contagem dnaA box: {abundancia_dnaA_box[pontuacao.index(max(pontuacao))]}")
    return origem_estimada


buscar_origem("encontros/encontro_29_09_26/data/GCF_000005845.2_ASM584v2_genomic.fna", "TTATCCACA", 3923883)
