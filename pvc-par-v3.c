#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"

#define SEED 132

typedef struct caminho_
{
    int* caminho;
    int peso;
} caminho;


// Variáveis globais
//ATENÇÃO, AO EXECUTAR O MPI, AS VÁRIAVEIS GLOBAIS NÃO SÃO TRANSFERIDAS, É NECESSÁRIO INICIALIZÁ-LAS EXPLICITAMENTE
int N;
int** pesoglobal;

int* alocavetorint(int size);
int** alocamatrizint(int linhas, int colunas);
void liberavetorint(int** vec);
void liberamatrizint(int** vec, int size);
void gerapesos(int** pesos, int size);
caminho* encontraCaminho(caminho* ateaqui,int quantosvisitados, int noatual, int peso_encontrado);
caminho* alocacaminho();
void liberacaminho(caminho* cam);
void swap (int *i, int *j);
int *generate_permutation(int N, int n_visitadosfixos, int qnt_caminhos);
void permute(int esquerda, int* arr, int direita, int **cam_permu, int *counter, int size);
caminho *melhor_caminho(int tam_visistadosfixos, int n_visitadosfixos, int *visitadosfixos);
void aumenta_caminhos(int *qnt_caminhos, int *n_visitadosfixos, int nodes, int n_vertices);

// Retonar um ponteiro para um vetor de inteiros do tamanho solicitado
int* alocavetorint(int size){

    int* result;

    result=(int*)malloc(sizeof(int)*size);
    if(result==NULL) {
        printf("Erro de alocação de memória no vetor de inteiros\n");
        exit(1);
    }

    return result;
}

// Retonar um ponteiro para uma matriz de inteiros do tamanho solicitado
int** alocamatrizint(int linhas, int colunas){

    int i;
    int** result;

    result=(int**)malloc(sizeof(int*)*linhas);
    if(result==NULL) {
        printf("Erro de alocação de memória no vetor de inteiros\n");
        exit(1);
    }

    for(i=0;i<linhas;i++){
       result[i]=alocavetorint(colunas); 
    }

    return result;
}

// Libera a memória alocada por um vetor de inteiros ao passálo como referência &vet
void liberavetorint(int** vec){

    if(vec!=NULL){
    free(*vec);
    *vec=NULL;
    }
    
    return;
}

// Aloca uma estrutura caminho, que contém a ordem e o peso total de um caminho no grafo
caminho* alocacaminho(){

    caminho* cam;
    cam=(caminho*)malloc(sizeof(caminho));
    cam->caminho=alocavetorint(N);
	for(int i = 0; i < N; i++)
		cam->caminho[i] = 0;
    return cam;
}

// Libera a memória alocada por uma matriz de inteiros ao informar o número de linhas da tal
void liberamatrizint(int** vec, int size){

    int i;

    for(i=0;i<size;i++){
        liberavetorint(&vec[i]);
    }

    free(vec);
    
    return;
}

// Libera a memória de um caminho solicitado
void liberacaminho(caminho* cam){

    liberavetorint(&(cam->caminho));
    free(cam);

    return;
}

// Gera de maneira pseudo-aleatória os pesos para a matriz de conectividade do grafo
void gerapesos(int** pesos,int size){

    int i,j;

    srand(SEED);

    for(i=0;i<size;i++){
        for(j=0;j<size;j++){
            if(i==j) pesos[i][j]=0;
            else pesos[i][j]=rand()%101;// 0 a 100
            if(pesos[i][j]==0 && i!=j) pesos[i][j]=(int)INFINITY; // O 0 significa não conexo (peso infinito)
        }
    }

    return;
}

// Função recursiva que encontra o melhor caminho através da força bruta (verifica os N! caminhos restantes)
// Dentro dela ela atualiza o caminho enviado como parâmetro, o desalocando e substituindo por um atualizado, isso é transparente ao usuário
// Basta notar que é necessário enviar e receber um caminho
caminho* encontraCaminho(caminho* ateaqui,int quantosvisitados, int noatual, int peso_encontrado){

    int i,j;
    int menorpeso,menorindice;
    int* aux;
    caminho* cam;

    menorpeso= (int)INFINITY; // Peso do pior caso

    cam=alocacaminho();

    for(i=0;i<N;i++){ //Copia o caminho e o peso até aqui
        cam->caminho[i]=ateaqui->caminho[i];
    }
    cam->peso=ateaqui->peso;

    aux=alocavetorint(N);

    if(quantosvisitados<N){ // Ainda preciso avançar
        for(i=0;i<N;i++){ //Para todos os nós
			
            if(ateaqui->caminho[i]==0 && cam->peso < peso_encontrado){ // Se eu não visitei, eu visito
                cam->caminho[i]=quantosvisitados+1; // Modifico o caminho
                cam->peso+=pesoglobal[noatual][i];
                cam=encontraCaminho(cam,quantosvisitados+1,i, peso_encontrado); // Avanço
                if(cam->peso<menorpeso) { // Comparo resultados
                    menorpeso= cam->peso;
                    for(j=0;j<N;j++){
                        aux[j]=cam->caminho[j];
                    }
                }
                // Volto
                for(j=0;j<N;j++){ 
                    cam->caminho[j]=ateaqui->caminho[j];
                }
                cam->peso=ateaqui->peso;
                
                
            }
        } 
        // Pego o Melhor
        for(j=0;j<N;j++){
            cam->caminho[j]=aux[j];           
        }
        cam->peso=menorpeso;
    }
    else{ // Completei um caminho, só falta voltar p/ o 0
			cam->peso+=pesoglobal[noatual][0];
    }

    liberacaminho(ateaqui);
    liberavetorint(&aux);

    return cam;

}

// Função auxiliar que imprime um vetor de inteiros dado seu tamanho
void print_vec(int *vec, int n){
	for(int i = 0; i < n; i++){
		printf("%d ", vec[i]);
	}
	printf("\n");
}

// Função auxiliar que imprime uma matriz de inteiros dado seu tamanho
void print_mat(int **mat, int rows, int col){
	for(int i = 0; i < rows; i++){
		print_vec(mat[i], col);
	}
	printf("\n");
}


// dada uma sequencia de [1, ..., n], gera (n_visitadosfixos) permutacoes diferentes 
// de tamanho = qnt_caminhos
// onde cada máquina irá considerar essas permutacoes como caminhos iniciais
int *generate_permutation(int N, int n_visitadosfixos, int qnt_caminhos){
		int size_permutations = qnt_caminhos * n_visitadosfixos;

		int *vec_aux = alocavetorint(N);
		int **cam_permutations = alocamatrizint(qnt_caminhos, n_visitadosfixos);
		int *permutations = alocavetorint(qnt_caminhos * n_visitadosfixos);

		for(int i = 0; i < N; i++){
			vec_aux[i] = i+1;
		}
		for(int i = 0; i < size_permutations; i++){
			permutations[i] = i+1;
		}
		int counter = 0;
		if(n_visitadosfixos > 1){
			permute(0, vec_aux, N-1, cam_permutations, &counter, n_visitadosfixos-1);
			for(int i = 0; i < qnt_caminhos; i++){
				memcpy(permutations + i*n_visitadosfixos, cam_permutations[i], sizeof(int) * n_visitadosfixos);
			}
		}
		liberavetorint(&vec_aux);
		liberamatrizint(cam_permutations, qnt_caminhos);

		return permutations;
}
// Função auxiliar chamada pela funcao 'generate_permutation' 
// que gera todas permutacoes e seleciona as desejaveis
void permute(int esquerda, int* arr, int direita, int **cam_permu, int *counter, int size) {
	if (esquerda == direita-1) {
		int antecessor;
		if(*counter == 0)
			 antecessor = *counter;
	 	else
			antecessor = *counter - 1;
		if(cam_permu[antecessor][size] != arr[size]){
				memcpy(cam_permu[*counter], arr, sizeof(arr));
				(*counter)++;
		 }
		 return;
   }
   	for (int i = esquerda; i < direita; i++) {
		swap (arr + esquerda, arr + i);
	   	permute(esquerda+1, arr, direita, cam_permu, counter, size);
       	swap (arr + i, arr + esquerda);
   }
   return;
}

// Função auxiliar Faz troca de dois valores inteiros
void swap (int *i, int *j) {
    int temp = *i;
    *i = *j;
    *j = temp;
}

// Toda maquina chama essa funcao 1 vez.
// marcando *visitadosfixos comos caminhos iniciais
// e fazendo a busca apartir disto, retorna o melhor caminho encontrado 

caminho *melhor_caminho(int tam_visistadosfixos, int n_visitadosfixos, int *visitadosfixos){
	caminho *resultado = alocacaminho();
	int *caminho_fixo = alocavetorint(n_visitadosfixos);
	int visitados_p_node = tam_visistadosfixos / n_visitadosfixos;

	resultado->peso = (int) INFINITY;
	for(int i = 0; i < visitados_p_node; i++){
		memcpy(caminho_fixo, visitadosfixos+(i*n_visitadosfixos), sizeof(int)*n_visitadosfixos);
		caminho *cam_atual = alocacaminho();
		cam_atual->peso = 0;
		cam_atual->caminho[0] = 1;
		int src = 0;

		for(int j = 0; j < n_visitadosfixos; j++){
			cam_atual->caminho[caminho_fixo[j]] = j+2;
			cam_atual->peso += pesoglobal[src][caminho_fixo[j]];
			src = caminho_fixo[j];
		}
		cam_atual = encontraCaminho(cam_atual,n_visitadosfixos+1,src, resultado->peso);

		if(cam_atual->peso < resultado->peso){
			memcpy(resultado->caminho, cam_atual->caminho, sizeof(int) * N);
			resultado->peso = cam_atual->peso;
		}
		liberacaminho(cam_atual);
	}
	return resultado;
}

// Se existem mais maquinas do que qnt caminhos (iniciado por n_vertices-1), 
// expande os caminhos possiveis em (vertices * vertices-1)
// Ate a qnt de caminhos for superior ao numero de maquinas,
// afim de garantir uma distribuição melhor de trabalho.
			
void aumenta_caminhos(int *qnt_caminhos, int *n_visitadosfixos, int nodes, int n_vertices){
	if(nodes <= *qnt_caminhos || n_vertices==*n_visitadosfixos+1){
		return;
	}
	*qnt_caminhos *= (n_vertices-*n_visitadosfixos-1);		// 1 it = N
	(*n_visitadosfixos)++;
	aumenta_caminhos(qnt_caminhos, n_visitadosfixos, nodes, n_vertices);
}

int main (int argc, char* argv[]){

	int Mprocessadores, myrank, n_visitadosfixos, qnt_caminhos,resto, ret, tam_visistadosfixos;
	int root = 0;
	// Inicializa o processo MPI, daqui para baixo todas as M máquinas executarão as instruções
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &Mprocessadores); //Recebe o  número máximo de processadores usados no comando MPIrun
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank); // Recebe o número da máquina atual (lembrando que todas executam, cada uma terá um rank diferente)

    int i,j,sum;
    int** pesos;
    caminho* resultado;
	int* caminhofinal;
	int *vetor_sizes, *vetor_dsp;
	int *distribuicao;
	int *acumulada;

	//Necessários para juntar os caminhos no final
	int* melhores_caminhos;
	int *caminhos_iniciais;


	//Faz o parse do argumento N (numero de nós no grafo)
	if (argc != 2){
			printf("Erro na execução do código\n");
			return EXIT_FAILURE;
	}

	else N = atoi(argv[1]);

	if(N<=1){
			printf("Erro na execução do código\n");
			return EXIT_FAILURE;
	}

	//Aloca espaço em memória para as estruturas necessárias
	pesos = alocamatrizint(N,N);
	gerapesos(pesos,N);

	if(myrank == root){ // Root é a máquina 0, ou seja, nesse if apenas uma máquina entra, e ela será responsável pelo controle da aplicação
		//Apenas ela gera os pesos pseudo-aleatórios (Caso contrário as máquinas poderiam gerar pesos diferentes)
			gerapesos(pesos,N);
			qnt_caminhos = N-1;
			n_visitadosfixos = 1;
			
			// Expande a quantidade de caminhos ate que ela seja maior que Mprocessadores
			aumenta_caminhos(&qnt_caminhos, &n_visitadosfixos, Mprocessadores, N);

			//Vector_sizes define quantos dados a máquina i vai receber
			// Vector_dsp representa o deslocamento que deve ser feito para ler o valor
			//Por exemplo se eu precisasse mandar 1 processo para a máquina 0, 3 para a máquina 1 e 5 para a máquina 2:
			//Vector sizes seria: 1 3 5 e vector_dsp seria 0 1 2 
			vetor_sizes = alocavetorint(Mprocessadores);
			vetor_dsp = alocavetorint(Mprocessadores);

			tam_visistadosfixos = (qnt_caminhos / Mprocessadores) * n_visitadosfixos;
			resto = qnt_caminhos  % Mprocessadores;

			for(i = 0; i < Mprocessadores; i++){
				vetor_sizes[i] = tam_visistadosfixos;
				if(i < resto){
					vetor_sizes[i] += n_visitadosfixos;
				}
			}
			int sum = 0;
			for(i = 0; i < Mprocessadores; i++){
				vetor_dsp[i] = sum;
				sum += vetor_sizes[i];
			}
			// gera permutacoes de tamanho n_visitadosfixos para marcar como visitado
			// com objetivo de fazer uma busca mais curta
			caminhos_iniciais = generate_permutation(N, n_visitadosfixos, qnt_caminhos);
	}

	//Aqui, a máquina 0 vai espalhar os pesos para todas as outras, que precisam desse dado para calcular os caminhos
	for(i = 0; i < N; i++) MPI_Bcast(pesos[i], N, MPI_INT, root, MPI_COMM_WORLD);	
	pesoglobal=pesos; // Cada máquina terá um ponteiro diferente como "pesoglobal", mas todos apontam para matrizes de mesmo valor

	// faz broadcast das informacoes que foram calculadas na maquina raiz.
	MPI_Bcast(&tam_visistadosfixos, 1, MPI_INT, root, MPI_COMM_WORLD);
	MPI_Bcast(&resto, 1, MPI_INT, root, MPI_COMM_WORLD);
	MPI_Bcast(&n_visitadosfixos, 1, MPI_INT, root, MPI_COMM_WORLD);

	// Inicialmente os caminhos iniciais sao divididos igualmente entre as maquinas
	// Depois sera distribuido o resto entre os primeiros processos
	if(myrank < resto){
		tam_visistadosfixos += n_visitadosfixos;
	}

	int *visitadosfixos = alocavetorint(tam_visistadosfixos);

	//caminhos_iniciais existia apenas na maquina 0, o scatterv vai espalhar os caminhos
	//cada maquina recebera uma parcela em visitadosfixos
	MPI_Scatterv(caminhos_iniciais,vetor_sizes,vetor_dsp,MPI_INT,visitadosfixos,tam_visistadosfixos,MPI_INT,root,MPI_COMM_WORLD);
	// Apos cada maquina saber o seu ponto de partida, aplica o algoritmo comparando o melhor
	resultado = alocacaminho();


	resultado = melhor_caminho(tam_visistadosfixos, n_visitadosfixos, visitadosfixos);

	int global_resultado_rank[2];
	int result_rank[2];

	result_rank[0] = resultado->peso;
	result_rank[1] = myrank;

	// Reducao MINLOC
	// Ira retornar em global_resultado_rank[0] = melhor resultado
	// e em global_resultado_rank[1] = rank onde isso ocorreu
	MPI_Allreduce(result_rank, global_resultado_rank, 1, MPI_2INT, MPI_MINLOC, MPI_COMM_WORLD);

	if(myrank == global_resultado_rank[1]){
		//Imprime o resultado final
		printf("Matriz de pesos:\n");
		print_mat(pesos,N,N);

		printf("Caminho encontrado: \n");
		for(i=1;i<=N;i++){
			for(j=0;j<N;j++){
				if(resultado->caminho[j]==i){
					printf("%d ",j);
					break;
				} 
			}
		}
		printf("0\n");
		printf("Custo do caminho encontrado = %d\n", resultado->peso);
	}
	if(myrank == root){ // Junta os caminhos
		liberavetorint(&caminhos_iniciais);
		liberavetorint(&vetor_sizes);
		liberavetorint(&vetor_dsp);
	}
	// Libera as estruturas alocadas
	liberacaminho(resultado);
	liberamatrizint(pesos,N);
	//Verifica retorno
	//MPI_Barrier(MPI_COMM_WORLD);

	ret = MPI_Finalize();
	if (ret == MPI_SUCCESS){
		printf("MPI_Finalize success! \n");
        return EXIT_SUCCESS;
    }
	else{
        printf("MPI_Finalize failure! \n");
        return EXIT_FAILURE;
    }
}
