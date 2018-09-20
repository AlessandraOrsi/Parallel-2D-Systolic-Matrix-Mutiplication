/*
 =====================================================================================================================================
 Name        : MatrixMultiplication4S.c
 Author      : Alessandra Orsi, Claudia Pipino, Daniele Vitale
 Version     : 2.0
 Description : Bisogna effettuare il calcolo del prodotto C = A x B con A[m][h], B[h][n] e C[m][n].
 			   Bisogna suddividere la matrici A e B in p blocchi di righe e p blocchi di colonne, ottenendo
 			   p^2 blocchi di A e B.
 
 Strategy	 : TORO
 =====================================================================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h> 
#include <unistd.h>

/*	FUNZIONI UTILIZZATE	*/

void initMatrix(float*, int, int, int);

void printMatrix(float*, int, int);

float prod_scal(float*, float*, int);

void mulMatrix(float*, float*, float*, int, int);


/*	MAIN	*/

int main(int argc, char** argv){

/*	MPI VARIABLES	*/
	int menum;
	int nproc;
	int tag;
	int row, col;
	int dim, *ndim,reorder, *period;
	int coordinate[2];
	MPI_Comm comm_grid;
	MPI_Status status;

/*	INPUT(A, B), LOCAL(local_A, local_B) AND RESULT(C, local_C) MATRIX	*/
	float *A;
	float *B;
	float *local_A;
	float *local_B;
	float *C;
	float *local_C;
	float *tmpA;
	float *tmpB;

/*	MATRIX VARIABLES	*/
	int n, m;
	int i, j, k, w;
	int local_n, local_m;
	int x1A, y1A, x2A, y2A;
	int x1B, y1B, x2B, y2B;

/*	FLAG PER LA STAMPA E LA RANDOMIZZAZIONE DELLA MATRICE*/
	int flag, rand;

/*	VARIABILI CATTURA DEL TEMPO	*/
	double start, finish;


/*	COORDS	*/
	int north[2], south[2], east[2], west[2];
/*	MPI INITIALIZE	*/
	MPI_Init(&argc, &argv);
/*	NUMBER OF PROCESSORS	*/
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
/*	GIVE PERSONAL ID TO ALL PROCESSORS	*/
	MPI_Comm_rank(MPI_COMM_WORLD, &menum);

	int p = sqrt(nproc);
/*	Controllo numero di processi quadrato	*/
	if((p * p) != nproc){
		if(menum == 0)
			printf("ERRORE! IL NUMERO DI PROCESSORI DEVE ESSERE TALE DA POTER CREARE UNA GRIGLIA QUADRATA\nInsight:\n\tInserisci ad esempio 4 o 9 processori\n\n");
		MPI_Finalize();
		exit(0);
	}
/*	Controllo valori da riga di comando		*/
	if(argc<5){
		if(menum == 0)
			printf("Verrà eseguito il programma con dimensioni standard (matrici %d*%d) \n\n", p, p);
		n = p;
		m = p;
		flag = 1;
		rand = 1;
	}else{
		n = atoi(argv[1]);	//A's rows, B's columns
		m = atoi(argv[2]);	//B's rows, A's columns
		flag = atoi(argv[3]);	//Flag per la stampa
		rand = atoi(argv[4]);	//Flag per la randomizzazione delle matrici
	/*	Controllo */
		if((m%p !=0)||(n%p !=0)){
			if(menum == 0)
				printf("Inserire un numero di righe e colonne\n\tdivisibile per la radice del numero di processi\n");		
			MPI_Finalize();

			exit(0);
		}
	}

	row = p;
	col = p;
/*	Creazione griglia toroidale	*/
/*	Spedizione di row da parte di 0 a tutti i processo	*/
	MPI_Bcast(&row,1,MPI_INT,0,MPI_COMM_WORLD);
	
/*	Numero di diimensioni della griglia	*/
	dim = 2;
	
/*	Vettore contenente le lunghezze di ciascuna dimensione	*/
	ndim = (int*)calloc(dim,sizeof(int));
	ndim[0] = row;
	ndim[1] = col;
	
/*	Vettore contenente la periodicità delle dimensioni	*/
	period = (int*)calloc(dim,sizeof(int));
	period[0] =1; period[1] = 1;
	reorder = 0;
	
/*	Definizione della griglia bidimensionale	*/
	MPI_Cart_create(MPI_COMM_WORLD,dim,ndim,period,reorder,&comm_grid);
	
/*	Assegnazione delle coordinate a ciascun processo nella griglia bidimensionale	*/
	MPI_Cart_coords(comm_grid,menum,dim,coordinate);
	
	local_n = n/p;
	local_m = m/p;

	local_A = (float*) malloc(local_n * local_m * sizeof(float));
	local_B = (float*) malloc(local_m * local_n * sizeof(float));

	local_C = (float*) calloc(local_n * local_n, sizeof(float));

	int sendCoords[2];
	int recvCoords[2];

	if(menum == 0){

/*	Inizio cattura tempo	*/
		start = MPI_Wtime();
		A = (float*) malloc(n * m * sizeof(float));
		B = (float*) malloc(m * n * sizeof(float));
		
		tmpA = (float*) malloc(n/p * m/p * sizeof(float));
		tmpB = (float*) malloc(m/p * n/p * sizeof(float));

		initMatrix(A, n, m, 0);
		initMatrix(B, m, n, 0);
		if(flag){
			printf("***\tSTAMPA MATRICE A\t***\n");
			printMatrix(A, n, m);
			printf("***\tSTAMPA MATRICE B\t***\n");
			printMatrix(B, m, n);
		}
		
		for(i = 0; i < local_n; i++)
			for(j = 0; j < local_m; j++)
				local_A[i * local_n + j] = A[i * n + j];


		for(i = 0; i < local_m; i++)
			for(j = 0; j < local_n; j++)
				local_B[i * local_m + j] = B[i * m + j];
			
		for(i = 1; i < nproc; i++){
			MPI_Cart_coords(comm_grid, i, 2, sendCoords);
			x1A = sendCoords[0] * (n/p);
			x2A = x1A + (n/p);
			y1A = (sendCoords[0] + sendCoords[1]) % p * (n/p);
			y2A = y1A + (n/p);
			x1B = y1A;
			x2B = y2A;
			y1B = sendCoords[1] * (n/p);
			y2B = y1B + (n/p);
			for(j = x1A; j < x2A; j++){
				for(w = y1A; w < y2A; w++){
					tmpA[(j % local_n) * local_n + w % local_m] = A[j * n + w];
				}
			}
			
			MPI_Send(tmpA, local_n*local_m, MPI_FLOAT, i, 0, MPI_COMM_WORLD);

			for(j = x1B; j < x2B; j++){
				for(w = y1B; w < y2B; w++){
					tmpB[(j % local_m) * local_m + w % local_n] = B[j * m + w];
				}
			}
	
			MPI_Send(tmpB, local_m*local_n, MPI_FLOAT, i, 1, MPI_COMM_WORLD);

		}
		
	}else{
		MPI_Recv(local_A, local_n*local_m, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(local_B, local_m*local_n, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, &status);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if(menum == 0){
		free(A);
		free(B);
		free(tmpA);
		free(tmpB);
	}

	north[0] = (coordinate[0]-1)%p;
	north[1] = coordinate[1];
	south[0] = (coordinate[0]+1)%p;
	south[1] = coordinate[1];
	east[0] = coordinate[0];
	east[1] = (coordinate[1]+1)%p;
	west[0] = coordinate[0];
	west[1] = (coordinate[1]-1)%p;

	mulMatrix(local_C, local_A, local_B, local_n, local_m);

	int rank;
	int rank1;
	for(i = 0; i < p-1; i++){

		MPI_Cart_rank(comm_grid, west, &rank);
		MPI_Cart_rank(comm_grid, east, &rank1);
	
		MPI_Sendrecv_replace(local_A, local_n*local_m,MPI_FLOAT,rank,0,rank1,0,MPI_COMM_WORLD,&status);

		MPI_Cart_rank(comm_grid, north, &rank);
		MPI_Cart_rank(comm_grid, south, &rank1);
		
		MPI_Sendrecv_replace(local_B, local_m*local_n,MPI_FLOAT,rank,0,rank1,0,MPI_COMM_WORLD,&status);

		mulMatrix(local_C, local_A, local_B, local_n, local_m);
		MPI_Barrier(MPI_COMM_WORLD);
	}


	if(menum == 0){
		float *tmpC = (float*) malloc(n/p * n/p * sizeof(float));
		C = (float*) calloc(n * n, sizeof(float));

		for(i = 0; i < local_n; i++)
			for(j = 0; j < local_n; j++)
				C[i*n+j] = local_C[i*local_n+j];

		for(i = 0; i < p; i++)
			for(j = 0; j < p; j++){
				if(!((i == 0)&&(j ==0))){
					recvCoords[0] = i;
					recvCoords[1] = j;

					MPI_Cart_rank(comm_grid, recvCoords, &rank);
					MPI_Recv(tmpC, local_n*local_n, MPI_FLOAT, rank, 0, MPI_COMM_WORLD, &status);

					for(k = recvCoords[0]*(n/p); k < n; k++)
						for(w = recvCoords[1]*(n/p); w < n; w++)
							C[k*n + w] = tmpC[(k%local_n)*local_n + w%local_n];
				}

			}

				
	}else{
		MPI_Send(local_C, local_n*local_n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
	}
	if(menum == 0)
		finish = MPI_Wtime();	
	MPI_Finalize();

	if(flag){
		if(menum == 0){
			printf("***\tSTAMPA MATRICE C\t***\n");
			fflush(stdout);
			printMatrix(C, n, n);
		}
	}


	if(menum == 0){
		printf("\n\n***\nTempo trascorso: %f\n***\n\n", finish-start);
		free(C);
	}	
	free(local_A);
	free(local_B);
	free(local_C);	
	return(0);
}


/*	USED FUNCTION	*/

void initMatrix(float* M, int n, int m, int how){

	int i,j;
	if(how == 0){
		for(i = 0; i < n; i++)
			for(j = 0; j < m; j++)
				M[i * n + j] = (float)rand()/((float)RAND_MAX+(float)1); 
	}else{
		for(i = 0; i < n; i++)
			for(j = 0; j < m; j++)
				M[i * n + j] = i * n + j; 
	}
}


void printMatrix(float* M, int n, int m){

	int i,j;
	for(i = 0; i < n; i++){
		for(j = 0;j < m; j++)
			printf("%f\t", M[i * m + j]);
		printf("\n");
	}
	printf("\n\n");

}

void mulMatrix(float* C, float* A, float* B, int n, int m){
	
	int i,j;
	float* Bt = (float*)malloc(n * m * sizeof(float));

	for(i=0; i<n; i++)
		for(j=0; j<m; j++)
			Bt[i*n+j] = B[j*m+i];

	for(i = 0; i < n; i++)
		for(j = 0;j < n; j++)
			C[i * m + j] += prod_scal(&A[i*n], &Bt[j*n], n);

	free(Bt);
}


float prod_scal(float* A, float* B, int n){
	
	float toRet = 0;
	int i;

	for(i=0; i<n; i++){
		toRet += A[i]* B[i];
	}
	return toRet;
			
}






















