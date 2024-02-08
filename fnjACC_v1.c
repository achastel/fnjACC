#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENACC
#include <openacc.h>
#endif
#define MIN(a,b) ((a) < (b) ? (a) : (b))
//#include <iostream.h>
struct p_node {
	long int d;
	int i,j;
};

struct otuName {
   char *name;
   int ordem;
};

//**************************************
// Driver Code
//**************************************
int main(int argc, char *argv[])
{
    FILE* ptr = fopen(argv[1], "r");	//arquivo contendo matriz de distancia triangular
    double f = strtod(argv[2],NULL);   //f=p/q percentual de pares a serem selecionados de uma vez
    clock_t t_start, t_end;
    double dt;

    if (ptr == NULL) {
       printf("no such file.\n");
       return 0;
    }

    int  size;					// No. OTU's
    fscanf(ptr,"%d", &size);
    printf("SIZE=%d\n",size);
    float s = 0.0;
    s = (size/2);
    s = (size+1)*s;
    int triangularSize = s-size;	//tamanho da matriz triangular
    int i,j,k;
    long int **m;				// MATRIZ DE DISTANCIAS
    long int *Q[size];			// MATRIZ Q
    long int D[size],Dseq[size];   // SOMA DAS DISTANCIAS DE UMA OTU(linha)
    int indices[size];
    int t=(size*f)+1;
    int limite=t;
    struct p_node L[t];	          // LISTA DE MENORES
    struct otuName otus[size];	// armazena nome/nós/ramos da árvore
    m=(long int**)malloc(size*sizeof(long int));
    int t_vetor;
    int menor,maior;
    long int mij, Di, Dj;
    double t1 = 0.00;
    t1=(size)/2;
    t1=t1*size;
    t_vetor = t1;
    struct p_node vetor[t_vetor];

    printf ("t1=%f t_vetor=%d size=%d \n",t1,t_vetor, size);
    printf(" ALOCAÇÃO LINHAS DAS MATRIZES E INICIALIZAÇÃO\n");

    for(int i=0;i<size;i++){
       D[i]=0;Dseq[i]=0;
	     m[i]= (long int*) malloc((i+1)*sizeof(long int));
	     Q[i]= (long int*) malloc(i*sizeof(long int));
	     for(int j=0;j<i;j++) m[i][j]=-1;
	     otus[i].name = (char *) malloc(8*sizeof(char));
	     sprintf(otus[i].name,"%d",i);
       otus[i].ordem=i;
	     indices[i]=0;
    }
    printf("Alocando Q\n");
    for(int i=1;i<size;i++)
	     Q[i]= (long int *) malloc(i*sizeof(long int));

    printf("Leirura da Matriz Triangular:\n");
    i=0;
    while (i<size){									// LENDO DISTANCIAS DO ARQUIVO...
       j=0;
       while(j<=i){
	       fscanf(ptr, "%ld ", &m[i][j]);					// ...inserindo na matriz
	       j++;
	     }
	     i++;
    }
    fclose(ptr);

    limite=(size*f)+1;
    struct p_node pares[limite];						//HEAP COM PARES DE VIZINHOS

    t_start = clock();
    //printf("Soma das Distacias na liha i:\n SEQUENCIAL \n");
    //for(i=0;i<size;i++){		// calculando soma das distancias das OTUS - D[i]
    //   for(j=0;  j<i;   j++) Dseq[i]=Dseq[i]+m[i][j];
    //   for(j=i+1;j<size;j++) Dseq[i]=Dseq[i]+m[j][i];
    //}

    //for(i=0;i<size;i++)
    //   printf("D[%d] = %ld \n",i,Dseq[i]);

    printf("Calculando soma das distancias das OTUS Paralelo- D[i]\n");
    #pragma acc data copy(D,m)
    #pragma acc parallel
    {
	  #pragma acc loop
       for(i=0;i<size;i++)
          #pragma acc loop
          for(j=0;  j<size;   j++){
             if (j<i ) D[i]=D[i]+m[i][j];
             if (j>i)  D[i]=D[i]+m[j][i];
       }
    }
    //printf("Soma das Distacias na liha i:\n");
    //for(i=0;i<size;i++) printf("D[%d] = %ld == %ld \n",i,D[i],Dseq[i]);

    printf("N. OTUs=%d f=%f\n",size,f);

    //**************************************** LAÇO PRINCIPAL

    while(size>=2){
       limite = (size*f);
       if (limite==0) limite=1;
       k=0;
       if (size>2){
          //********************************** CALCULO DA MATRIZ Q
          #pragma acc data copy(Q,m,D,vetor,indices,k)
          #pragma acc parallel private(i,j)
          {
             #pragma acc loop
             for(i=0;i<size;i++){
                indices[i]=0;
                #pragma acc loop
                for(j=0;j<i;j++){
                   vetor[k].i=i;
                   vetor[k].j=j;
                   Q[i][j]=((long int)(size - 2) * m[i][j]) - (D[i] + D[j]);
                   vetor[k].d = Q[i][j];
                   k++;
                }
             }
          }

          //printf("\n MATRIZ Q \n");
          //for(i=0;i<size;i++) {
          //   for(j=0;j<i;j++) printf("[%d,%d] %ld ",i,j,Q[i][j]);
          //      printf("\n");
          //}
          #pragma acc update self(vetor)
          //k=0;
          //printf("\n VETOR \n");
          //for(i=0;i<size;i++) {
          //   for(j=0;j<i;j++){
          //      printf("[%d] %d,%d,%ld ",k,vetor[k].i,vetor[k].j,vetor[k].d);
          //      k++;
          //   }
          //   printf("\n");
          //}

          t_vetor=k;
          printf ("k=%d t_vetor=%d \n",k,t_vetor);
          //********************************** ORDENANDO O VETOR
          struct p_node aux;
          #pragma acc data copy(pares, vetor, t_vetor)
          #pragma acc parallel
          {
             #pragma acc loop
             for(i=0;i<t_vetor;i++)
                #pragma acc loop
                for(j=0;j<t_vetor-1;j++)
                   if(vetor[j].d > vetor[j+1].d ){
                      aux=vetor[j];
                      vetor[j]=vetor[j+1];
                      vetor[j+1]=aux;
                   }
          }

          k=0;
          #pragma acc update self(vetor)
          while(k<limite){
             if(k==0){
                //printf("%d, ",k);
                pares[k]=vetor[k];
                indices[vetor[k].i]=1;
                indices[vetor[k].j]=1; }
             else {
                int j=0;
                while((j<(t_vetor)) && !( indices[vetor[j].i]==0 && indices[vetor[j].j]==0))
                   j++;
                //printf("%d, ",j);
                pares[k]=vetor[j];
                indices[vetor[j].i]=1;
                indices[vetor[j].j]=1;
             }
             k++;
          }//while(k<limite)
          //printf("\n PARES\n");
          //for(i=0;i<limite;i++) {
          //    printf("[%d] (%d,%d) %ld \n",i, pares[i].i, pares[i].j, pares[i].d);
          //}
          //printf("SIZE=%d\n",size);

          printf("Juntando %d Pares de OTU's\n",limite);
          for(int c=0; c<limite; c++) {
             if (size >2){
                double d_i_novo,d_j_novo=0.0;
	              double t1, t2;
                if (vetor[c].i > vetor[c].j) {
                   mij=m[pares[c].i][pares[c].j]; Di=D[pares[c].j]; Dj=D[pares[c].i];
                   menor=pares[c].j;
                   maior=pares[c].i; }
                else {
                   mij = m[pares[c].j][pares[c].i]; Di = D[pares[c].i]; Dj = D[pares[c].j];
                   menor=pares[c].i;
                   maior=pares[c].j; }
                //printf("menor=%d=>%d maior=%d=>%d \n",menor, otus[menor].ordem, maior, otus[maior].ordem);

                int y;
                char aux[(strlen(otus[menor].name)+strlen(otus[maior].name)+4)];
                if(otus[menor].ordem < otus[maior].ordem)
                    sprintf(aux,"(%s,%s)",otus[menor].name,otus[maior].name);
                else
                   sprintf(aux,"(%s,%s)",otus[maior].name,otus[menor].name);

                if(otus[menor].ordem >= otus[maior].ordem){
                   int ax=otus[menor].ordem;
                   otus[menor].ordem = otus[maior].ordem;
                   otus[maior].ordem = ax; }
                //printf("menor=%d=>%d maior=%d=>%d \n",menor, otus[menor].ordem, maior, otus[maior].ordem);


                y = 8 - (strlen(aux)%8) + 1;
                if (otus[menor].name) free(otus[menor].name);

                if(!(otus[menor].name=(char*) malloc((strlen(aux)+y)*sizeof(char))))
                   printf("Não alocou*************\n");
                //printf("y=%d %s \n",y,aux);
	              strcpy(otus[menor].name, aux);
                //printf("menor=%d maior=%d %s \n",menor, maior, aux);

                t1= (size - 2) * 2;
                t2= (Di - Dj);
                t2 = t2/t1;
                //printf("t1 = %f  t2=%f \n",t1,t2);
                d_i_novo = (0.5*mij);
                d_i_novo = d_i_novo + t2 ;
                d_j_novo = mij - d_i_novo;

                //printf("calculando distancia do novo aos demais -- calculando a linha-coluna\n");
                //#pragma omp for nowait
                for(int k=0;k<menor;k++){
		         D[menor] = D[menor] - (m[menor][k]);
		         D[k] = D[k] - (m[menor][k]+m[maior][k]);
		         m[menor][k] = (int)((m[menor][k]+m[maior][k])-mij)*0.5;
		         D[k] = D[k] + m[menor][k];
		         D[menor] = D[menor] + m[menor][k];
		      }
		      //#pragma omp for nowait
		      for(int k=menor+1;k<maior;k++){
                   D[menor] = D[menor]-m[k][menor];
		         D[k] = D[k] - (m[k][menor]+m[maior][k]);
		         m[k][menor]= (int)((m[k][menor]+m[maior][k])-mij)*0.5;
		         D[k] = D[k] + m[k][menor];
		         D[menor]=D[menor]+m[k][menor];
		      }
                //#pragma omp for nowait
		      for(int k=maior;k<size;k++){
		         D[menor] = D[menor]-m[k][menor];
		         D[k] = D[k] - (m[k][menor]+m[k][maior]);
		         m[k][menor] = (int)((m[k][menor]+m[k][maior])-mij)*0.5;
		         D[k] = D[k] + m[k][menor];
		         D[menor]=D[menor]+m[k][menor];
		      }

		         //printf("Eliminando linha-coluna com cópia .\n");
	              //#pragma omp for nowait
	              for(int k=0;k<size;k++)
	                 if (k>maior)
	                    m[k][maior]=m[size-1][k];
	                 else
	                    m[maior][k]=m[size-1][k];

                //#pragma omp single
	              //{
	               if((otus[maior].name) && (maior != (size-1))) free(otus[maior].name);
	               otus[maior].name = otus[size-1].name;
	               otus[maior].ordem=otus[size-1].ordem;
	            //}
              //for(int x1=0;x1<size;x1++)
              //   printf("[%d]=%d %s\n",x1,otus[x1].ordem,otus[x1].name);

              //printf("Corrigindo coordenadas.\n");
              //for(int x1=0;x1<size-1;x1++)
              //   printf("[%d]=(%d,%d)\n",x1,pares[x1].i,pares[x1].j);
	            //#pragma omp for nowait
	            for(int e=c+1; e<limite; e++){
	               if(pares[e].i==size-1)
	                  if (pares[e].j > maior) {
		                pares[e].i=pares[e].j;
		                pares[e].j=maior;
	                  }
	                  else pares[e].i=maior;
	               if (pares[e].j==size-1)
	                  pares[e].j=maior;
	            }
              //printf("CORRIGIDO\n");
              //for(int x1=0;x1<size-1;x1++)
              //   printf("[%d]=(%d,%d)\n",x1,pares[x1].i,pares[x1].j);

	            D[maior]=0;
	            //#pragma omp for nowait
                 for(int k=0;k<size-1;k++) D[maior]=D[maior]+m[size-1][k];
             }//if (size >2)
             size--;
          }//for(int c=0; c<limite; c++)
       }//if(size>2)
       else size--;
   }//while(size>=2)
   char aux[(strlen(otus[0].name)+strlen(otus[1].name)+4)];

   if(otus[0].ordem  < otus[1].ordem)
     sprintf(aux,"(%s,%s)",otus[0].name,otus[1].name);
   else
     sprintf(aux,"(%s,%s)",otus[1].name,otus[0].name);

   t_end = clock();
   dt = ((double)(t_end-t_start))/ (double)CLOCKS_PER_SEC;
   printf("%s\n",aux);
   printf("\nTime =%18.4f sec\n",dt);

   return 0;
}
