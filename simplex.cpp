// Guilherme Pessoa Silva
// Ciência da Computação - Programação Linear

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cstdio>
#include <cstdlib>
#define MAX 10

typedef struct{
/*1*/	float *a;		// Matriz dos Coef. Restricoes
/*2*/	float *b;		// Vetor das Igualdades das Rest.
/*3*/	float *c;		// Vetor dos Coef. da FC
/*4*/	float *CB;		// Vetor Calculo do custo 
/*5*/	float *x;		// Vetor das incognitas das Rest. (XB)=[x1,x2,x3,x4]
/*6*/	float *Aj;		// Vetor Coluna de a
/*7*/	float *CBt;		// Vetor coluna de c
/* */	float *U;		// Vetor coluna Arm. Binv*Aj
/* */	float *O;		// Vetor coluna 
/* */	float Cj;		// 
/*11*/	double *B;		// Matriz Base escolhida da Mat. a (Valores)
/*12*/	double *Binv;	// Matriz Inversa de B
/*13*/	int *bj;		// Vetor que armaz. os Indices das col. da Mat. B
/*  */	int *cj;		// Vetor que armaz. os Indices das col.(nao esc.) da Mat. B
		
		int i;
		int j;	
}splx;

void preeMat(float *x,int l,int c);
void printMat(splx *si,int op);
void escolheBase(splx *si);
void baseInicial(splx *si);
int  verificaPos(float *x,int l);
void custo(splx *si,int l,int c);
void teta(splx *si,int i);
void MatrizInversa( double *r, const double *m, int ordem);

int main(){
	splx si;
	int pos,pos2,j;
	srand((unsigned int)time(NULL));
	si.i=2;
	printf("\n Matriz A, Digite a linha:");
	scanf("%d",&si.i);
	printf("\n Matriz A, Digite a coluna:");
	scanf("%d",&si.j);//4
	si.a	= (float*)calloc(si.i*si.j,sizeof(float));//Ex.: 2x4
	si.b	= (float*)calloc(si.i,sizeof(float));//tam igual as linhas de a: 2
	si.c	= (float*)calloc(si.j,sizeof(float));//tam igual as colunas de a: 4
	si.CB	= (float*)calloc(si.j,sizeof(float));//tam igual as colunas de a: 4
	si.x	= (float*)calloc(si.j,sizeof(float));//tam igual as colunas de a: 4
	si.Aj  	= (float*)calloc(si.i,sizeof(float));//tam igual as linhas de a: 2
	si.CBt 	= (float*)calloc(si.i,sizeof(float));//tam igual as linhas de a: 2
	si.U   	= (float*)calloc(si.i,sizeof(float));	
	si.O   	= (float*)calloc(si.i,sizeof(float));		
	si.B	= (double*)calloc(si.i*si.i,sizeof(double));//Ordem igual as linhas de a: 2x2
	si.Binv	= (double*)calloc(si.i*si.i,sizeof(double));//Ordem igual as linhas de a: 2x2
	si.bj	= (int*)calloc(si.i,sizeof(int));//tam igual as linhas de a: 2
	si.cj	= (int*)calloc(si.i,sizeof(int));//tam igual as linhas de a: 2
	
	if(!si.a||!si.b||!si.c||!si.CB||!si.x||!si.Aj||!si.CBt||!si.U||!si.B||!si.Binv||!si.bj||!si.cj){
		puts("Sem memoria");
		exit(1);
	}	
	printf("\nPreec. a");
	preeMat(si.a,si.i,si.j);
	printMat(&si,1);	
	printf("\nPreec. b");
	preeMat(si.b,si.i,1);
	printMat(&si,2);	
	printf("\nPreec. c");
	preeMat(si.c,1,si.j);	
	printMat(&si,3);	
	do{
		escolheBase(&si);	
		printMat(&si,11);	
		MatrizInversa(si.Binv,si.B, si.i);
		printMat(&si,12);	
		baseInicial(&si);
		printMat(&si,5);	
		pos=verificaPos(si.x,si.j);
		if(pos){
			printf("\nContinue, Calcular Custo");	
			for(j=0;j<si.i;j++){	
				custo(&si,j,si.cj[j]);	
				pos2=verificaPos(si.CB,si.i);//Tiver neg = 0
				if(!pos2){// !0, 0 = neg
					teta(&si,j);
				}
			}
		}else{
			printf("\nFIM");
		}
	}while(!pos);
	return 0;
}

void teta(splx *si, int j){
	int i,pos;
	float min;
	
	for(i=0;i<si->i;i++){
		pos=verificaPos(si->U,si->i);
		if(pos){
			si->O[i]=si->x[i]/si->U[i];
			printf("\nO[%d]:%.2f /%.2f: %.2f",i,si->x[i],si->U[i],si->O[i]);
		}
	}
	min=si->O[0];
	for(i=0;i<si->i;i++){
		if(min>si->O[i])
			min=si->O[i];
	}
	printf("\nO min: %.2f",min);
}
void custo(splx *si,int l,int c){
	int i,j,a=0;	
	//Pegando um Valor em c
	for(j=0;j<si->j;j++){
		if(j==c){
			si->Cj=si->c[j];//Cj
			printf("\nCj:%.2f",si->Cj);
			break;
		}
	}
	a=0;
	//Pegando os Valores em c
	for(j=0;j<si->j;j++){
		if(j==si->bj[a]){
			si->CBt[a]=si->c[j];//CBt
			printf("\nCBt[%d]:%.2f",a,si->CBt[a]);
			a++;
		}
	}
	a=0;
	//Pegando os Valores em a
	for(j=0;j<si->j;j++){
		if(j==c){
			for(i=0;i<si->i;i++){
				si->Aj[a] = si->a[i*si->j+j];
				printf("\nAj[%d]:%.2f",a,si->Aj[a]);
				a++;
			}
			break;
		}
	}

	//Multip. Binv*Aj=U
	for(i=0;i<si->i;i++){
		a=si->bj[i];
		for(j=0;j<si->i;j++){				
			si->U[i] += si->Binv[i*si->i+j] * si->Aj[j];
			printf("\nU[%d]:%.2f*%.2f=%.2f",i,si->Binv[i*si->i+j],si->Aj[j],si->U[i]);
		}
		//printf("\nx[%d]:%.2f",a,si->XB[a]);
	}

	//Multip. CBt*U
	for(j=0;j<si->i;j++){	
		si->CB[l] += si->CBt[j] * si->U[j];
		//printf("\nx[%d]:%.2f",a,si->XB[a]);
	}
	//CBj = Cj-(CBt*B^-1*Aj)
	si->CB[l] = si->Cj-si->CB[l];
	printf("\nCB[%d]:%.2f",l,si->CB[l]);
}

int verificaPos(float *x, int l){
	int i;
	for(i=0;i < l;i++){
		//printf("\nx[%d]:%.2f",i,si->x[i]);
		if(x[i]<0)
			return 0;//X, pelo menos um negativo (Nao positivo)
	}
	return 1;//X, todo positivo
}

void baseInicial(splx *si){//base Inicial
	int i,j,a=si->bj[0];
	for(i=0;i<si->i;i++){
		a=si->bj[i];
		for(j=0;j<si->i;j++){	
			//printf("\nB[%d]:%.2f*%.2f=%.2f",i,si->Bi[i*si->i+j],si->b[j],si->Bi[i*si->i+j]*si->b[j]);
			si->x[a] += si->Binv[i*si->i+j] * si->b[j];//B*b
		}
		//printf("\nx[%d]:%.2f",a,si->XB[a]);
	}
	
}

void escolheBase(splx *si){//esconhendo a base
	int i,j,p=0,q=0,x=0,y=0;
	
	printf("\nDigite a Base de 0 a %d",si->j-1);
	for(i=0;i<si->i;i++){
		printf("\nDigite a B[%d]: ",i);
		scanf("%d",(si->bj)+i);		
		//si-> bj[i]=rand()%c-1;
		printf("\n%d",si->bj[i]);
	}

	for(i=0;i<si->i;i++){
		for(j=0;j<si->j;j++){		
			if(j==si->bj[q]){
				//printf("\n%d==%d",j,si->bj[q]);
				si->B[x*si->i+y] = si->a[i*si->j+j];
				q++;y++;	
			}else{
				si->cj[p]=j;//Indice nao escolhido na Base
				p++;
			}
		}
		x++;y=0;q=0,p=0;
	}
}
	
void preeMat(float *x,int l,int c){
	int i;	
	
	for(i=0;i<l*c;i++){		
		//printf("\n [%d][%d]:",(i/c),(i%c));
		//scanf("%f",x+i);//4	
		x[i] = rand()%MAX;
	}		
}

void printMat(splx *si,int op){
	int i,l,c;	
	float *p;
	double *q;
	switch(op){
	case 1:
		p=si->a;l=si->i;c=si->j;
		printf("\n\nMat a: ");
	break;	
	case 2:
		p=si->b;l=si->i;c=1;		
		printf("\n\n b: ");
	break;
	case 3:
		p=si->c;l=1;c=si->j;	
		printf("\n\n c: ");	
	break;
	case 4:
		p=si->CB;l=1;c=si->j;	
		printf("\n\n CB: ");	
	break;		
	case 5:
		p=si->x;l=1;c=si->j;	
		printf("\n\n x: ");	
	break;	
	case 11:
		q=si->B;l=si->i;c=si->i;
		printf("\n\nMat B: ");		
	break;
	case 12:
		q=si->Binv;l=si->i;c=si->i;		
		printf("\n\nMat B^-1: ");
	break;			
	}
	
	for(i=0;i< l*c ;i++){		
		if(!(i%c)) printf("\n");
		if(op>=1&&op<=5)
			printf("%.2f ",p[i]);	
		else if (op==11||op==12)
			printf("%.2lf ",q[i]);		
	}		
}

void MatrizInversa( double *r, const double *m, int ordem){
    double *inv = new double[ordem * ordem];
    double *temp = new double[ordem * ordem];

    if (!inv || !temp)
    {
        printf("\n\nERRO: falha na alocação de memoria\n\n");
        exit(1);
    }

    // Copia para 'temp' o conteúdo de 'm'
    for( int i = 0; i < ordem * ordem; i++ )
    {
        temp[ i ] = m[ i ];
    }

    // Transforma 'inv' na matriz identidade
    int j = 0;

    for( int i = 0; i < ordem * ordem; i++ )
    {
        if( i == (j * (ordem + 1)) )
        {
            inv[i] = 1.0;
            j++;
        }
        else inv[i] = 0.0;
    }

    //////  Escalona a parte inferior
    j = 0;
    double pivo;

    for ( int i = 0; i < ordem; i++ )
    {
        if ( temp[ i * ( ordem + 1 ) ] == 0.0 ) /// Verifica se o pivo é nulo
        {
            for ( j = i + 1; j < ordem; j++ )   /// Procura o valor não nulo abaixo do pivo
            {
                if ( temp[ i + j * ordem ] != 0.0 )
                {
                    for ( int k = 0; k < ordem; k++ )
                    {
                        /// Soma a linha do pivo nulo com a linha abaixo
                        temp[ k + i * ordem ] += temp[ k + j * ordem ];
                        inv[ k + i * ordem ] += inv[ k + j * ordem ];
                    }
                    break;
                }
            }

            if ( j == ordem )   /// Se não achar retorna matriz nula
            {
                delete[] temp;
                delete[] inv;
                printf("\n\nERRO: Matriz nao possui inversa\n\n");
                return;
            }
        }

        for ( j = i + 1; j < ordem; j++ )   /// Zera os elementos abaixo do pivô
        {
            if ( temp[ i + j * ordem ] != 0.0 ) /// Ignora se elemento é nulo
            {
                pivo = temp[ i + j * ordem ] / temp[ i * ( ordem + 1 ) ];
                for( int k = 0; k < ordem; k++ )
                {
                    temp[ k + j * ordem ] -= temp[ k + i * ordem ] * pivo;
                    inv[ k + j * ordem ] -= inv[ k + i * ordem ] * pivo;
                }
            }
        }
    }

    //////// Escalona a parte superior
    for( int i = ordem - 1; i >= 0; i-- )
    {
        for( j = i - 1; j >= 0; j-- )   /// Zera os elementos abaixo do pivô
        {
            if( temp[ i + j * ordem ] != 0.0 )  /// Ignora se elemento é nulo
            {
                pivo = temp[ i + j * ordem ] / temp[ i * ( ordem + 1 ) ];
                for( int k = ordem - 1; k >= 0; k-- )
                {
                    temp[ k + j * ordem ] -= temp[ k + i * ordem ] * pivo;
                    inv[ k + j * ordem ] -= inv[ k + i * ordem ] * pivo;
                }
            }
        }
    }

    //////// Transformando os elementos da coluna principal de 'temp' em '1'
    for( int i = 0; i < ordem; i++ )
    {
        pivo = temp[ i * ( ordem + 1 ) ];
        for( j = 0; j < ordem; j++ )
        {
            temp[ j + i * ordem ] /= pivo;
            inv[ j + i * ordem ] /= pivo;
        }
    }

    // Copia os elementos de 'inv' para 'r'
    for ( int i = 0; i < ordem * ordem; i++ )
    {
        r[i] = inv[i];
    }

    delete[] inv;
    delete[] temp;
}

