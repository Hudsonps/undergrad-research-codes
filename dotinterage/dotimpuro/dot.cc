#include <cmath>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "giv.h"

//Esse é o programa em que há impureza.

//const double KB = 8.61734315e-5;
const double KB = 1.;
const double GAMMA = 1;
const double EPS = 1e-14;

using namespace std;

int ind(int i, int j)  //Função que dá os índices de uma matriz em uma forma Hermitiana
{
  int m;

  if (i >= j)
    {
      m = i*(i-1)/int(2) + j;
    }

  else
    {
      m = j*(j-1)/int(2) + i;
      //Podemos trocar a ordem dos índices porque a matriz manipulada é Hermitiana.
    }

  return m;
}

double fatorial(int N)
{
  double product = 1.;
  for(int i = N; i > 1; i--)
  {
    product = product * double(i);
  }
  
  return product;
}  

double matrix_trace(double matrix[], int dim)
{
  double tr = 0;
  for(int i = 1; i <= dim; i++)
  {
    tr = tr + matrix[ind(i, i)];
  }
 return tr;
}
    
  

double binomial(double N, double m)
{
  return fatorial(int(N))/(fatorial(int(N-m)) * fatorial(int(m)) );
}  
  

double number_of_states(int dim, int charge, int spin_z)
{
//Função utilizada para calcular o número de estados de um subespaço com carga e componente do spin no eixo z bem definidos.
//dim: dimensão do problema sem considerar o spin, no espaço de 1 partícula. A dimensão real do problema de 1 partícula será 2*dim
//charge: número de partículas presentes no sistema
//spin_up: valor da componente z do spin. O usuário deve entrar o valor de 2*Sz, supondo que, para um elétron, Sz = +- 1/2.

//Se os spins estivessem todos alinhados para cima, charge = Sz (fica subentendido o fator 2 multiplicando Sz).
//A cada partícula com spin para baixo, este número deve ser reduzido de 2.
                                   
//Seja m o número de partículas com spin up
//Seja n o número de partículas com spin down

//m + n = Q
//m - n = Sz

//m = (Q+Sz)/2
//n = (Q-Sz)/2

double m = (charge + spin_z)/2.;
double n = (charge - spin_z)/2.;

  double product1 = binomial(double(dim), m);
    
  double product2 = binomial(double(dim), n);
  
  return product1*product2;
}

double partition_function(double avaltotal[], int dimtotal, double temp)
{
  double Z = 0;
  for(int i = 1; i <= dimtotal; i++)
  {
    double aux = exp(-avaltotal[i]/(KB*temp));
    if(aux > EPS)
    {
    Z = Z + aux;
    }
   
  }
  return Z;
}

double media_spin2(double spinz[], double avaltotal[], int dimtotal, double temp)
{
  double sum = 0;
  for(int i = 1; i <= dimtotal; i++)
  {
    double aux = exp(-avaltotal[i]/(KB*temp));
    if(aux > EPS)
    {
    sum = sum + aux*spinz[i]*spinz[i];
    }
    //cout << avaltotal[i] << "\t" << spinz[i] << endl;
  }
  
  sum = sum/partition_function(avaltotal, dimtotal, temp);
  return sum;
}

int main(void)
{
  for(double W = -5; W <= 5; W = W + 0.01)
  {
  
  double Ed = -0.1;
  double U = -2.*Ed;
  double V = 0.05;
  double t = 1;
  int dim = 4;  //dimensão do problema de 1 corpo
  
  int dim11 = int (number_of_states(dim, 1, 1));
  int dim22 = int (number_of_states(dim, 2, 2));
  int dim20 = int (number_of_states(dim, 2, 0));
  int dim33 = int (number_of_states(dim, 3, 3));
  int dim31 = int (number_of_states(dim, 3, 1));
  int dim44 = int (number_of_states(dim, 4, 4));
  int dim42 = int (number_of_states(dim, 4, 2));
  int dim40 = int (number_of_states(dim, 4, 0));
  int dim53 = int (number_of_states(dim, 5, 3));
  int dim51 = int (number_of_states(dim, 5, 1));
  int dim62 = int (number_of_states(dim, 6, 2));
  int dim60 = int (number_of_states(dim, 6, 0));
  int dim71 = int (number_of_states(dim, 7, 1));
  int dim80 = int (number_of_states(dim, 8, 0));
  
  int dimtotal = 2*(dim11 + dim22 + dim33 + dim31 + dim44 + dim42 + dim53 + dim51 + dim62 + dim71) + (dim20 + dim40 + dim60 + dim80) + 1;
  cout << "A dimensão total do sistema eh " << dimtotal << "." << endl;
  
  double *avaltotal = new double[dimtotal + 1];
  
  double *spinz = new double[dimtotal + 1];
  //Cada autovalor tem um Spin_z correspondente, que será guardado nesta matriz.
  
  //Nas expressões acima, o primeiro índice indica a carga; o segundo índice indica 2*Spin_z.
  
  
  double *H11 = new double[ind(dim11, dim11) + 1]; //Essa matriz guardará os elementos do Hamiltoniano.
  double *aval11 = new double[int(dim11 + 1)]; //Essa matriz guardará os autovalores do Hamiltoniano.
  double **avet11 = new double*[dim11+1];
  for(int lin = 1; lin <= dim11; lin++)
  {
    avet11[lin] = new double[dim11+1];
  }
  
  double *H22 = new double[ind(dim22, dim22) + 1]; //Essa matriz guardará os elementos do Hamiltoniano.
  double *aval22 = new double[int(dim22 + 1)]; //Essa matriz guardará os autovalores do Hamiltoniano.
  double **avet22 = new double*[dim22+1];
  for(int lin = 1; lin <= dim22; lin++)
  {
    avet22[lin] = new double[dim22+1];
  }
  
  double *H20 = new double[ind(dim20, dim20) + 1]; //Essa matriz guardará os elementos do Hamiltoniano.
  double *aval20 = new double[int(dim20 + 1)]; //Essa matriz guardará os autovalores do Hamiltoniano.
  double **avet20 = new double*[dim20+1];
  for(int lin = 1; lin <= dim20; lin++)
  {
    avet20[lin] = new double[dim20+1];
  }
  
  double *H33 = new double[ind(dim33, dim33) + 1]; //Essa matriz guardará os elementos do Hamiltoniano.
  double *aval33 = new double[int(dim33 + 1)]; //Essa matriz guardará os autovalores do Hamiltoniano.
  double **avet33 = new double*[dim33+1];
  for(int lin = 1; lin <= dim33; lin++)
  {
    avet33[lin] = new double[dim33+1];
  }
  
  double *H31 = new double[ind(dim31, dim31) + 1]; //Essa matriz guardará os elementos do Hamiltoniano.
  double *aval31 = new double[int(dim31 + 1)]; //Essa matriz guardará os autovalores do Hamiltoniano.
  double **avet31 = new double*[dim31+1];
  for(int lin = 1; lin <= dim31; lin++)
  {
    avet31[lin] = new double[dim31+1];
  }
  
  double *H44 = new double[ind(dim44, dim44) + 1]; //Essa matriz guardará os elementos do Hamiltoniano.
  double *aval44 = new double[int(dim44 + 1)]; //Essa matriz guardará os autovalores do Hamiltoniano.
  double **avet44 = new double*[dim44+1];
  for(int lin = 1; lin <= dim44; lin++)
  {
    avet44[lin] = new double[dim44+1];
  }
  
  double *H42 = new double[ind(dim42, dim42) + 1]; //Essa matriz guardará os elementos do Hamiltoniano.
  double *aval42 = new double[int(dim42 + 1)]; //Essa matriz guardará os autovalores do Hamiltoniano.
  double **avet42 = new double*[dim42+1];
  for(int lin = 1; lin <= dim42; lin++)
  {
    avet42[lin] = new double[dim42+1];
  }
  
  double *H40 = new double[ind(dim40, dim40) + 1]; //Essa matriz guardará os elementos do Hamiltoniano.
  double *aval40 = new double[int(dim40 + 1)]; //Essa matriz guardará os autovalores do Hamiltoniano.
  double **avet40 = new double*[dim40+1];
  for(int lin = 1; lin <= dim40; lin++)
  {
    avet40[lin] = new double[dim40+1];
  }
  
  double *H53 = new double[ind(dim53, dim53) + 1]; //Essa matriz guardará os elementos do Hamiltoniano.
  double *aval53 = new double[int(dim53 + 1)]; //Essa matriz guardará os autovalores do Hamiltoniano.
  double **avet53 = new double*[dim53+1];
  for(int lin = 1; lin <= dim53; lin++)
  {
    avet53[lin] = new double[dim53+1];
  }
  
  double *H51 = new double[ind(dim51, dim51) + 1]; //Essa matriz guardará os elementos do Hamiltoniano.
  double *aval51 = new double[int(dim51 + 1)]; //Essa matriz guardará os autovalores do Hamiltoniano.
  double **avet51 = new double*[dim51+1];
  for(int lin = 1; lin <= dim51; lin++)
  {
    avet51[lin] = new double[dim51+1];
  }
  
  double *H62 = new double[ind(dim62, dim62) + 1]; //Essa matriz guardará os elementos do Hamiltoniano.
  double *aval62 = new double[int(dim62 + 1)]; //Essa matriz guardará os autovalores do Hamiltoniano.
  double **avet62 = new double*[dim62+1];
  for(int lin = 1; lin <= dim62; lin++)
  {
    avet62[lin] = new double[dim62+1];
  }
  
  double *H60 = new double[ind(dim60, dim60) + 1]; //Essa matriz guardará os elementos do Hamiltoniano.
  double *aval60 = new double[int(dim60 + 1)]; //Essa matriz guardará os autovalores do Hamiltoniano.
  double **avet60 = new double*[dim60+1];
  for(int lin = 1; lin <= dim60; lin++)
  {
    avet60[lin] = new double[dim60+1];
  }
  
  double *H71 = new double[ind(dim71, dim71) + 1]; //Essa matriz guardará os elementos do Hamiltoniano.
  double *aval71 = new double[int(dim71 + 1)]; //Essa matriz guardará os autovalores do Hamiltoniano.
  double **avet71 = new double*[dim71+1];
  for(int lin = 1; lin <= dim71; lin++)
  {
    avet71[lin] = new double[dim71+1];
  }
  
  double *H80 = new double[ind(dim80, dim80) + 1];
  double *aval80 = new double[int(dim80 + 1)];
  double **avet80 = new double*[dim80+1];
  for(int lin = 1; lin <= dim80; lin++)
  {
    avet80[lin] = new double[dim80+1];
  }

H11[ind(1, 1)] = Ed; 
H11[ind(1, 2)] = V; 
H11[ind(1, 3)] = 0; 
H11[ind(1, 4)] = 0; 
H11[ind(2, 2)] = 0; 
H11[ind(2, 3)] = t; 
H11[ind(2, 4)] = 0; 
H11[ind(3, 3)] = 0; 
H11[ind(3, 4)] = t; 
H11[ind(4, 4)] = W; 

H22[ind(1, 1)] = Ed; 
H22[ind(1, 2)] = t; 
H22[ind(1, 3)] = 0; 
H22[ind(1, 4)] = 0; 
H22[ind(1, 5)] = 0; 
H22[ind(1, 6)] = 0; 
H22[ind(2, 2)] = Ed; 
H22[ind(2, 3)] = t; 
H22[ind(2, 4)] = V; 
H22[ind(2, 5)] = 0; 
H22[ind(2, 6)] = 0; 
H22[ind(3, 3)] = Ed+W; 
H22[ind(3, 4)] = 0; 
H22[ind(3, 5)] = V; 
H22[ind(3, 6)] = 0; 
H22[ind(4, 4)] = 0; 
H22[ind(4, 5)] = t; 
H22[ind(4, 6)] = 0; 
H22[ind(5, 5)] = W; 
H22[ind(5, 6)] = t; 
H22[ind(6, 6)] = W; 

H20[ind(1, 1)] = 2*Ed+U; 
H20[ind(1, 2)] = V; 
H20[ind(1, 3)] = 0; 
H20[ind(1, 4)] = 0; 
H20[ind(1, 5)] = V; 
H20[ind(1, 6)] = 0; 
H20[ind(1, 7)] = 0; 
H20[ind(1, 8)] = 0; 
H20[ind(1, 9)] = 0; 
H20[ind(1, 10)] = 0; 
H20[ind(1, 11)] = 0; 
H20[ind(1, 12)] = 0; 
H20[ind(1, 13)] = 0; 
H20[ind(1, 14)] = 0; 
H20[ind(1, 15)] = 0; 
H20[ind(1, 16)] = 0; 
H20[ind(2, 2)] = Ed; 
H20[ind(2, 3)] = t; 
H20[ind(2, 4)] = 0; 
H20[ind(2, 5)] = 0; 
H20[ind(2, 6)] = V; 
H20[ind(2, 7)] = 0; 
H20[ind(2, 8)] = 0; 
H20[ind(2, 9)] = 0; 
H20[ind(2, 10)] = 0; 
H20[ind(2, 11)] = 0; 
H20[ind(2, 12)] = 0; 
H20[ind(2, 13)] = 0; 
H20[ind(2, 14)] = 0; 
H20[ind(2, 15)] = 0; 
H20[ind(2, 16)] = 0; 
H20[ind(3, 3)] = Ed; 
H20[ind(3, 4)] = t; 
H20[ind(3, 5)] = 0; 
H20[ind(3, 6)] = 0; 
H20[ind(3, 7)] = V; 
H20[ind(3, 8)] = 0; 
H20[ind(3, 9)] = 0; 
H20[ind(3, 10)] = 0; 
H20[ind(3, 11)] = 0; 
H20[ind(3, 12)] = 0; 
H20[ind(3, 13)] = 0; 
H20[ind(3, 14)] = 0; 
H20[ind(3, 15)] = 0; 
H20[ind(3, 16)] = 0; 
H20[ind(4, 4)] = Ed+W; 
H20[ind(4, 5)] = 0; 
H20[ind(4, 6)] = 0; 
H20[ind(4, 7)] = 0; 
H20[ind(4, 8)] = V; 
H20[ind(4, 9)] = 0; 
H20[ind(4, 10)] = 0; 
H20[ind(4, 11)] = 0; 
H20[ind(4, 12)] = 0; 
H20[ind(4, 13)] = 0; 
H20[ind(4, 14)] = 0; 
H20[ind(4, 15)] = 0; 
H20[ind(4, 16)] = 0; 
H20[ind(5, 5)] = Ed; 
H20[ind(5, 6)] = V; 
H20[ind(5, 7)] = 0; 
H20[ind(5, 8)] = 0; 
H20[ind(5, 9)] = t; 
H20[ind(5, 10)] = 0; 
H20[ind(5, 11)] = 0; 
H20[ind(5, 12)] = 0; 
H20[ind(5, 13)] = 0; 
H20[ind(5, 14)] = 0; 
H20[ind(5, 15)] = 0; 
H20[ind(5, 16)] = 0; 
H20[ind(6, 6)] = 0; 
H20[ind(6, 7)] = t; 
H20[ind(6, 8)] = 0; 
H20[ind(6, 9)] = 0; 
H20[ind(6, 10)] = t; 
H20[ind(6, 11)] = 0; 
H20[ind(6, 12)] = 0; 
H20[ind(6, 13)] = 0; 
H20[ind(6, 14)] = 0; 
H20[ind(6, 15)] = 0; 
H20[ind(6, 16)] = 0; 
H20[ind(7, 7)] = 0; 
H20[ind(7, 8)] = t; 
H20[ind(7, 9)] = 0; 
H20[ind(7, 10)] = 0; 
H20[ind(7, 11)] = t; 
H20[ind(7, 12)] = 0; 
H20[ind(7, 13)] = 0; 
H20[ind(7, 14)] = 0; 
H20[ind(7, 15)] = 0; 
H20[ind(7, 16)] = 0; 
H20[ind(8, 8)] = W; 
H20[ind(8, 9)] = 0; 
H20[ind(8, 10)] = 0; 
H20[ind(8, 11)] = 0; 
H20[ind(8, 12)] = t; 
H20[ind(8, 13)] = 0; 
H20[ind(8, 14)] = 0; 
H20[ind(8, 15)] = 0; 
H20[ind(8, 16)] = 0; 
H20[ind(9, 9)] = Ed; 
H20[ind(9, 10)] = V; 
H20[ind(9, 11)] = 0; 
H20[ind(9, 12)] = 0; 
H20[ind(9, 13)] = t; 
H20[ind(9, 14)] = 0; 
H20[ind(9, 15)] = 0; 
H20[ind(9, 16)] = 0; 
H20[ind(10, 10)] = 0; 
H20[ind(10, 11)] = t; 
H20[ind(10, 12)] = 0; 
H20[ind(10, 13)] = 0; 
H20[ind(10, 14)] = t; 
H20[ind(10, 15)] = 0; 
H20[ind(10, 16)] = 0; 
H20[ind(11, 11)] = 0; 
H20[ind(11, 12)] = t; 
H20[ind(11, 13)] = 0; 
H20[ind(11, 14)] = 0; 
H20[ind(11, 15)] = t; 
H20[ind(11, 16)] = 0; 
H20[ind(12, 12)] = W; 
H20[ind(12, 13)] = 0; 
H20[ind(12, 14)] = 0; 
H20[ind(12, 15)] = 0; 
H20[ind(12, 16)] = t; 
H20[ind(13, 13)] = Ed+W; 
H20[ind(13, 14)] = V; 
H20[ind(13, 15)] = 0; 
H20[ind(13, 16)] = 0; 
H20[ind(14, 14)] = W; 
H20[ind(14, 15)] = t; 
H20[ind(14, 16)] = 0; 
H20[ind(15, 15)] = W; 
H20[ind(15, 16)] = t; 
H20[ind(16, 16)] = 2*W; 

H33[ind(1, 1)] = Ed; 
H33[ind(1, 2)] = t; 
H33[ind(1, 3)] = 0; 
H33[ind(1, 4)] = 0; 
H33[ind(2, 2)] = Ed+W; 
H33[ind(2, 3)] = t; 
H33[ind(2, 4)] = 0; 
H33[ind(3, 3)] = Ed+W; 
H33[ind(3, 4)] = V; 
H33[ind(4, 4)] = W; 

H31[ind(1, 1)] = 2*Ed+U; 
H31[ind(1, 2)] = V; 
H31[ind(1, 3)] = 0; 
H31[ind(1, 4)] = 0; 
H31[ind(1, 5)] = t; 
H31[ind(1, 6)] = 0; 
H31[ind(1, 7)] = 0; 
H31[ind(1, 8)] = 0; 
H31[ind(1, 9)] = 0; 
H31[ind(1, 10)] = 0; 
H31[ind(1, 11)] = 0; 
H31[ind(1, 12)] = 0; 
H31[ind(1, 13)] = 0; 
H31[ind(1, 14)] = 0; 
H31[ind(1, 15)] = 0; 
H31[ind(1, 16)] = 0; 
H31[ind(1, 17)] = 0; 
H31[ind(1, 18)] = 0; 
H31[ind(1, 19)] = 0; 
H31[ind(1, 20)] = 0; 
H31[ind(1, 21)] = 0; 
H31[ind(1, 22)] = 0; 
H31[ind(1, 23)] = 0; 
H31[ind(1, 24)] = 0; 
H31[ind(2, 2)] = Ed; 
H31[ind(2, 3)] = t; 
H31[ind(2, 4)] = 0; 
H31[ind(2, 5)] = 0; 
H31[ind(2, 6)] = t; 
H31[ind(2, 7)] = 0; 
H31[ind(2, 8)] = 0; 
H31[ind(2, 9)] = 0; 
H31[ind(2, 10)] = 0; 
H31[ind(2, 11)] = 0; 
H31[ind(2, 12)] = 0; 
H31[ind(2, 13)] = 0; 
H31[ind(2, 14)] = 0; 
H31[ind(2, 15)] = 0; 
H31[ind(2, 16)] = 0; 
H31[ind(2, 17)] = 0; 
H31[ind(2, 18)] = 0; 
H31[ind(2, 19)] = 0; 
H31[ind(2, 20)] = 0; 
H31[ind(2, 21)] = 0; 
H31[ind(2, 22)] = 0; 
H31[ind(2, 23)] = 0; 
H31[ind(2, 24)] = 0; 
H31[ind(3, 3)] = Ed; 
H31[ind(3, 4)] = t; 
H31[ind(3, 5)] = 0; 
H31[ind(3, 6)] = 0; 
H31[ind(3, 7)] = t; 
H31[ind(3, 8)] = 0; 
H31[ind(3, 9)] = 0; 
H31[ind(3, 10)] = 0; 
H31[ind(3, 11)] = 0; 
H31[ind(3, 12)] = 0; 
H31[ind(3, 13)] = 0; 
H31[ind(3, 14)] = 0; 
H31[ind(3, 15)] = 0; 
H31[ind(3, 16)] = 0; 
H31[ind(3, 17)] = 0; 
H31[ind(3, 18)] = 0; 
H31[ind(3, 19)] = 0; 
H31[ind(3, 20)] = 0; 
H31[ind(3, 21)] = 0; 
H31[ind(3, 22)] = 0; 
H31[ind(3, 23)] = 0; 
H31[ind(3, 24)] = 0; 
H31[ind(4, 4)] = Ed+W; 
H31[ind(4, 5)] = 0; 
H31[ind(4, 6)] = 0; 
H31[ind(4, 7)] = 0; 
H31[ind(4, 8)] = t; 
H31[ind(4, 9)] = 0; 
H31[ind(4, 10)] = 0; 
H31[ind(4, 11)] = 0; 
H31[ind(4, 12)] = 0; 
H31[ind(4, 13)] = 0; 
H31[ind(4, 14)] = 0; 
H31[ind(4, 15)] = 0; 
H31[ind(4, 16)] = 0; 
H31[ind(4, 17)] = 0; 
H31[ind(4, 18)] = 0; 
H31[ind(4, 19)] = 0; 
H31[ind(4, 20)] = 0; 
H31[ind(4, 21)] = 0; 
H31[ind(4, 22)] = 0; 
H31[ind(4, 23)] = 0; 
H31[ind(4, 24)] = 0; 
H31[ind(5, 5)] = 2*Ed+U; 
H31[ind(5, 6)] = V; 
H31[ind(5, 7)] = 0; 
H31[ind(5, 8)] = 0; 
H31[ind(5, 9)] = t; 
H31[ind(5, 10)] = 0; 
H31[ind(5, 11)] = 0; 
H31[ind(5, 12)] = 0; 
H31[ind(5, 13)] = V; 
H31[ind(5, 14)] = 0; 
H31[ind(5, 15)] = 0; 
H31[ind(5, 16)] = 0; 
H31[ind(5, 17)] = 0; 
H31[ind(5, 18)] = 0; 
H31[ind(5, 19)] = 0; 
H31[ind(5, 20)] = 0; 
H31[ind(5, 21)] = 0; 
H31[ind(5, 22)] = 0; 
H31[ind(5, 23)] = 0; 
H31[ind(5, 24)] = 0; 
H31[ind(6, 6)] = Ed; 
H31[ind(6, 7)] = t; 
H31[ind(6, 8)] = 0; 
H31[ind(6, 9)] = 0; 
H31[ind(6, 10)] = t; 
H31[ind(6, 11)] = 0; 
H31[ind(6, 12)] = 0; 
H31[ind(6, 13)] = 0; 
H31[ind(6, 14)] = V; 
H31[ind(6, 15)] = 0; 
H31[ind(6, 16)] = 0; 
H31[ind(6, 17)] = 0; 
H31[ind(6, 18)] = 0; 
H31[ind(6, 19)] = 0; 
H31[ind(6, 20)] = 0; 
H31[ind(6, 21)] = 0; 
H31[ind(6, 22)] = 0; 
H31[ind(6, 23)] = 0; 
H31[ind(6, 24)] = 0; 
H31[ind(7, 7)] = Ed; 
H31[ind(7, 8)] = t; 
H31[ind(7, 9)] = 0; 
H31[ind(7, 10)] = 0; 
H31[ind(7, 11)] = t; 
H31[ind(7, 12)] = 0; 
H31[ind(7, 13)] = 0; 
H31[ind(7, 14)] = 0; 
H31[ind(7, 15)] = V; 
H31[ind(7, 16)] = 0; 
H31[ind(7, 17)] = 0; 
H31[ind(7, 18)] = 0; 
H31[ind(7, 19)] = 0; 
H31[ind(7, 20)] = 0; 
H31[ind(7, 21)] = 0; 
H31[ind(7, 22)] = 0; 
H31[ind(7, 23)] = 0; 
H31[ind(7, 24)] = 0; 
H31[ind(8, 8)] = Ed+W; 
H31[ind(8, 9)] = 0; 
H31[ind(8, 10)] = 0; 
H31[ind(8, 11)] = 0; 
H31[ind(8, 12)] = t; 
H31[ind(8, 13)] = 0; 
H31[ind(8, 14)] = 0; 
H31[ind(8, 15)] = 0; 
H31[ind(8, 16)] = V; 
H31[ind(8, 17)] = 0; 
H31[ind(8, 18)] = 0; 
H31[ind(8, 19)] = 0; 
H31[ind(8, 20)] = 0; 
H31[ind(8, 21)] = 0; 
H31[ind(8, 22)] = 0; 
H31[ind(8, 23)] = 0; 
H31[ind(8, 24)] = 0; 
H31[ind(9, 9)] = 2*Ed+U+W; 
H31[ind(9, 10)] = V; 
H31[ind(9, 11)] = 0; 
H31[ind(9, 12)] = 0; 
H31[ind(9, 13)] = 0; 
H31[ind(9, 14)] = 0; 
H31[ind(9, 15)] = 0; 
H31[ind(9, 16)] = 0; 
H31[ind(9, 17)] = V; 
H31[ind(9, 18)] = 0; 
H31[ind(9, 19)] = 0; 
H31[ind(9, 20)] = 0; 
H31[ind(9, 21)] = 0; 
H31[ind(9, 22)] = 0; 
H31[ind(9, 23)] = 0; 
H31[ind(9, 24)] = 0; 
H31[ind(10, 10)] = Ed+W; 
H31[ind(10, 11)] = t; 
H31[ind(10, 12)] = 0; 
H31[ind(10, 13)] = 0; 
H31[ind(10, 14)] = 0; 
H31[ind(10, 15)] = 0; 
H31[ind(10, 16)] = 0; 
H31[ind(10, 17)] = 0; 
H31[ind(10, 18)] = V; 
H31[ind(10, 19)] = 0; 
H31[ind(10, 20)] = 0; 
H31[ind(10, 21)] = 0; 
H31[ind(10, 22)] = 0; 
H31[ind(10, 23)] = 0; 
H31[ind(10, 24)] = 0; 
H31[ind(11, 11)] = Ed+W; 
H31[ind(11, 12)] = t; 
H31[ind(11, 13)] = 0; 
H31[ind(11, 14)] = 0; 
H31[ind(11, 15)] = 0; 
H31[ind(11, 16)] = 0; 
H31[ind(11, 17)] = 0; 
H31[ind(11, 18)] = 0; 
H31[ind(11, 19)] = V; 
H31[ind(11, 20)] = 0; 
H31[ind(11, 21)] = 0; 
H31[ind(11, 22)] = 0; 
H31[ind(11, 23)] = 0; 
H31[ind(11, 24)] = 0; 
H31[ind(12, 12)] = Ed+2*W; 
H31[ind(12, 13)] = 0; 
H31[ind(12, 14)] = 0; 
H31[ind(12, 15)] = 0; 
H31[ind(12, 16)] = 0; 
H31[ind(12, 17)] = 0; 
H31[ind(12, 18)] = 0; 
H31[ind(12, 19)] = 0; 
H31[ind(12, 20)] = V; 
H31[ind(12, 21)] = 0; 
H31[ind(12, 22)] = 0; 
H31[ind(12, 23)] = 0; 
H31[ind(12, 24)] = 0; 
H31[ind(13, 13)] = Ed; 
H31[ind(13, 14)] = V; 
H31[ind(13, 15)] = 0; 
H31[ind(13, 16)] = 0; 
H31[ind(13, 17)] = t; 
H31[ind(13, 18)] = 0; 
H31[ind(13, 19)] = 0; 
H31[ind(13, 20)] = 0; 
H31[ind(13, 21)] = 0; 
H31[ind(13, 22)] = 0; 
H31[ind(13, 23)] = 0; 
H31[ind(13, 24)] = 0; 
H31[ind(14, 14)] = 0; 
H31[ind(14, 15)] = t; 
H31[ind(14, 16)] = 0; 
H31[ind(14, 17)] = 0; 
H31[ind(14, 18)] = t; 
H31[ind(14, 19)] = 0; 
H31[ind(14, 20)] = 0; 
H31[ind(14, 21)] = 0; 
H31[ind(14, 22)] = 0; 
H31[ind(14, 23)] = 0; 
H31[ind(14, 24)] = 0; 
H31[ind(15, 15)] = 0; 
H31[ind(15, 16)] = t; 
H31[ind(15, 17)] = 0; 
H31[ind(15, 18)] = 0; 
H31[ind(15, 19)] = t; 
H31[ind(15, 20)] = 0; 
H31[ind(15, 21)] = 0; 
H31[ind(15, 22)] = 0; 
H31[ind(15, 23)] = 0; 
H31[ind(15, 24)] = 0; 
H31[ind(16, 16)] = W; 
H31[ind(16, 17)] = 0; 
H31[ind(16, 18)] = 0; 
H31[ind(16, 19)] = 0; 
H31[ind(16, 20)] = t; 
H31[ind(16, 21)] = 0; 
H31[ind(16, 22)] = 0; 
H31[ind(16, 23)] = 0; 
H31[ind(16, 24)] = 0; 
H31[ind(17, 17)] = Ed+W; 
H31[ind(17, 18)] = V; 
H31[ind(17, 19)] = 0; 
H31[ind(17, 20)] = 0; 
H31[ind(17, 21)] = t; 
H31[ind(17, 22)] = 0; 
H31[ind(17, 23)] = 0; 
H31[ind(17, 24)] = 0; 
H31[ind(18, 18)] = W; 
H31[ind(18, 19)] = t; 
H31[ind(18, 20)] = 0; 
H31[ind(18, 21)] = 0; 
H31[ind(18, 22)] = t; 
H31[ind(18, 23)] = 0; 
H31[ind(18, 24)] = 0; 
H31[ind(19, 19)] = W; 
H31[ind(19, 20)] = t; 
H31[ind(19, 21)] = 0; 
H31[ind(19, 22)] = 0; 
H31[ind(19, 23)] = t; 
H31[ind(19, 24)] = 0; 
H31[ind(20, 20)] = 2*W; 
H31[ind(20, 21)] = 0; 
H31[ind(20, 22)] = 0; 
H31[ind(20, 23)] = 0; 
H31[ind(20, 24)] = t; 
H31[ind(21, 21)] = Ed+W; 
H31[ind(21, 22)] = V; 
H31[ind(21, 23)] = 0; 
H31[ind(21, 24)] = 0; 
H31[ind(22, 22)] = W; 
H31[ind(22, 23)] = t; 
H31[ind(22, 24)] = 0; 
H31[ind(23, 23)] = W; 
H31[ind(23, 24)] = t; 
H31[ind(24, 24)] = 2*W; 

H44[ind(1, 1)] = Ed+W; 

H42[ind(1, 1)] = 2*Ed+U; 
H42[ind(1, 2)] = V; 
H42[ind(1, 3)] = 0; 
H42[ind(1, 4)] = 0; 
H42[ind(1, 5)] = t; 
H42[ind(1, 6)] = 0; 
H42[ind(1, 7)] = 0; 
H42[ind(1, 8)] = 0; 
H42[ind(1, 9)] = 0; 
H42[ind(1, 10)] = 0; 
H42[ind(1, 11)] = 0; 
H42[ind(1, 12)] = 0; 
H42[ind(1, 13)] = 0; 
H42[ind(1, 14)] = 0; 
H42[ind(1, 15)] = 0; 
H42[ind(1, 16)] = 0; 
H42[ind(2, 2)] = Ed; 
H42[ind(2, 3)] = t; 
H42[ind(2, 4)] = 0; 
H42[ind(2, 5)] = 0; 
H42[ind(2, 6)] = t; 
H42[ind(2, 7)] = 0; 
H42[ind(2, 8)] = 0; 
H42[ind(2, 9)] = 0; 
H42[ind(2, 10)] = 0; 
H42[ind(2, 11)] = 0; 
H42[ind(2, 12)] = 0; 
H42[ind(2, 13)] = 0; 
H42[ind(2, 14)] = 0; 
H42[ind(2, 15)] = 0; 
H42[ind(2, 16)] = 0; 
H42[ind(3, 3)] = Ed; 
H42[ind(3, 4)] = t; 
H42[ind(3, 5)] = 0; 
H42[ind(3, 6)] = 0; 
H42[ind(3, 7)] = t; 
H42[ind(3, 8)] = 0; 
H42[ind(3, 9)] = 0; 
H42[ind(3, 10)] = 0; 
H42[ind(3, 11)] = 0; 
H42[ind(3, 12)] = 0; 
H42[ind(3, 13)] = 0; 
H42[ind(3, 14)] = 0; 
H42[ind(3, 15)] = 0; 
H42[ind(3, 16)] = 0; 
H42[ind(4, 4)] = Ed+W; 
H42[ind(4, 5)] = 0; 
H42[ind(4, 6)] = 0; 
H42[ind(4, 7)] = 0; 
H42[ind(4, 8)] = t; 
H42[ind(4, 9)] = 0; 
H42[ind(4, 10)] = 0; 
H42[ind(4, 11)] = 0; 
H42[ind(4, 12)] = 0; 
H42[ind(4, 13)] = 0; 
H42[ind(4, 14)] = 0; 
H42[ind(4, 15)] = 0; 
H42[ind(4, 16)] = 0; 
H42[ind(5, 5)] = 2*Ed+U+W; 
H42[ind(5, 6)] = V; 
H42[ind(5, 7)] = 0; 
H42[ind(5, 8)] = 0; 
H42[ind(5, 9)] = t; 
H42[ind(5, 10)] = 0; 
H42[ind(5, 11)] = 0; 
H42[ind(5, 12)] = 0; 
H42[ind(5, 13)] = 0; 
H42[ind(5, 14)] = 0; 
H42[ind(5, 15)] = 0; 
H42[ind(5, 16)] = 0; 
H42[ind(6, 6)] = Ed+W; 
H42[ind(6, 7)] = t; 
H42[ind(6, 8)] = 0; 
H42[ind(6, 9)] = 0; 
H42[ind(6, 10)] = t; 
H42[ind(6, 11)] = 0; 
H42[ind(6, 12)] = 0; 
H42[ind(6, 13)] = 0; 
H42[ind(6, 14)] = 0; 
H42[ind(6, 15)] = 0; 
H42[ind(6, 16)] = 0; 
H42[ind(7, 7)] = Ed+W; 
H42[ind(7, 8)] = t; 
H42[ind(7, 9)] = 0; 
H42[ind(7, 10)] = 0; 
H42[ind(7, 11)] = t; 
H42[ind(7, 12)] = 0; 
H42[ind(7, 13)] = 0; 
H42[ind(7, 14)] = 0; 
H42[ind(7, 15)] = 0; 
H42[ind(7, 16)] = 0; 
H42[ind(8, 8)] = Ed+2*W; 
H42[ind(8, 9)] = 0; 
H42[ind(8, 10)] = 0; 
H42[ind(8, 11)] = 0; 
H42[ind(8, 12)] = t; 
H42[ind(8, 13)] = 0; 
H42[ind(8, 14)] = 0; 
H42[ind(8, 15)] = 0; 
H42[ind(8, 16)] = 0; 
H42[ind(9, 9)] = 2*Ed+U+W; 
H42[ind(9, 10)] = V; 
H42[ind(9, 11)] = 0; 
H42[ind(9, 12)] = 0; 
H42[ind(9, 13)] = V; 
H42[ind(9, 14)] = 0; 
H42[ind(9, 15)] = 0; 
H42[ind(9, 16)] = 0; 
H42[ind(10, 10)] = Ed+W; 
H42[ind(10, 11)] = t; 
H42[ind(10, 12)] = 0; 
H42[ind(10, 13)] = 0; 
H42[ind(10, 14)] = V; 
H42[ind(10, 15)] = 0; 
H42[ind(10, 16)] = 0; 
H42[ind(11, 11)] = Ed+W; 
H42[ind(11, 12)] = t; 
H42[ind(11, 13)] = 0; 
H42[ind(11, 14)] = 0; 
H42[ind(11, 15)] = V; 
H42[ind(11, 16)] = 0; 
H42[ind(12, 12)] = Ed+2*W; 
H42[ind(12, 13)] = 0; 
H42[ind(12, 14)] = 0; 
H42[ind(12, 15)] = 0; 
H42[ind(12, 16)] = V; 
H42[ind(13, 13)] = Ed+W; 
H42[ind(13, 14)] = V; 
H42[ind(13, 15)] = 0; 
H42[ind(13, 16)] = 0; 
H42[ind(14, 14)] = W; 
H42[ind(14, 15)] = t; 
H42[ind(14, 16)] = 0; 
H42[ind(15, 15)] = W; 
H42[ind(15, 16)] = t; 
H42[ind(16, 16)] = 2*W; 

H40[ind(1, 1)] = 2*Ed+U; 
H40[ind(1, 2)] = t; 
H40[ind(1, 3)] = 0; 
H40[ind(1, 4)] = 0; 
H40[ind(1, 5)] = 0; 
H40[ind(1, 6)] = 0; 
H40[ind(1, 7)] = t; 
H40[ind(1, 8)] = 0; 
H40[ind(1, 9)] = 0; 
H40[ind(1, 10)] = 0; 
H40[ind(1, 11)] = 0; 
H40[ind(1, 12)] = 0; 
H40[ind(1, 13)] = 0; 
H40[ind(1, 14)] = 0; 
H40[ind(1, 15)] = 0; 
H40[ind(1, 16)] = 0; 
H40[ind(1, 17)] = 0; 
H40[ind(1, 18)] = 0; 
H40[ind(1, 19)] = 0; 
H40[ind(1, 20)] = 0; 
H40[ind(1, 21)] = 0; 
H40[ind(1, 22)] = 0; 
H40[ind(1, 23)] = 0; 
H40[ind(1, 24)] = 0; 
H40[ind(1, 25)] = 0; 
H40[ind(1, 26)] = 0; 
H40[ind(1, 27)] = 0; 
H40[ind(1, 28)] = 0; 
H40[ind(1, 29)] = 0; 
H40[ind(1, 30)] = 0; 
H40[ind(1, 31)] = 0; 
H40[ind(1, 32)] = 0; 
H40[ind(1, 33)] = 0; 
H40[ind(1, 34)] = 0; 
H40[ind(1, 35)] = 0; 
H40[ind(1, 36)] = 0; 
H40[ind(2, 2)] = 2*Ed+U; 
H40[ind(2, 3)] = t; 
H40[ind(2, 4)] = V; 
H40[ind(2, 5)] = 0; 
H40[ind(2, 6)] = 0; 
H40[ind(2, 7)] = 0; 
H40[ind(2, 8)] = t; 
H40[ind(2, 9)] = 0; 
H40[ind(2, 10)] = 0; 
H40[ind(2, 11)] = 0; 
H40[ind(2, 12)] = 0; 
H40[ind(2, 13)] = 0; 
H40[ind(2, 14)] = 0; 
H40[ind(2, 15)] = 0; 
H40[ind(2, 16)] = 0; 
H40[ind(2, 17)] = 0; 
H40[ind(2, 18)] = 0; 
H40[ind(2, 19)] = 0; 
H40[ind(2, 20)] = 0; 
H40[ind(2, 21)] = 0; 
H40[ind(2, 22)] = 0; 
H40[ind(2, 23)] = 0; 
H40[ind(2, 24)] = 0; 
H40[ind(2, 25)] = 0; 
H40[ind(2, 26)] = 0; 
H40[ind(2, 27)] = 0; 
H40[ind(2, 28)] = 0; 
H40[ind(2, 29)] = 0; 
H40[ind(2, 30)] = 0; 
H40[ind(2, 31)] = 0; 
H40[ind(2, 32)] = 0; 
H40[ind(2, 33)] = 0; 
H40[ind(2, 34)] = 0; 
H40[ind(2, 35)] = 0; 
H40[ind(2, 36)] = 0; 
H40[ind(3, 3)] = 2*Ed+U+W; 
H40[ind(3, 4)] = 0; 
H40[ind(3, 5)] = V; 
H40[ind(3, 6)] = 0; 
H40[ind(3, 7)] = 0; 
H40[ind(3, 8)] = 0; 
H40[ind(3, 9)] = t; 
H40[ind(3, 10)] = 0; 
H40[ind(3, 11)] = 0; 
H40[ind(3, 12)] = 0; 
H40[ind(3, 13)] = 0; 
H40[ind(3, 14)] = 0; 
H40[ind(3, 15)] = 0; 
H40[ind(3, 16)] = 0; 
H40[ind(3, 17)] = 0; 
H40[ind(3, 18)] = 0; 
H40[ind(3, 19)] = 0; 
H40[ind(3, 20)] = 0; 
H40[ind(3, 21)] = 0; 
H40[ind(3, 22)] = 0; 
H40[ind(3, 23)] = 0; 
H40[ind(3, 24)] = 0; 
H40[ind(3, 25)] = 0; 
H40[ind(3, 26)] = 0; 
H40[ind(3, 27)] = 0; 
H40[ind(3, 28)] = 0; 
H40[ind(3, 29)] = 0; 
H40[ind(3, 30)] = 0; 
H40[ind(3, 31)] = 0; 
H40[ind(3, 32)] = 0; 
H40[ind(3, 33)] = 0; 
H40[ind(3, 34)] = 0; 
H40[ind(3, 35)] = 0; 
H40[ind(3, 36)] = 0; 
H40[ind(4, 4)] = Ed; 
H40[ind(4, 5)] = t; 
H40[ind(4, 6)] = 0; 
H40[ind(4, 7)] = 0; 
H40[ind(4, 8)] = 0; 
H40[ind(4, 9)] = 0; 
H40[ind(4, 10)] = t; 
H40[ind(4, 11)] = 0; 
H40[ind(4, 12)] = 0; 
H40[ind(4, 13)] = 0; 
H40[ind(4, 14)] = 0; 
H40[ind(4, 15)] = 0; 
H40[ind(4, 16)] = 0; 
H40[ind(4, 17)] = 0; 
H40[ind(4, 18)] = 0; 
H40[ind(4, 19)] = 0; 
H40[ind(4, 20)] = 0; 
H40[ind(4, 21)] = 0; 
H40[ind(4, 22)] = 0; 
H40[ind(4, 23)] = 0; 
H40[ind(4, 24)] = 0; 
H40[ind(4, 25)] = 0; 
H40[ind(4, 26)] = 0; 
H40[ind(4, 27)] = 0; 
H40[ind(4, 28)] = 0; 
H40[ind(4, 29)] = 0; 
H40[ind(4, 30)] = 0; 
H40[ind(4, 31)] = 0; 
H40[ind(4, 32)] = 0; 
H40[ind(4, 33)] = 0; 
H40[ind(4, 34)] = 0; 
H40[ind(4, 35)] = 0; 
H40[ind(4, 36)] = 0; 
H40[ind(5, 5)] = Ed+W; 
H40[ind(5, 6)] = t; 
H40[ind(5, 7)] = 0; 
H40[ind(5, 8)] = 0; 
H40[ind(5, 9)] = 0; 
H40[ind(5, 10)] = 0; 
H40[ind(5, 11)] = t; 
H40[ind(5, 12)] = 0; 
H40[ind(5, 13)] = 0; 
H40[ind(5, 14)] = 0; 
H40[ind(5, 15)] = 0; 
H40[ind(5, 16)] = 0; 
H40[ind(5, 17)] = 0; 
H40[ind(5, 18)] = 0; 
H40[ind(5, 19)] = 0; 
H40[ind(5, 20)] = 0; 
H40[ind(5, 21)] = 0; 
H40[ind(5, 22)] = 0; 
H40[ind(5, 23)] = 0; 
H40[ind(5, 24)] = 0; 
H40[ind(5, 25)] = 0; 
H40[ind(5, 26)] = 0; 
H40[ind(5, 27)] = 0; 
H40[ind(5, 28)] = 0; 
H40[ind(5, 29)] = 0; 
H40[ind(5, 30)] = 0; 
H40[ind(5, 31)] = 0; 
H40[ind(5, 32)] = 0; 
H40[ind(5, 33)] = 0; 
H40[ind(5, 34)] = 0; 
H40[ind(5, 35)] = 0; 
H40[ind(5, 36)] = 0; 
H40[ind(6, 6)] = Ed+W; 
H40[ind(6, 7)] = 0; 
H40[ind(6, 8)] = 0; 
H40[ind(6, 9)] = 0; 
H40[ind(6, 10)] = 0; 
H40[ind(6, 11)] = 0; 
H40[ind(6, 12)] = t; 
H40[ind(6, 13)] = 0; 
H40[ind(6, 14)] = 0; 
H40[ind(6, 15)] = 0; 
H40[ind(6, 16)] = 0; 
H40[ind(6, 17)] = 0; 
H40[ind(6, 18)] = 0; 
H40[ind(6, 19)] = 0; 
H40[ind(6, 20)] = 0; 
H40[ind(6, 21)] = 0; 
H40[ind(6, 22)] = 0; 
H40[ind(6, 23)] = 0; 
H40[ind(6, 24)] = 0; 
H40[ind(6, 25)] = 0; 
H40[ind(6, 26)] = 0; 
H40[ind(6, 27)] = 0; 
H40[ind(6, 28)] = 0; 
H40[ind(6, 29)] = 0; 
H40[ind(6, 30)] = 0; 
H40[ind(6, 31)] = 0; 
H40[ind(6, 32)] = 0; 
H40[ind(6, 33)] = 0; 
H40[ind(6, 34)] = 0; 
H40[ind(6, 35)] = 0; 
H40[ind(6, 36)] = 0; 
H40[ind(7, 7)] = 2*Ed+U; 
H40[ind(7, 8)] = t; 
H40[ind(7, 9)] = 0; 
H40[ind(7, 10)] = 0; 
H40[ind(7, 11)] = 0; 
H40[ind(7, 12)] = 0; 
H40[ind(7, 13)] = t; 
H40[ind(7, 14)] = 0; 
H40[ind(7, 15)] = 0; 
H40[ind(7, 16)] = 0; 
H40[ind(7, 17)] = 0; 
H40[ind(7, 18)] = 0; 
H40[ind(7, 19)] = V; 
H40[ind(7, 20)] = 0; 
H40[ind(7, 21)] = 0; 
H40[ind(7, 22)] = 0; 
H40[ind(7, 23)] = 0; 
H40[ind(7, 24)] = 0; 
H40[ind(7, 25)] = 0; 
H40[ind(7, 26)] = 0; 
H40[ind(7, 27)] = 0; 
H40[ind(7, 28)] = 0; 
H40[ind(7, 29)] = 0; 
H40[ind(7, 30)] = 0; 
H40[ind(7, 31)] = 0; 
H40[ind(7, 32)] = 0; 
H40[ind(7, 33)] = 0; 
H40[ind(7, 34)] = 0; 
H40[ind(7, 35)] = 0; 
H40[ind(7, 36)] = 0; 
H40[ind(8, 8)] = 2*Ed+U; 
H40[ind(8, 9)] = t; 
H40[ind(8, 10)] = V; 
H40[ind(8, 11)] = 0; 
H40[ind(8, 12)] = 0; 
H40[ind(8, 13)] = 0; 
H40[ind(8, 14)] = t; 
H40[ind(8, 15)] = 0; 
H40[ind(8, 16)] = 0; 
H40[ind(8, 17)] = 0; 
H40[ind(8, 18)] = 0; 
H40[ind(8, 19)] = 0; 
H40[ind(8, 20)] = V; 
H40[ind(8, 21)] = 0; 
H40[ind(8, 22)] = 0; 
H40[ind(8, 23)] = 0; 
H40[ind(8, 24)] = 0; 
H40[ind(8, 25)] = 0; 
H40[ind(8, 26)] = 0; 
H40[ind(8, 27)] = 0; 
H40[ind(8, 28)] = 0; 
H40[ind(8, 29)] = 0; 
H40[ind(8, 30)] = 0; 
H40[ind(8, 31)] = 0; 
H40[ind(8, 32)] = 0; 
H40[ind(8, 33)] = 0; 
H40[ind(8, 34)] = 0; 
H40[ind(8, 35)] = 0; 
H40[ind(8, 36)] = 0; 
H40[ind(9, 9)] = 2*Ed+U+W; 
H40[ind(9, 10)] = 0; 
H40[ind(9, 11)] = V; 
H40[ind(9, 12)] = 0; 
H40[ind(9, 13)] = 0; 
H40[ind(9, 14)] = 0; 
H40[ind(9, 15)] = t; 
H40[ind(9, 16)] = 0; 
H40[ind(9, 17)] = 0; 
H40[ind(9, 18)] = 0; 
H40[ind(9, 19)] = 0; 
H40[ind(9, 20)] = 0; 
H40[ind(9, 21)] = V; 
H40[ind(9, 22)] = 0; 
H40[ind(9, 23)] = 0; 
H40[ind(9, 24)] = 0; 
H40[ind(9, 25)] = 0; 
H40[ind(9, 26)] = 0; 
H40[ind(9, 27)] = 0; 
H40[ind(9, 28)] = 0; 
H40[ind(9, 29)] = 0; 
H40[ind(9, 30)] = 0; 
H40[ind(9, 31)] = 0; 
H40[ind(9, 32)] = 0; 
H40[ind(9, 33)] = 0; 
H40[ind(9, 34)] = 0; 
H40[ind(9, 35)] = 0; 
H40[ind(9, 36)] = 0; 
H40[ind(10, 10)] = Ed; 
H40[ind(10, 11)] = t; 
H40[ind(10, 12)] = 0; 
H40[ind(10, 13)] = 0; 
H40[ind(10, 14)] = 0; 
H40[ind(10, 15)] = 0; 
H40[ind(10, 16)] = t; 
H40[ind(10, 17)] = 0; 
H40[ind(10, 18)] = 0; 
H40[ind(10, 19)] = 0; 
H40[ind(10, 20)] = 0; 
H40[ind(10, 21)] = 0; 
H40[ind(10, 22)] = V; 
H40[ind(10, 23)] = 0; 
H40[ind(10, 24)] = 0; 
H40[ind(10, 25)] = 0; 
H40[ind(10, 26)] = 0; 
H40[ind(10, 27)] = 0; 
H40[ind(10, 28)] = 0; 
H40[ind(10, 29)] = 0; 
H40[ind(10, 30)] = 0; 
H40[ind(10, 31)] = 0; 
H40[ind(10, 32)] = 0; 
H40[ind(10, 33)] = 0; 
H40[ind(10, 34)] = 0; 
H40[ind(10, 35)] = 0; 
H40[ind(10, 36)] = 0; 
H40[ind(11, 11)] = Ed+W; 
H40[ind(11, 12)] = t; 
H40[ind(11, 13)] = 0; 
H40[ind(11, 14)] = 0; 
H40[ind(11, 15)] = 0; 
H40[ind(11, 16)] = 0; 
H40[ind(11, 17)] = t; 
H40[ind(11, 18)] = 0; 
H40[ind(11, 19)] = 0; 
H40[ind(11, 20)] = 0; 
H40[ind(11, 21)] = 0; 
H40[ind(11, 22)] = 0; 
H40[ind(11, 23)] = V; 
H40[ind(11, 24)] = 0; 
H40[ind(11, 25)] = 0; 
H40[ind(11, 26)] = 0; 
H40[ind(11, 27)] = 0; 
H40[ind(11, 28)] = 0; 
H40[ind(11, 29)] = 0; 
H40[ind(11, 30)] = 0; 
H40[ind(11, 31)] = 0; 
H40[ind(11, 32)] = 0; 
H40[ind(11, 33)] = 0; 
H40[ind(11, 34)] = 0; 
H40[ind(11, 35)] = 0; 
H40[ind(11, 36)] = 0; 
H40[ind(12, 12)] = Ed+W; 
H40[ind(12, 13)] = 0; 
H40[ind(12, 14)] = 0; 
H40[ind(12, 15)] = 0; 
H40[ind(12, 16)] = 0; 
H40[ind(12, 17)] = 0; 
H40[ind(12, 18)] = t; 
H40[ind(12, 19)] = 0; 
H40[ind(12, 20)] = 0; 
H40[ind(12, 21)] = 0; 
H40[ind(12, 22)] = 0; 
H40[ind(12, 23)] = 0; 
H40[ind(12, 24)] = V; 
H40[ind(12, 25)] = 0; 
H40[ind(12, 26)] = 0; 
H40[ind(12, 27)] = 0; 
H40[ind(12, 28)] = 0; 
H40[ind(12, 29)] = 0; 
H40[ind(12, 30)] = 0; 
H40[ind(12, 31)] = 0; 
H40[ind(12, 32)] = 0; 
H40[ind(12, 33)] = 0; 
H40[ind(12, 34)] = 0; 
H40[ind(12, 35)] = 0; 
H40[ind(12, 36)] = 0; 
H40[ind(13, 13)] = 2*Ed+U+W; 
H40[ind(13, 14)] = t; 
H40[ind(13, 15)] = 0; 
H40[ind(13, 16)] = 0; 
H40[ind(13, 17)] = 0; 
H40[ind(13, 18)] = 0; 
H40[ind(13, 19)] = 0; 
H40[ind(13, 20)] = 0; 
H40[ind(13, 21)] = 0; 
H40[ind(13, 22)] = 0; 
H40[ind(13, 23)] = 0; 
H40[ind(13, 24)] = 0; 
H40[ind(13, 25)] = V; 
H40[ind(13, 26)] = 0; 
H40[ind(13, 27)] = 0; 
H40[ind(13, 28)] = 0; 
H40[ind(13, 29)] = 0; 
H40[ind(13, 30)] = 0; 
H40[ind(13, 31)] = 0; 
H40[ind(13, 32)] = 0; 
H40[ind(13, 33)] = 0; 
H40[ind(13, 34)] = 0; 
H40[ind(13, 35)] = 0; 
H40[ind(13, 36)] = 0; 
H40[ind(14, 14)] = 2*Ed+U+W; 
H40[ind(14, 15)] = t; 
H40[ind(14, 16)] = V; 
H40[ind(14, 17)] = 0; 
H40[ind(14, 18)] = 0; 
H40[ind(14, 19)] = 0; 
H40[ind(14, 20)] = 0; 
H40[ind(14, 21)] = 0; 
H40[ind(14, 22)] = 0; 
H40[ind(14, 23)] = 0; 
H40[ind(14, 24)] = 0; 
H40[ind(14, 25)] = 0; 
H40[ind(14, 26)] = V; 
H40[ind(14, 27)] = 0; 
H40[ind(14, 28)] = 0; 
H40[ind(14, 29)] = 0; 
H40[ind(14, 30)] = 0; 
H40[ind(14, 31)] = 0; 
H40[ind(14, 32)] = 0; 
H40[ind(14, 33)] = 0; 
H40[ind(14, 34)] = 0; 
H40[ind(14, 35)] = 0; 
H40[ind(14, 36)] = 0; 
H40[ind(15, 15)] = 2*Ed+U+2*W; 
H40[ind(15, 16)] = 0; 
H40[ind(15, 17)] = V; 
H40[ind(15, 18)] = 0; 
H40[ind(15, 19)] = 0; 
H40[ind(15, 20)] = 0; 
H40[ind(15, 21)] = 0; 
H40[ind(15, 22)] = 0; 
H40[ind(15, 23)] = 0; 
H40[ind(15, 24)] = 0; 
H40[ind(15, 25)] = 0; 
H40[ind(15, 26)] = 0; 
H40[ind(15, 27)] = V; 
H40[ind(15, 28)] = 0; 
H40[ind(15, 29)] = 0; 
H40[ind(15, 30)] = 0; 
H40[ind(15, 31)] = 0; 
H40[ind(15, 32)] = 0; 
H40[ind(15, 33)] = 0; 
H40[ind(15, 34)] = 0; 
H40[ind(15, 35)] = 0; 
H40[ind(15, 36)] = 0; 
H40[ind(16, 16)] = Ed+W; 
H40[ind(16, 17)] = t; 
H40[ind(16, 18)] = 0; 
H40[ind(16, 19)] = 0; 
H40[ind(16, 20)] = 0; 
H40[ind(16, 21)] = 0; 
H40[ind(16, 22)] = 0; 
H40[ind(16, 23)] = 0; 
H40[ind(16, 24)] = 0; 
H40[ind(16, 25)] = 0; 
H40[ind(16, 26)] = 0; 
H40[ind(16, 27)] = 0; 
H40[ind(16, 28)] = V; 
H40[ind(16, 29)] = 0; 
H40[ind(16, 30)] = 0; 
H40[ind(16, 31)] = 0; 
H40[ind(16, 32)] = 0; 
H40[ind(16, 33)] = 0; 
H40[ind(16, 34)] = 0; 
H40[ind(16, 35)] = 0; 
H40[ind(16, 36)] = 0; 
H40[ind(17, 17)] = Ed+2*W; 
H40[ind(17, 18)] = t; 
H40[ind(17, 19)] = 0; 
H40[ind(17, 20)] = 0; 
H40[ind(17, 21)] = 0; 
H40[ind(17, 22)] = 0; 
H40[ind(17, 23)] = 0; 
H40[ind(17, 24)] = 0; 
H40[ind(17, 25)] = 0; 
H40[ind(17, 26)] = 0; 
H40[ind(17, 27)] = 0; 
H40[ind(17, 28)] = 0; 
H40[ind(17, 29)] = V; 
H40[ind(17, 30)] = 0; 
H40[ind(17, 31)] = 0; 
H40[ind(17, 32)] = 0; 
H40[ind(17, 33)] = 0; 
H40[ind(17, 34)] = 0; 
H40[ind(17, 35)] = 0; 
H40[ind(17, 36)] = 0; 
H40[ind(18, 18)] = Ed+2*W; 
H40[ind(18, 19)] = 0; 
H40[ind(18, 20)] = 0; 
H40[ind(18, 21)] = 0; 
H40[ind(18, 22)] = 0; 
H40[ind(18, 23)] = 0; 
H40[ind(18, 24)] = 0; 
H40[ind(18, 25)] = 0; 
H40[ind(18, 26)] = 0; 
H40[ind(18, 27)] = 0; 
H40[ind(18, 28)] = 0; 
H40[ind(18, 29)] = 0; 
H40[ind(18, 30)] = V; 
H40[ind(18, 31)] = 0; 
H40[ind(18, 32)] = 0; 
H40[ind(18, 33)] = 0; 
H40[ind(18, 34)] = 0; 
H40[ind(18, 35)] = 0; 
H40[ind(18, 36)] = 0; 
H40[ind(19, 19)] = Ed; 
H40[ind(19, 20)] = t; 
H40[ind(19, 21)] = 0; 
H40[ind(19, 22)] = 0; 
H40[ind(19, 23)] = 0; 
H40[ind(19, 24)] = 0; 
H40[ind(19, 25)] = t; 
H40[ind(19, 26)] = 0; 
H40[ind(19, 27)] = 0; 
H40[ind(19, 28)] = 0; 
H40[ind(19, 29)] = 0; 
H40[ind(19, 30)] = 0; 
H40[ind(19, 31)] = 0; 
H40[ind(19, 32)] = 0; 
H40[ind(19, 33)] = 0; 
H40[ind(19, 34)] = 0; 
H40[ind(19, 35)] = 0; 
H40[ind(19, 36)] = 0; 
H40[ind(20, 20)] = Ed; 
H40[ind(20, 21)] = t; 
H40[ind(20, 22)] = V; 
H40[ind(20, 23)] = 0; 
H40[ind(20, 24)] = 0; 
H40[ind(20, 25)] = 0; 
H40[ind(20, 26)] = t; 
H40[ind(20, 27)] = 0; 
H40[ind(20, 28)] = 0; 
H40[ind(20, 29)] = 0; 
H40[ind(20, 30)] = 0; 
H40[ind(20, 31)] = 0; 
H40[ind(20, 32)] = 0; 
H40[ind(20, 33)] = 0; 
H40[ind(20, 34)] = 0; 
H40[ind(20, 35)] = 0; 
H40[ind(20, 36)] = 0; 
H40[ind(21, 21)] = Ed+W; 
H40[ind(21, 22)] = 0; 
H40[ind(21, 23)] = V; 
H40[ind(21, 24)] = 0; 
H40[ind(21, 25)] = 0; 
H40[ind(21, 26)] = 0; 
H40[ind(21, 27)] = t; 
H40[ind(21, 28)] = 0; 
H40[ind(21, 29)] = 0; 
H40[ind(21, 30)] = 0; 
H40[ind(21, 31)] = 0; 
H40[ind(21, 32)] = 0; 
H40[ind(21, 33)] = 0; 
H40[ind(21, 34)] = 0; 
H40[ind(21, 35)] = 0; 
H40[ind(21, 36)] = 0; 
H40[ind(22, 22)] = 0; 
H40[ind(22, 23)] = t; 
H40[ind(22, 24)] = 0; 
H40[ind(22, 25)] = 0; 
H40[ind(22, 26)] = 0; 
H40[ind(22, 27)] = 0; 
H40[ind(22, 28)] = t; 
H40[ind(22, 29)] = 0; 
H40[ind(22, 30)] = 0; 
H40[ind(22, 31)] = 0; 
H40[ind(22, 32)] = 0; 
H40[ind(22, 33)] = 0; 
H40[ind(22, 34)] = 0; 
H40[ind(22, 35)] = 0; 
H40[ind(22, 36)] = 0; 
H40[ind(23, 23)] = W; 
H40[ind(23, 24)] = t; 
H40[ind(23, 25)] = 0; 
H40[ind(23, 26)] = 0; 
H40[ind(23, 27)] = 0; 
H40[ind(23, 28)] = 0; 
H40[ind(23, 29)] = t; 
H40[ind(23, 30)] = 0; 
H40[ind(23, 31)] = 0; 
H40[ind(23, 32)] = 0; 
H40[ind(23, 33)] = 0; 
H40[ind(23, 34)] = 0; 
H40[ind(23, 35)] = 0; 
H40[ind(23, 36)] = 0; 
H40[ind(24, 24)] = W; 
H40[ind(24, 25)] = 0; 
H40[ind(24, 26)] = 0; 
H40[ind(24, 27)] = 0; 
H40[ind(24, 28)] = 0; 
H40[ind(24, 29)] = 0; 
H40[ind(24, 30)] = t; 
H40[ind(24, 31)] = 0; 
H40[ind(24, 32)] = 0; 
H40[ind(24, 33)] = 0; 
H40[ind(24, 34)] = 0; 
H40[ind(24, 35)] = 0; 
H40[ind(24, 36)] = 0; 
H40[ind(25, 25)] = Ed+W; 
H40[ind(25, 26)] = t; 
H40[ind(25, 27)] = 0; 
H40[ind(25, 28)] = 0; 
H40[ind(25, 29)] = 0; 
H40[ind(25, 30)] = 0; 
H40[ind(25, 31)] = t; 
H40[ind(25, 32)] = 0; 
H40[ind(25, 33)] = 0; 
H40[ind(25, 34)] = 0; 
H40[ind(25, 35)] = 0; 
H40[ind(25, 36)] = 0; 
H40[ind(26, 26)] = Ed+W; 
H40[ind(26, 27)] = t; 
H40[ind(26, 28)] = V; 
H40[ind(26, 29)] = 0; 
H40[ind(26, 30)] = 0; 
H40[ind(26, 31)] = 0; 
H40[ind(26, 32)] = t; 
H40[ind(26, 33)] = 0; 
H40[ind(26, 34)] = 0; 
H40[ind(26, 35)] = 0; 
H40[ind(26, 36)] = 0; 
H40[ind(27, 27)] = Ed+2*W; 
H40[ind(27, 28)] = 0; 
H40[ind(27, 29)] = V; 
H40[ind(27, 30)] = 0; 
H40[ind(27, 31)] = 0; 
H40[ind(27, 32)] = 0; 
H40[ind(27, 33)] = t; 
H40[ind(27, 34)] = 0; 
H40[ind(27, 35)] = 0; 
H40[ind(27, 36)] = 0; 
H40[ind(28, 28)] = W; 
H40[ind(28, 29)] = t; 
H40[ind(28, 30)] = 0; 
H40[ind(28, 31)] = 0; 
H40[ind(28, 32)] = 0; 
H40[ind(28, 33)] = 0; 
H40[ind(28, 34)] = t; 
H40[ind(28, 35)] = 0; 
H40[ind(28, 36)] = 0; 
H40[ind(29, 29)] = 2*W; 
H40[ind(29, 30)] = t; 
H40[ind(29, 31)] = 0; 
H40[ind(29, 32)] = 0; 
H40[ind(29, 33)] = 0; 
H40[ind(29, 34)] = 0; 
H40[ind(29, 35)] = t; 
H40[ind(29, 36)] = 0; 
H40[ind(30, 30)] = 2*W; 
H40[ind(30, 31)] = 0; 
H40[ind(30, 32)] = 0; 
H40[ind(30, 33)] = 0; 
H40[ind(30, 34)] = 0; 
H40[ind(30, 35)] = 0; 
H40[ind(30, 36)] = t; 
H40[ind(31, 31)] = Ed+W; 
H40[ind(31, 32)] = t; 
H40[ind(31, 33)] = 0; 
H40[ind(31, 34)] = 0; 
H40[ind(31, 35)] = 0; 
H40[ind(31, 36)] = 0; 
H40[ind(32, 32)] = Ed+W; 
H40[ind(32, 33)] = t; 
H40[ind(32, 34)] = V; 
H40[ind(32, 35)] = 0; 
H40[ind(32, 36)] = 0; 
H40[ind(33, 33)] = Ed+2*W; 
H40[ind(33, 34)] = 0; 
H40[ind(33, 35)] = V; 
H40[ind(33, 36)] = 0; 
H40[ind(34, 34)] = W; 
H40[ind(34, 35)] = t; 
H40[ind(34, 36)] = 0; 
H40[ind(35, 35)] = 2*W; 
H40[ind(35, 36)] = t; 
H40[ind(36, 36)] = 2*W; 

H53[ind(1, 1)] = 2*Ed+U+W; 
H53[ind(1, 2)] = V; 
H53[ind(1, 3)] = 0; 
H53[ind(1, 4)] = 0; 
H53[ind(2, 2)] = Ed+W; 
H53[ind(2, 3)] = t; 
H53[ind(2, 4)] = 0; 
H53[ind(3, 3)] = Ed+W; 
H53[ind(3, 4)] = t; 
H53[ind(4, 4)] = Ed+2*W; 

H51[ind(1, 1)] = 2*Ed+U; 
H51[ind(1, 2)] = t; 
H51[ind(1, 3)] = 0; 
H51[ind(1, 4)] = 0; 
H51[ind(1, 5)] = 0; 
H51[ind(1, 6)] = 0; 
H51[ind(1, 7)] = t; 
H51[ind(1, 8)] = 0; 
H51[ind(1, 9)] = 0; 
H51[ind(1, 10)] = 0; 
H51[ind(1, 11)] = 0; 
H51[ind(1, 12)] = 0; 
H51[ind(1, 13)] = 0; 
H51[ind(1, 14)] = 0; 
H51[ind(1, 15)] = 0; 
H51[ind(1, 16)] = 0; 
H51[ind(1, 17)] = 0; 
H51[ind(1, 18)] = 0; 
H51[ind(1, 19)] = 0; 
H51[ind(1, 20)] = 0; 
H51[ind(1, 21)] = 0; 
H51[ind(1, 22)] = 0; 
H51[ind(1, 23)] = 0; 
H51[ind(1, 24)] = 0; 
H51[ind(2, 2)] = 2*Ed+U; 
H51[ind(2, 3)] = t; 
H51[ind(2, 4)] = V; 
H51[ind(2, 5)] = 0; 
H51[ind(2, 6)] = 0; 
H51[ind(2, 7)] = 0; 
H51[ind(2, 8)] = t; 
H51[ind(2, 9)] = 0; 
H51[ind(2, 10)] = 0; 
H51[ind(2, 11)] = 0; 
H51[ind(2, 12)] = 0; 
H51[ind(2, 13)] = 0; 
H51[ind(2, 14)] = 0; 
H51[ind(2, 15)] = 0; 
H51[ind(2, 16)] = 0; 
H51[ind(2, 17)] = 0; 
H51[ind(2, 18)] = 0; 
H51[ind(2, 19)] = 0; 
H51[ind(2, 20)] = 0; 
H51[ind(2, 21)] = 0; 
H51[ind(2, 22)] = 0; 
H51[ind(2, 23)] = 0; 
H51[ind(2, 24)] = 0; 
H51[ind(3, 3)] = 2*Ed+U+W; 
H51[ind(3, 4)] = 0; 
H51[ind(3, 5)] = V; 
H51[ind(3, 6)] = 0; 
H51[ind(3, 7)] = 0; 
H51[ind(3, 8)] = 0; 
H51[ind(3, 9)] = t; 
H51[ind(3, 10)] = 0; 
H51[ind(3, 11)] = 0; 
H51[ind(3, 12)] = 0; 
H51[ind(3, 13)] = 0; 
H51[ind(3, 14)] = 0; 
H51[ind(3, 15)] = 0; 
H51[ind(3, 16)] = 0; 
H51[ind(3, 17)] = 0; 
H51[ind(3, 18)] = 0; 
H51[ind(3, 19)] = 0; 
H51[ind(3, 20)] = 0; 
H51[ind(3, 21)] = 0; 
H51[ind(3, 22)] = 0; 
H51[ind(3, 23)] = 0; 
H51[ind(3, 24)] = 0; 
H51[ind(4, 4)] = Ed; 
H51[ind(4, 5)] = t; 
H51[ind(4, 6)] = 0; 
H51[ind(4, 7)] = 0; 
H51[ind(4, 8)] = 0; 
H51[ind(4, 9)] = 0; 
H51[ind(4, 10)] = t; 
H51[ind(4, 11)] = 0; 
H51[ind(4, 12)] = 0; 
H51[ind(4, 13)] = 0; 
H51[ind(4, 14)] = 0; 
H51[ind(4, 15)] = 0; 
H51[ind(4, 16)] = 0; 
H51[ind(4, 17)] = 0; 
H51[ind(4, 18)] = 0; 
H51[ind(4, 19)] = 0; 
H51[ind(4, 20)] = 0; 
H51[ind(4, 21)] = 0; 
H51[ind(4, 22)] = 0; 
H51[ind(4, 23)] = 0; 
H51[ind(4, 24)] = 0; 
H51[ind(5, 5)] = Ed+W; 
H51[ind(5, 6)] = t; 
H51[ind(5, 7)] = 0; 
H51[ind(5, 8)] = 0; 
H51[ind(5, 9)] = 0; 
H51[ind(5, 10)] = 0; 
H51[ind(5, 11)] = t; 
H51[ind(5, 12)] = 0; 
H51[ind(5, 13)] = 0; 
H51[ind(5, 14)] = 0; 
H51[ind(5, 15)] = 0; 
H51[ind(5, 16)] = 0; 
H51[ind(5, 17)] = 0; 
H51[ind(5, 18)] = 0; 
H51[ind(5, 19)] = 0; 
H51[ind(5, 20)] = 0; 
H51[ind(5, 21)] = 0; 
H51[ind(5, 22)] = 0; 
H51[ind(5, 23)] = 0; 
H51[ind(5, 24)] = 0; 
H51[ind(6, 6)] = Ed+W; 
H51[ind(6, 7)] = 0; 
H51[ind(6, 8)] = 0; 
H51[ind(6, 9)] = 0; 
H51[ind(6, 10)] = 0; 
H51[ind(6, 11)] = 0; 
H51[ind(6, 12)] = t; 
H51[ind(6, 13)] = 0; 
H51[ind(6, 14)] = 0; 
H51[ind(6, 15)] = 0; 
H51[ind(6, 16)] = 0; 
H51[ind(6, 17)] = 0; 
H51[ind(6, 18)] = 0; 
H51[ind(6, 19)] = 0; 
H51[ind(6, 20)] = 0; 
H51[ind(6, 21)] = 0; 
H51[ind(6, 22)] = 0; 
H51[ind(6, 23)] = 0; 
H51[ind(6, 24)] = 0; 
H51[ind(7, 7)] = 2*Ed+U+W; 
H51[ind(7, 8)] = t; 
H51[ind(7, 9)] = 0; 
H51[ind(7, 10)] = 0; 
H51[ind(7, 11)] = 0; 
H51[ind(7, 12)] = 0; 
H51[ind(7, 13)] = t; 
H51[ind(7, 14)] = 0; 
H51[ind(7, 15)] = 0; 
H51[ind(7, 16)] = 0; 
H51[ind(7, 17)] = 0; 
H51[ind(7, 18)] = 0; 
H51[ind(7, 19)] = 0; 
H51[ind(7, 20)] = 0; 
H51[ind(7, 21)] = 0; 
H51[ind(7, 22)] = 0; 
H51[ind(7, 23)] = 0; 
H51[ind(7, 24)] = 0; 
H51[ind(8, 8)] = 2*Ed+U+W; 
H51[ind(8, 9)] = t; 
H51[ind(8, 10)] = V; 
H51[ind(8, 11)] = 0; 
H51[ind(8, 12)] = 0; 
H51[ind(8, 13)] = 0; 
H51[ind(8, 14)] = t; 
H51[ind(8, 15)] = 0; 
H51[ind(8, 16)] = 0; 
H51[ind(8, 17)] = 0; 
H51[ind(8, 18)] = 0; 
H51[ind(8, 19)] = 0; 
H51[ind(8, 20)] = 0; 
H51[ind(8, 21)] = 0; 
H51[ind(8, 22)] = 0; 
H51[ind(8, 23)] = 0; 
H51[ind(8, 24)] = 0; 
H51[ind(9, 9)] = 2*Ed+U+2*W; 
H51[ind(9, 10)] = 0; 
H51[ind(9, 11)] = V; 
H51[ind(9, 12)] = 0; 
H51[ind(9, 13)] = 0; 
H51[ind(9, 14)] = 0; 
H51[ind(9, 15)] = t; 
H51[ind(9, 16)] = 0; 
H51[ind(9, 17)] = 0; 
H51[ind(9, 18)] = 0; 
H51[ind(9, 19)] = 0; 
H51[ind(9, 20)] = 0; 
H51[ind(9, 21)] = 0; 
H51[ind(9, 22)] = 0; 
H51[ind(9, 23)] = 0; 
H51[ind(9, 24)] = 0; 
H51[ind(10, 10)] = Ed+W; 
H51[ind(10, 11)] = t; 
H51[ind(10, 12)] = 0; 
H51[ind(10, 13)] = 0; 
H51[ind(10, 14)] = 0; 
H51[ind(10, 15)] = 0; 
H51[ind(10, 16)] = t; 
H51[ind(10, 17)] = 0; 
H51[ind(10, 18)] = 0; 
H51[ind(10, 19)] = 0; 
H51[ind(10, 20)] = 0; 
H51[ind(10, 21)] = 0; 
H51[ind(10, 22)] = 0; 
H51[ind(10, 23)] = 0; 
H51[ind(10, 24)] = 0; 
H51[ind(11, 11)] = Ed+2*W; 
H51[ind(11, 12)] = t; 
H51[ind(11, 13)] = 0; 
H51[ind(11, 14)] = 0; 
H51[ind(11, 15)] = 0; 
H51[ind(11, 16)] = 0; 
H51[ind(11, 17)] = t; 
H51[ind(11, 18)] = 0; 
H51[ind(11, 19)] = 0; 
H51[ind(11, 20)] = 0; 
H51[ind(11, 21)] = 0; 
H51[ind(11, 22)] = 0; 
H51[ind(11, 23)] = 0; 
H51[ind(11, 24)] = 0; 
H51[ind(12, 12)] = Ed+2*W; 
H51[ind(12, 13)] = 0; 
H51[ind(12, 14)] = 0; 
H51[ind(12, 15)] = 0; 
H51[ind(12, 16)] = 0; 
H51[ind(12, 17)] = 0; 
H51[ind(12, 18)] = t; 
H51[ind(12, 19)] = 0; 
H51[ind(12, 20)] = 0; 
H51[ind(12, 21)] = 0; 
H51[ind(12, 22)] = 0; 
H51[ind(12, 23)] = 0; 
H51[ind(12, 24)] = 0; 
H51[ind(13, 13)] = 2*Ed+U+W; 
H51[ind(13, 14)] = t; 
H51[ind(13, 15)] = 0; 
H51[ind(13, 16)] = 0; 
H51[ind(13, 17)] = 0; 
H51[ind(13, 18)] = 0; 
H51[ind(13, 19)] = V; 
H51[ind(13, 20)] = 0; 
H51[ind(13, 21)] = 0; 
H51[ind(13, 22)] = 0; 
H51[ind(13, 23)] = 0; 
H51[ind(13, 24)] = 0; 
H51[ind(14, 14)] = 2*Ed+U+W; 
H51[ind(14, 15)] = t; 
H51[ind(14, 16)] = V; 
H51[ind(14, 17)] = 0; 
H51[ind(14, 18)] = 0; 
H51[ind(14, 19)] = 0; 
H51[ind(14, 20)] = V; 
H51[ind(14, 21)] = 0; 
H51[ind(14, 22)] = 0; 
H51[ind(14, 23)] = 0; 
H51[ind(14, 24)] = 0; 
H51[ind(15, 15)] = 2*Ed+U+2*W; 
H51[ind(15, 16)] = 0; 
H51[ind(15, 17)] = V; 
H51[ind(15, 18)] = 0; 
H51[ind(15, 19)] = 0; 
H51[ind(15, 20)] = 0; 
H51[ind(15, 21)] = V; 
H51[ind(15, 22)] = 0; 
H51[ind(15, 23)] = 0; 
H51[ind(15, 24)] = 0; 
H51[ind(16, 16)] = Ed+W; 
H51[ind(16, 17)] = t; 
H51[ind(16, 18)] = 0; 
H51[ind(16, 19)] = 0; 
H51[ind(16, 20)] = 0; 
H51[ind(16, 21)] = 0; 
H51[ind(16, 22)] = V; 
H51[ind(16, 23)] = 0; 
H51[ind(16, 24)] = 0; 
H51[ind(17, 17)] = Ed+2*W; 
H51[ind(17, 18)] = t; 
H51[ind(17, 19)] = 0; 
H51[ind(17, 20)] = 0; 
H51[ind(17, 21)] = 0; 
H51[ind(17, 22)] = 0; 
H51[ind(17, 23)] = V; 
H51[ind(17, 24)] = 0; 
H51[ind(18, 18)] = Ed+2*W; 
H51[ind(18, 19)] = 0; 
H51[ind(18, 20)] = 0; 
H51[ind(18, 21)] = 0; 
H51[ind(18, 22)] = 0; 
H51[ind(18, 23)] = 0; 
H51[ind(18, 24)] = V; 
H51[ind(19, 19)] = Ed+W; 
H51[ind(19, 20)] = t; 
H51[ind(19, 21)] = 0; 
H51[ind(19, 22)] = 0; 
H51[ind(19, 23)] = 0; 
H51[ind(19, 24)] = 0; 
H51[ind(20, 20)] = Ed+W; 
H51[ind(20, 21)] = t; 
H51[ind(20, 22)] = V; 
H51[ind(20, 23)] = 0; 
H51[ind(20, 24)] = 0; 
H51[ind(21, 21)] = Ed+2*W; 
H51[ind(21, 22)] = 0; 
H51[ind(21, 23)] = V; 
H51[ind(21, 24)] = 0; 
H51[ind(22, 22)] = W; 
H51[ind(22, 23)] = t; 
H51[ind(22, 24)] = 0; 
H51[ind(23, 23)] = 2*W; 
H51[ind(23, 24)] = t; 
H51[ind(24, 24)] = 2*W; 

H62[ind(1, 1)] = 2*Ed+U+W; 
H62[ind(1, 2)] = t; 
H62[ind(1, 3)] = 0; 
H62[ind(1, 4)] = 0; 
H62[ind(1, 5)] = 0; 
H62[ind(1, 6)] = 0; 
H62[ind(2, 2)] = 2*Ed+U+W; 
H62[ind(2, 3)] = t; 
H62[ind(2, 4)] = V; 
H62[ind(2, 5)] = 0; 
H62[ind(2, 6)] = 0; 
H62[ind(3, 3)] = 2*Ed+U+2*W; 
H62[ind(3, 4)] = 0; 
H62[ind(3, 5)] = V; 
H62[ind(3, 6)] = 0; 
H62[ind(4, 4)] = Ed+W; 
H62[ind(4, 5)] = t; 
H62[ind(4, 6)] = 0; 
H62[ind(5, 5)] = Ed+2*W; 
H62[ind(5, 6)] = t; 
H62[ind(6, 6)] = Ed+2*W; 

H60[ind(1, 1)] = 2*Ed+U; 
H60[ind(1, 2)] = t; 
H60[ind(1, 3)] = 0; 
H60[ind(1, 4)] = 0; 
H60[ind(1, 5)] = t; 
H60[ind(1, 6)] = 0; 
H60[ind(1, 7)] = 0; 
H60[ind(1, 8)] = 0; 
H60[ind(1, 9)] = 0; 
H60[ind(1, 10)] = 0; 
H60[ind(1, 11)] = 0; 
H60[ind(1, 12)] = 0; 
H60[ind(1, 13)] = 0; 
H60[ind(1, 14)] = 0; 
H60[ind(1, 15)] = 0; 
H60[ind(1, 16)] = 0; 
H60[ind(2, 2)] = 2*Ed+U+W; 
H60[ind(2, 3)] = t; 
H60[ind(2, 4)] = 0; 
H60[ind(2, 5)] = 0; 
H60[ind(2, 6)] = t; 
H60[ind(2, 7)] = 0; 
H60[ind(2, 8)] = 0; 
H60[ind(2, 9)] = 0; 
H60[ind(2, 10)] = 0; 
H60[ind(2, 11)] = 0; 
H60[ind(2, 12)] = 0; 
H60[ind(2, 13)] = 0; 
H60[ind(2, 14)] = 0; 
H60[ind(2, 15)] = 0; 
H60[ind(2, 16)] = 0; 
H60[ind(3, 3)] = 2*Ed+U+W; 
H60[ind(3, 4)] = V; 
H60[ind(3, 5)] = 0; 
H60[ind(3, 6)] = 0; 
H60[ind(3, 7)] = t; 
H60[ind(3, 8)] = 0; 
H60[ind(3, 9)] = 0; 
H60[ind(3, 10)] = 0; 
H60[ind(3, 11)] = 0; 
H60[ind(3, 12)] = 0; 
H60[ind(3, 13)] = 0; 
H60[ind(3, 14)] = 0; 
H60[ind(3, 15)] = 0; 
H60[ind(3, 16)] = 0; 
H60[ind(4, 4)] = Ed+W; 
H60[ind(4, 5)] = 0; 
H60[ind(4, 6)] = 0; 
H60[ind(4, 7)] = 0; 
H60[ind(4, 8)] = t; 
H60[ind(4, 9)] = 0; 
H60[ind(4, 10)] = 0; 
H60[ind(4, 11)] = 0; 
H60[ind(4, 12)] = 0; 
H60[ind(4, 13)] = 0; 
H60[ind(4, 14)] = 0; 
H60[ind(4, 15)] = 0; 
H60[ind(4, 16)] = 0; 
H60[ind(5, 5)] = 2*Ed+U+W; 
H60[ind(5, 6)] = t; 
H60[ind(5, 7)] = 0; 
H60[ind(5, 8)] = 0; 
H60[ind(5, 9)] = t; 
H60[ind(5, 10)] = 0; 
H60[ind(5, 11)] = 0; 
H60[ind(5, 12)] = 0; 
H60[ind(5, 13)] = 0; 
H60[ind(5, 14)] = 0; 
H60[ind(5, 15)] = 0; 
H60[ind(5, 16)] = 0; 
H60[ind(6, 6)] = 2*Ed+U+2*W; 
H60[ind(6, 7)] = t; 
H60[ind(6, 8)] = 0; 
H60[ind(6, 9)] = 0; 
H60[ind(6, 10)] = t; 
H60[ind(6, 11)] = 0; 
H60[ind(6, 12)] = 0; 
H60[ind(6, 13)] = 0; 
H60[ind(6, 14)] = 0; 
H60[ind(6, 15)] = 0; 
H60[ind(6, 16)] = 0; 
H60[ind(7, 7)] = 2*Ed+U+2*W; 
H60[ind(7, 8)] = V; 
H60[ind(7, 9)] = 0; 
H60[ind(7, 10)] = 0; 
H60[ind(7, 11)] = t; 
H60[ind(7, 12)] = 0; 
H60[ind(7, 13)] = 0; 
H60[ind(7, 14)] = 0; 
H60[ind(7, 15)] = 0; 
H60[ind(7, 16)] = 0; 
H60[ind(8, 8)] = Ed+2*W; 
H60[ind(8, 9)] = 0; 
H60[ind(8, 10)] = 0; 
H60[ind(8, 11)] = 0; 
H60[ind(8, 12)] = t; 
H60[ind(8, 13)] = 0; 
H60[ind(8, 14)] = 0; 
H60[ind(8, 15)] = 0; 
H60[ind(8, 16)] = 0; 
H60[ind(9, 9)] = 2*Ed+U+W; 
H60[ind(9, 10)] = t; 
H60[ind(9, 11)] = 0; 
H60[ind(9, 12)] = 0; 
H60[ind(9, 13)] = V; 
H60[ind(9, 14)] = 0; 
H60[ind(9, 15)] = 0; 
H60[ind(9, 16)] = 0; 
H60[ind(10, 10)] = 2*Ed+U+2*W; 
H60[ind(10, 11)] = t; 
H60[ind(10, 12)] = 0; 
H60[ind(10, 13)] = 0; 
H60[ind(10, 14)] = V; 
H60[ind(10, 15)] = 0; 
H60[ind(10, 16)] = 0; 
H60[ind(11, 11)] = 2*Ed+U+2*W; 
H60[ind(11, 12)] = V; 
H60[ind(11, 13)] = 0; 
H60[ind(11, 14)] = 0; 
H60[ind(11, 15)] = V; 
H60[ind(11, 16)] = 0; 
H60[ind(12, 12)] = Ed+2*W; 
H60[ind(12, 13)] = 0; 
H60[ind(12, 14)] = 0; 
H60[ind(12, 15)] = 0; 
H60[ind(12, 16)] = V; 
H60[ind(13, 13)] = Ed+W; 
H60[ind(13, 14)] = t; 
H60[ind(13, 15)] = 0; 
H60[ind(13, 16)] = 0; 
H60[ind(14, 14)] = Ed+2*W; 
H60[ind(14, 15)] = t; 
H60[ind(14, 16)] = 0; 
H60[ind(15, 15)] = Ed+2*W; 
H60[ind(15, 16)] = V; 
H60[ind(16, 16)] = 2*W; 

H71[ind(1, 1)] = 2*Ed+U+W; 
H71[ind(1, 2)] = t; 
H71[ind(1, 3)] = 0; 
H71[ind(1, 4)] = 0; 
H71[ind(2, 2)] = 2*Ed+U+2*W; 
H71[ind(2, 3)] = t; 
H71[ind(2, 4)] = 0; 
H71[ind(3, 3)] = 2*Ed+U+2*W; 
H71[ind(3, 4)] = V; 
H71[ind(4, 4)] = Ed+2*W;

H80[ind(1, 1)] = 2*Ed+U+2*W;

  givens(dim11, dim11, dim11, H11, aval11, avet11);
  givens(dim22, dim22, dim22, H22, aval22, avet22);
  givens(dim20, dim20, dim20, H20, aval20, avet20);
  givens(dim33, dim33, dim33, H33, aval33, avet33);
  givens(dim31, dim31, dim31, H31, aval31, avet31);
  givens(dim44, dim44, dim44, H44, aval44, avet44);
  givens(dim42, dim42, dim42, H42, aval42, avet42);
  givens(dim40, dim40, dim40, H40, aval40, avet40);
  givens(dim53, dim53, dim53, H53, aval53, avet53);
  givens(dim51, dim51, dim51, H51, aval51, avet51);
  givens(dim62, dim62, dim62, H62, aval62, avet62);
  givens(dim60, dim60, dim60, H60, aval60, avet60);
  givens(dim71, dim71, dim71, H71, aval71, avet71);
  givens(dim80, dim80, dim80, H80, aval80, avet80);
  
  avaltotal[1] = 0;
  spinz[1] = 0;
  
  int inttotal = 2;
  for(int i = 1; i <= dim11; i = i+1)
  {
    avaltotal[inttotal] = aval11[i];
    avaltotal[inttotal+1] = aval11[i];
    spinz[inttotal] = 1./2.;
    spinz[inttotal+1] = -1./2.;
    inttotal = inttotal + 2;
  }
  
  for(int i = 1; i <= dim22; i = i+1)
  {
    avaltotal[inttotal] = aval22[i];
    avaltotal[inttotal+1] = aval22[i];
    spinz[inttotal] = 2./2.;
    spinz[inttotal+1] = -2./2.;
    inttotal = inttotal + 2;
  }
  
  for(int i = 1; i <= dim33; i = i+1)
  {
    avaltotal[inttotal] = aval33[i];
    avaltotal[inttotal+1] = aval33[i];
    spinz[inttotal] = 3./2.;
    spinz[inttotal+1] = -3./2.;
    inttotal = inttotal + 2;
  }
  
  for(int i = 1; i <= dim31; i = i+1)
  {
    avaltotal[inttotal] = aval31[i];
    avaltotal[inttotal+1] = aval31[i];
    spinz[inttotal] = 1./2.;
    spinz[inttotal+1] = -1./2.;
    inttotal = inttotal + 2;
  }
  
  for(int i = 1; i <= dim44; i = i+1)
  {
    avaltotal[inttotal] = aval44[i];
    avaltotal[inttotal+1] = aval44[i];
    spinz[inttotal] = 4./2.;
    spinz[inttotal+1] = -4./2.;
    inttotal = inttotal + 2;
  }
  
  for(int i = 1; i <= dim42; i = i+1)
  {
    avaltotal[inttotal] = aval42[i];
    avaltotal[inttotal+1] = aval42[i];
    spinz[inttotal] = 2./2.;
    spinz[inttotal+1] = -2./2.;
    inttotal = inttotal + 2;
  }
  
  for(int i = 1; i <= dim53; i = i+1)
  {
    avaltotal[inttotal] = aval53[i];
    avaltotal[inttotal+1] = aval53[i];
    spinz[inttotal] = 3./2.;
    spinz[inttotal+1] = -3./2.;
    inttotal = inttotal + 2;
  }
  
  for(int i = 1; i <= dim51; i = i+1)
  {
    avaltotal[inttotal] = aval51[i];
    avaltotal[inttotal+1] = aval51[i];
    spinz[inttotal] = 1./2.;
    spinz[inttotal+1] = -1./2.;
    inttotal = inttotal + 2;
  }
  
  for(int i = 1; i <= dim62; i = i+1)
  {
    avaltotal[inttotal] = aval62[i];
    avaltotal[inttotal+1] = aval62[i];
    spinz[inttotal] = 2./2.;
    spinz[inttotal+1] = -2./2.;
    inttotal = inttotal + 2;
  }
  
  for(int i = 1; i <= dim71; i = i+1)
  {
    avaltotal[inttotal] = aval71[i];
    avaltotal[inttotal+1] = aval71[i];
    spinz[inttotal] = 1./2.;
    spinz[inttotal+1] = -1./2.;
    inttotal = inttotal + 2;
  }
  
  for(int i = 1; i <= dim20; i = i+1)
  {
    avaltotal[inttotal] = aval20[i];
    spinz[inttotal] = 0;
    inttotal = inttotal + 1;
  }
  
  for(int i = 1; i <= dim40; i = i+1)
  {
    avaltotal[inttotal] = aval40[i];
    spinz[inttotal] = 0;
    inttotal = inttotal + 1;
  }
  
  for(int i = 1; i <= dim60; i = i+1)
  {
    avaltotal[inttotal] = aval60[i];
    spinz[inttotal] = 0;
    inttotal = inttotal + 1;
  }
  
  for(int i = 1; i <= dim80; i = i+1)
  {
    avaltotal[inttotal] = aval80[i];
    spinz[inttotal] = 0;
    inttotal = inttotal + 1;
  }
  
  
  const int MAXCHAR = 100;
  char nome[MAXCHAR];
  
    sprintf(nome, "mag_susceptimp%.2f.plt", W);
  
  ofstream mag_suscept(nome);
  
  for(double temp = 0.006; temp <= 5; temp = temp + 0.001)
  {
  double magnetic_suscept = media_spin2(spinz, avaltotal, dimtotal, temp);
  
  mag_suscept << temp << "\t" << magnetic_suscept << endl;
  }
  
  }
  //Diagonalizamos o Hamiltoniano com a impureza. Agora devemos diagonalizar o Hamiltoniano sem a impureza.
  
}