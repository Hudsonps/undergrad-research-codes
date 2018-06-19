#include <cmath>
#include <iostream>
#include <cstdlib>
#include <limits.h>
#include <float.h>
#include <fstream>
#include <iomanip>
#include "giv.h"

const double D = 1.0; //largura da banda
const double LAMB = 10; //parâmetro de discretização logarítmica
const double ED = 0.002; //energia associada ao sítio do ponto quântico
const int DIM = 40; //dimensão da matriz que será diagonalizada

using namespace std;

int ind(int i, int j)
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

double sinal(double x)
{
  return x/sqrt(x*x);   //Essa função retorna o sinal de x.
}

double Em(double m, double signal)
  {
    return sinal(signal)*D*pow(LAMB, -m);
    //Nesta função, são considerados os casos de energia positivas E de negativas
  }

double V(double E, double alfa)
{
  //Nessa função, aparece o termo de acoplamento dos sítios de fio e ponto, sem considerar ainda a discretização logarítmica.
  double gamma = 0.1+0.1*sin(M_PI*E);
  return (gamma/M_PI);
}

double Vred(double m, double signal, double alfa)
{
  //Aqui está o termo de acoplamento da versão logarítmica, seguindo a nossa receita.
  double erg = Em(m, signal);
  return sqrt(log(LAMB)*V(erg, alfa)*V(erg, alfa)*sqrt(erg*erg));
}

double delta(double ergaval, double z)
{
  //Em nossa proposta, os coeficientes de Anderson da diagonal são da forma D Lambda ^ - (m + z).
  //Isso significa que um autovalor deve ser reescrito na forma D Lambda ^ - ( m + z + gamma).
  //Podemos, portanto, descobrir em qual intervalo (m, m+1) o autovalor está e, a partir disso, determinar o gamma correspondente.
  //A convenção usual, contudo, é escolhr gamma variando entre -0.5 a 0.5, por isso é melhor trabalhar com intervalos da forma (m - 0.5, m + 0.5).
  
  //Módulo de um autovalor pode ser escrito na forma D Lambda ^ -(s + z)
  //Então s é dado por:
  
  double s;
  s = -log(sqrt(ergaval*ergaval)/(D))/log(LAMB) - z;
  
  //s = m + gamma, com m ainda indeterminado e gamma também indeterminado
  // Como m é sempre um inteiro, basta descobrir qual é o objeto fracionário correspondente.
  //Por exemplo, se s = 5.4, então m = 5 e gamma = 0.4.
  //Já se s = 6.8, então m = 7 e gamma = -0.2
  //Em geral, devemos olhar se a parte fracionária de um número é maior ou menor que 0.5. 
  
  double m;
  
  if( (s - int(s)) >= 0.5) //Se a condição neste if é verdadeira, estamos mais próximos do inteiro superior, e gamma é negativo.
  {
    m = int(s) + 1;
  }
  else //Do contrário, então estamos mais próximos do inteiro inferior.
  {
    m = int(s);
  }
  
  //Determinado m, gamma está determinado. gamma = s - m.
  double gamma;
  gamma = s - m;
  
  //gamma e delta, na malha logarítmica, estão relacionados por:  delta = pi gamma * sinal(ergaval)
  
  double deltinha;
  deltinha = sinal(ergaval)*gamma*M_PI;
  
  return deltinha;
}

int main(int argc, char* argv[])
{
  double alfa = 0.0;
  switch(argc){
  case 2:
    alfa = atof(argv[1]);
    break;

  default:
    cerr << "***No. argumentos incorreto: " << argv[0] << " alfa" << endl;
    exit(1);
  }
  const int MAXCHAR = 100;
  char nome[MAXCHAR];
  int ialf = int(1000. * alfa);
  cout << ialf << endl;
  sprintf(nome, "deltas.%04d.plt", ialf);
  ofstream delts(nome);
  
  ofstream avetd("avetd.plt");
  
  ofstream spectrad("spectrad.plt");
  
  ofstream spectraf0("spectraf0.plt");
  

  int dim = DIM % 2? DIM : DIM + 1;
  
  //O termo da ponta da direita da matriz é da forma H[dim][dim]. Na forma de um vetor, isso corresponde a H[ind(dim, dim)].
  
  double *ham = new double[ind(dim, dim) + 1]; //Essa matriz guardará os elementos do Hamiltoniano.
  double *aval = new double[int(dim + 1)]; //Essa matriz guardará os autovalores do Hamiltoniano.

  double *avalaux = new double [int(dim + 1)];
  //Essa matriz será necessária para calcular derivada dos autovalores com relação a z.
  //Ela guardará os autovalores provenientes de uma interação com z levemente menor que outro, para que se possa derivar.
  
  double *vm = new double [int(dim + 1)]; //Ema certa etapa do programa, precisaremos dos V_m. eles ficarão guardados nessa matriz
  //Poderíamos acessá-los da matriz Ham, mas após ela ser diagonalizada, os coeficientes são perdidos.
  
  double *em = new double [int(dim + 1)];
  //Do mesmo modo, também precisaremos dos coeficientes da diagonal presentes antes da diagonalização.

  double **avet = new double*[dim+1];
  for(int lin = 1; lin <=dim; lin++){
    avet[lin] = new double[dim+1];
  }
 
  int i = 0;
  int j = 0; //índices a serem usados para percorrer os vários elementos de matrizes
  
  for(double w = -1; w <= -0; w += 0.001){
    //esse w corresponde ao z convencional. Mas como queremos derivada de z de algumas grandezas, precisamos, para cada w, pegar 2 z próximos.
    for(double z = w; z < (w + 1e-4 + 1e-8); z += 1e-4){ //a razão de termos dois z tão próximos é que precisamos da derivada de z dos autovalores.

//Note que só faz sentido tomar a derivada quando estamos na segunda 'rodada' deste for, pois precisamos de, pelo menos, dois pontos.
//O critério para saber se estamos na segunda rodada é se z é ou não igual a w. Na segunda rodada, z é diferente de w por 1e-4.
//Então sabemos se estamos no segundo turno quando z - w > 0.5e-4.

//Em síntese:   z - w  < 0.5e-4 ->  PRimeiro turno.               z - w > 0.5e-4 ->  Segundo turno.
    
//Esse for produz a matriz nula.
    for(i = 1; i <= ind(dim, dim); i++)
      {
	ham[i] = 0.0;  //Com esse for, anulamos todos os elementos da matriz.
      }
  
    ham[ind(1,1)] = ED; //Na primeira posição, fica a energia do dot.
  
    //O Hamiltoniano possui elementos na diagonal da forma Em- e Em+;  e possui elementos nas primeiras linha e coluna da forma Vm- e Vm+.
    //Para manter as energias em ordem crescente, começamos com energias negativas.
    //Como escolhemos a dimensão como ímpar, há um número par disponível de posições para esses coeficientes.
    //Metade para energias negativas. Metade para energias positivas.
    // (dim-1)/2 para cada. O índice m, portanto, irá variar de 0 a (dim-1)/2 - 1
    //Por exemplo: se escolhemos dim = 5, dim-1 = 4.  (dim-1)/2 = 2. Então temos m = 0 e m = 1. O m máximo é 2 - 1.
    //Como m tem um valor máximo, vamos definir uma variável para explicitar isso, para deixar o programa mais claro.
  
    double maxm = double((dim - 1)/2 - 1); //máximo valor para o índice m.
  
    //Começamos pelas energias negativas. Neste regime, incremento de m leva a incremento na energia, por isso começamos do menor m.
    //Há termos da forma ham[i,i] e ham[i,1] ou ham[1, i]  não nulos.
  
    i = 2; //i guarda a primeira linha em que aparecem os elementos negativos; nesse contexto, é a segunda.
    for(double m = 0; m <= maxm; m++)
      {
	ham[ind(i, i)] = Em(m+z, -1);
	ham[ind(i, 1)] = Vred(m+z, -1, alfa)/sqrt(2.0*D);
	i++; //para mudar de linha.
      }
  
    //Isso deu conta dos termos de energia negativa. Se tudo foi feito corretamente, i sofreu maxm + 1 incrementos (1 para cada m possível).
    //O +1 é porque também consideramos m = 0, totalizando maxm + 1 números.
    
    //O valor inicial de i era 2. Então i está em (dim -1)/2 + 2.
    //No exemplo de dim = 5, isso é 4.
  
    for(double m = maxm; m >= 0; m--)
      {
	ham[ind(i, i)] = Em(m+z, 1);
	ham[ind(i, 1)] = Vred(m+z, 1, alfa)/sqrt(2.0*D); //isso é Vm(Em+)/sqrt(2D)
	i++; //para mudar de linha.
      }
  
    //Se i estava em 2 + maxm + 1, agora recebeu outro incremento de maxm + 1, totalizando 4 + 2 maxm = 4 - 2 + dim - 1 =  dim + 1
    //Isso é esperado. Se tudo foi feito corretamente, o índice i deve ser tal que estamos na primeira linha/coluna fora da matriz
  
    // cout << i - dim << endl;
 
 
 
    //Para checar se a matriz foi corretamente escrita, podemos fazer com que o programa a escreva na tela.
/*     for(i = 1; i <= dim; i++)
       {
 	for(j = 1; j <= dim; j++)
 	  {
 	    cout << ham[ind(i, j)] << " " ;
 	  }
 	cout << endl;
       }
*/

//Antes de diagonalizar a matriz, vamos guardar os coeficientes V_m, que serão perdidos após a diagonalização.
  for(int m = 2; m <= dim; m++)
  {
    vm[m] = sqrt(2.*D)*ham[ind(m, 1)]; //esses coeficientes foram previamente divididos por sqrt(2D), mas queremos eles na forma convencional.
    em[m] = ham[ind(m, m)];
  }
  
    
    givens(dim, dim, dim, ham, aval, avet);

    for(int ell = 1; ell <= dim; ell++){
      double del = delta(aval[ell], z);
      double delTeor = atan(M_PI * pow(V(aval[ell], alfa), 2.0)/(2.0*D) / (ED - aval[ell]));

      double Y1 = M_PI*V(aval[ell], alfa)*V(aval[ell], alfa)/tan(del)/2 - ED;
      double Y2 = M_PI*V(aval[ell], alfa)*V(aval[ell], alfa)/tan(delTeor)/2 - ED;
      double Y3 = Y1 + aval[ell];

    if(abs(aval[ell]) > 1e-7)
      {
      delts << setw(11) << aval[ell] << "\t" << Y1 << "\t" << Y2 << "\t" << Y3 << endl;
      
      double dotcomponent = avet[ell][1];
      
     
      avetd << aval[ell] << "\t" << dotcomponent << endl;
      }

      }

//Antes de encerrarmos este for, temos que guardar os autovalores em uma matriz. Quando a próxima rodada começar, a matriz aval será perdida.
//Para isso, usaremos o vetor avalaux, sobre o qual os autovalores serão guardados.
//Além disso, para o cálculo da densidade espectral de c_d, precisaremos guardar a componente do dot.

   if((z - w) < 0.5e-4)  //Esse if diz que estamos no primeiro turno.
   {
    for(int ell = 1; ell <= dim; ell++)
    {
      avalaux[ell] = aval[ell]; //Vamos guardar os autovalores do primeiro turno para, posteriormente, poder calcular a derivada dos aval em relação a z.
    }
   }
   
//Já se estivermos na segunda rodada, é hora de calcular a densidade espectral
  if((z - w) > 0.5e-4)
  {
    for(int ell = 1; ell <= dim; ell++)
    {
      double diffze = (z - w)/(aval[ell] - avalaux[ell]); //derivada de z com relação a varepsilon
      
      //cout << diffze << endl;
      
      double spectraldot = avet[ell][1]*avet[ell][1]*diffze*sinal(-aval[ell]); //densidade espectral de c_d, de acordo com o relatório 5.
      
      double delTeor = atan(M_PI * pow(V(aval[ell], alfa), 2.0)/(2.0*D) / (ED - aval[ell]));  //fórmula teórica da densidade espectral de c_d
      
      double spectraldteor = 2.0*D*sin(delTeor)*sin(delTeor)/(M_PI*M_PI*V(aval[ell], alfa)*V(aval[ell], alfa));
      
	if(abs(aval[ell]) > 1e-7)
	{
      spectrad << aval[ell] << "\t" << spectraldot << "\t" << spectraldteor << "\t" << spectraldot/spectraldteor << endl;
      	}
      //Podemos, ainda neste for, calcular também a densidade espectral do operador f_0
      
      double spectralf0 = 0;
      for(int m = 2; m <= dim; m++)
      {
	spectralf0 = spectralf0 + avet[ell][m]*vm[m]/sqrt(2.0*D);
      }
      spectralf0 = spectralf0*spectralf0;
      spectralf0 = spectralf0*diffze*sinal(-aval[ell]);
      
		
      double spectralf0teor = V(aval[ell], alfa)*V(aval[ell], alfa)*cos(delTeor)*cos(delTeor)/(2.*D);
	
      if(abs(aval[ell]) > 1e-7)
	{
      spectraf0 << aval[ell] << "\t" << spectralf0 << "\t" << spectralf0teor << "\t" << spectralf0/spectralf0teor << endl;
	}
     //spectraf0 << aval[ell] << "\t" << (avet[ell][2]/avet[ell][1])*(avet[ell][2]/avet[ell][1])*D << "\t" << vm[2]*vm[2]/((aval[ell] - em[2])*(aval[ell] - em[2])) << endl ;
    }
  }
   
}
}
cout << vm[2]*vm[2] << endl;
}
      

