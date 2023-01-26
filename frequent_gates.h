//classe che raccoglie gates di uso frequente, in file separato per essere maneggiato senza rischi in caso di necessità di aggiungere gates alla mappa. 
#include<iostream>
#include<cstdlib>
#include<fstream>
#include<cmath>
#include<vector>
#include<string>
#include<complex>
#include<algorithm>
#include<random>
#include<ctime>
#include<typeinfo>
#include<set>
#include<map>
using namespace std;

//importo una matrice quadrata a valori double o complex (in formato "(real,imag)" ) da un matrix_file esterno 
vector<vector<complex<double>>> import_complex_matrix(string matrix_file){
	ifstream INPUT(matrix_file);
	vector<complex<double>> prel_vec;
	complex<double> z;
	while(INPUT>>z){
		prel_vec.push_back(z);
	}
	int N=sqrt(prel_vec.size());
	vector<vector<complex<double>>> M;
	int k=0;	
	for(int i=0; i<N; i++){
		vector<complex<double>> riga;
		for (int j=0; j<N; j++){
			riga.push_back(prel_vec[k]);
			k++;	
		}
		M.push_back(riga);
	}
	return M;
}


class frequent_gates {
	public: 
		map<string, vector<vector<complex<double>>>> gate{
			{	"NOT"	,	{{0,1},{1, 0}}															},
			{	"CNOT"	,	{{1,0,0,0},{0,1,0,0},{0,0,0,1},{0,0,1,0}} 								},
			{	"H"		,	{{1/sqrt(2), 1/sqrt(2)},{1/sqrt(2), -1/sqrt(2)}} 						},	
			{	"SX"	, 	{{{0.5, 0.5},{0.5, -0.5}},{{0.5, -0.5},{0.5,0.5}}}						},
			{	"SY"	, 	{{{0.5, 0.5},{-0.5, -0.5}},{{0.5, 0.5},{0.5,0.5}}}						},
			{	"SW"	, 	{{{0.5, 0.5},{0, -0.5*sqrt(2)}},{{0.5*sqrt(2), 0},{0.5,0.5}}}			},
			{	"iSWAP"	, 	{{1,0,0,0},{0,0,{0,1},0},{0,{0,1},0,0},{0,0,0,1}}						},
			{	"CZ"	, 	{{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,-1}}								},
			{	"ID1"	,	{{1,0},{0,1}}															},
			{	"GOO"	, 	{{1,0,0,0},{0,0,{0,-1},0},{0,{0,-1},0,0},{0,0,0,{0.5*sqrt(3),-0.5}}}	},					
		};
		
		//funzione per aggiungere da file esterno un gate alla classe dato il nome del gate e del file
		void add_gate(string gate_name, string file_name){
			ifstream INPUT(file_name);
			gate[gate_name]=import_complex_matrix(file_name);
		}
		//funzione per aggiungere un gate alla classe dato il nome del gate e la matrice corrispondente
		void add_gate(string gate_name, vector<vector<complex<double>>> G){
			gate[gate_name]=G;
		}
		
};
