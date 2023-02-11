
		
//Prima versione libreria tesi infoQ
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
#include<set>
#include<armadillo>
#include "frequent_gates.h"
using namespace std;
using namespace arma;

//funzione che restituisce un vettore con i valori degli n qubit dato il label dello stato di base computazionale e il numero di qubits (non fa nient'altro che convertire in binario insomma...)
vector<int> bits(int m, int n){
	vector<int> bits_;
	int i=n-1;
	while (i>=0){
		bits_.push_back(m/int(pow(2,i)));
		m%=int(pow(2,i));
		i--;
	}
	return (bits_);
}

//funzione che dato un vettore con i valori degli n qubit restituisce il label dello stato di base computazionale (non fa nient'altro che convertire dal binario insomma...)
int comp(vector<int> bits){
	int sum=0;
	for(int i=0; i<bits.size(); i++){
		sum+=bits[i]*pow(2,bits.size()-1-i);
	}
	return (sum);
}

//funzione che dato un elemento di un vettore restituisce (il minor) intero posizione corrispondente
template <typename mytype>
int posizione(mytype el, vector<mytype> v){
	int i=0;
	for (auto c : v){
		if(c==el) { break; }
		else i++;
	}
	return i;
}

//funzione che restituisce la norma quadra di un ket
double psi_norm (vector<complex<double>> psi){
	double sum=0.;
	for(int i=0; i<psi.size(); i++){
		sum+=norm(psi[i]);
	}
	return (sum);
}

//funzione che normalizza un histo di double dato come reference e restituisce la norma come return
double distr_norm(vector<double> &histo){
	double sum=0;
	for(auto el : histo){
		sum+=el;
	}
	for(auto el : histo){
		el/=sum;
	}
	return sum;
}

//funzione che presi due vettori li riordina secondo la stessa permutazione, tale che il primo sia in ordine decrescente
class my_duple{
	public:
		double elA;
		double elB;
};
bool my_duplecompare(const my_duple &a, const my_duple &b)
{
    return (a.elA < b.elA);
}
void my_duplesort(vector<double> &vA, vector<double> &vB){
	if(vA.size()!=vB.size()) cout<<"dimensioni vettori in my_duplesort incompatibili!"<<endl;
	vector<my_duple> duplevec;
	for(int i=0; i<vA.size(); i++){
		my_duple Dup;
		Dup.elA=vA[i];
		Dup.elB=vB[i];
		duplevec.push_back(Dup);
	}
	sort(duplevec.rbegin(), duplevec.rend(), my_duplecompare);
	for(int i=0; i<duplevec.size(); i++){
		vA[i]=duplevec[i].elA;
		vB[i]=duplevec[i].elB;
	}
}
//funzione che, date due frazioni (double) e una reference a un istogramma, rimuove dal dominio dell'istogramma una porzione iniziale e finale corrispondenti alle due frazioni ed eventualmente (argomento opzionale "rinorm") rinormalizza il contenuto
void trim_histo(vector<double>& histo,string mode="", double leftfrac=1./5, double rightfrac=2./5){
	int lenhisto=histo.size();
	int leftlim=floor(lenhisto*leftfrac);
	int rightlim=floor(lenhisto*rightfrac);
	histo.erase(histo.begin(),histo.begin()+leftlim);
	histo.erase(histo.end()-rightlim, histo.end());
	if(mode=="rinorm"){
		distr_norm(histo);	
	}	
}

//distribuzione esponenziale, normalizzata
double exp_distro(int i, int N){
	double integral=(exp(1)-exp(1-N))/(exp(1)-1);
	return(exp(-i)/integral);
}
//Porter-Thomas distribution, normalizzata
double PT_distro(int i, int N){
	double integral=exp(N)*(1-exp(-N*N))*N/(exp(N)-1);
	return(N*exp(-N*i)/integral);
}

//funzione che restiruisce la somma element-wise pesata di due vettori
template <typename mytype>
vector<mytype> weight_sum(complex<double> w1, vector<mytype> v1, complex<double> w2, vector<mytype> v2){
	if(v1.size()!=v2.size()) cout<<"dimensioni vettori incompatibili!"<<endl;
	vector<mytype> v3;
	for(int i=0; i<v1.size(); i++){
		v3.push_back(w1*v1[i]+w2*v2[i]);
	}
	return(v3);
}
//funzione che date due matrici ne fa la somma
template <typename mytypeA, typename mytypeB>
vector<vector<complex<double>>> mat_sum(vector<vector<mytypeA>> A, vector<vector<mytypeB>> B ){
	vector<vector<complex<double>>> C;
	if(A.size()!=B.size() || A[0].size()!=B[0].size()) {
		cout<<"dimensioni matrici da sommare sono incompatibili"<<endl;
	}else{
		for(int i=0; i< A.size(); i++){
			vector<complex<double>> riga;
			for(int j=0; j<A[0].size(); j++){
				riga.push_back((complex<double>)A[i][j]+(complex<double>)B[i][j]);
			}
			C.push_back(riga);
		}
	}
	return C;
}

//funzione che applica matrice a vettore restituendo vettore di complessi (DA CONTROLLARE)
template<typename mytypeM, typename mytypeV>
vector<complex<double>> apply_mat(vector<vector<mytypeM>> M, vector<mytypeV> V){
	if(M[0].size()!=V.size()) cout<<"dimensioni matriec e vettore incompatibili!"<<endl;
	vector<complex<double>> result;
	for(int i=0; i<V.size(); i++){
		complex<double> sum(0, 0);
		for(int j=0; j<V.size();j++){
			sum+=M[i][j]*V[j];
		}
		result.push_back(sum);
	}
	return(result);
}

//funzione che date due matrice ne fa il prodotto consueto, ed eventualmente col mode== "Badjoint" calcola A * B trasposto coniugato
template <typename mytypeA, typename mytypeB>
vector<vector<complex<double>>> mat_prod(vector<vector<mytypeA>> A, vector<vector<mytypeB>> B , string mode="normal"){
	vector<vector<complex<double>>>C;
	if(mode=="normal"){
		if(A[0].size()!=B.size()) cout<<"dimensioni interne matrici differenti!"<<endl;
		for(int ai=0; ai< A.size(); ai++){
			vector<complex<double>> riga;
			for(int bj=0; bj<B[0].size();bj++){
				complex<double> sum=0;
				for(int bi=0; bi<B.size(); bi++){
					sum+=(complex<double>)A[ai][bi]*(complex<double>)B[bi][bj];
				}
				riga.push_back(sum);
			}
			C.push_back(riga);
		}
	}else if(mode=="Badjoint"){
		if(A[0].size()!=B[0].size()) cout<<"dimensioni interne matrici differenti!"<<endl;
		for(int ai=0; ai< A.size(); ai++){
			vector<complex<double>> riga;
			for(int bi=0; bi<B.size(); bi++){
				complex<double> sum=0;
				for(int bj=0; bj<B[0].size(); bj++){
					sum+=(complex<double>)A[ai][bj]*conj((complex<double>)B[bi][bj]);
				}
				riga.push_back(sum);
			}
			C.push_back(riga);
		}
	}
	return C;	
}

//funzione che date due matrici ne fa il prodotto di kronecker
template <typename mytypeA, typename mytypeB>
vector<vector<complex<double>>> kron_prod(vector<vector<mytypeA>> A, vector<vector<mytypeB>> B ){
	//inizializzo a matrice di 0 con opportune dimensioni
	vector<vector<complex<double>>> C(A.size()*B.size(), vector<complex<double>>(A[0].size()*B[0].size(), 0.));
	for(int ai=0; ai< A.size(); ai++){
		for(int aj=0; aj< A[0].size(); aj++){
			for(int bi=0; bi< B.size(); bi++){
				for(int bj=0; bj< B[0].size(); bj++){
					C[ai*B.size()+bi][ aj*B[0].size()+bj]=(complex<double>)A[ai][aj]*(complex<double>)B[bi][bj];
				}
			}
		}	
	}
	return C;
}

//funzione che genera una matrice identità NxN dato N
vector<vector<int>> ID(int N){
	vector<vector<int>> myID;
	for(int i=0; i<N; i++){
		vector<int> riga;
		for(int j=0; j<N; j++){
			if(i==j){
				riga.push_back(1);
			}else{
				riga.push_back(0);
			}			
		}
		myID.push_back(riga);
	}
	return myID;
}

//funzione che genera una matrice SWAP tra i due qubit pos1 e pos2 dato il numero di qubit 
vector<vector<int>> SWAP(int pos1, int pos2, int n){
	vector<vector<int>> mySWAP;
	for(int i=0; i<pow(2, n); i++){
		vector<int> riga;
		vector<int> swapped_bits=bits(i, n);
		int a =swapped_bits[pos1];
		swapped_bits[pos1]=swapped_bits[pos2];
		swapped_bits[pos2]=a;
		for(int l=0; l<pow(2, n);l++){
			if(l==comp(swapped_bits)) {
				riga.push_back(1);
			}else{
				riga.push_back(0);
			}
		}
		mySWAP.push_back(riga);
	}
	return mySWAP;
}

//funzione che genera una matrice SWAP qubits dato un vettore che mi rappresenta la permutazione (in ogni posizione la posizione immagine)
vector<vector<int>> SWAP(vector<int> perm){
	int n=perm.size(); 
	int N=pow(2, n);
	vector<vector<int>> mySWAP;
	for(int i=0; i<N; i++){
		vector<int> riga;
		vector<int> swapped_bits=bits(i, n);
		int temp;
		for(int k=0; k<perm.size(); k++){
			if(perm[k]>k){
				temp=swapped_bits[k];
				swapped_bits[k]=swapped_bits[perm[k]];
				swapped_bits[perm[k]]=temp;
			}	
		}
		for(int l=0; l<N;l++){
			if(l==comp(swapped_bits)) {
				riga.push_back(1);
			}else{
				riga.push_back(0);
			}
		}
		mySWAP.push_back(riga);
	}
	return mySWAP;
}

//funzione che stampa un complesso in forma real+i*imag
void print_complex (complex<double> z){
	cout<<real(z);
	if(imag(z)!=0) {
		if(imag(z)>0) {
			cout<<"+"<<"i"<<imag(z);
		}else{
			cout<<"i"<<imag(z);
		}
	}else{
		cout<<" ";
	}
}

//funzione che stampa su file matrice di complessi, i files sono prodotti in automatico concatenando title con il numero p_c che aumenta ad ogni stampa
//eventualmente l'argomento opzionale mode== "easy_read" moltiplica i valori per 10N e li arrotonda all'intero superiore (in modulo), così da occupare poco spazio (cioè "(##,##)") ma dare un 'idea della matrice, inoltre impagina in modo piu largo/leggibile
template<typename el_type>
void print_mat(vector<vector<el_type>> A, int &p_c, string title="",  string mode="normal"){ 
		string output_file=title+"_"+to_string(p_c)+".txt";
		ofstream OUTPUT(output_file);
		if (mode=="normal"){											
			for(vector<el_type> v : A){											
				for (auto z : v){
					OUTPUT<<z<<" ";
				}																
				OUTPUT<<endl;												
			}															
		}else if(mode=="easy_read"){
			for(auto v : A){
				for (auto x : v){																								
					complex<double> z=complex<double>(x)*complex<double>(10*A.size());																																
					z={real(z)>=0.?ceil(real(z)):floor(real(z)) , imag(z)>=0.?ceil(imag(z)):floor(imag(z))};
					OUTPUT<<z<<"	";
				}
				OUTPUT<<endl<<endl;;
			}
		}
		p_c++;
}

//funzione che stampa su un unico ostream dato in reference matrice di complessi, 
//eventualmente l'argomento opzionale mode== "easy_read" moltiplica i valori per 10 e li arrotonda all'intero superiore (in modulo), così da occupare poco spazio (cioè "(##,##)") ma dare un 'idea della matrice, inoltre impagina in modo piu largo/leggibile
template<typename el_type>
void print_mat(vector<vector<el_type>> A, ostream & OUTPUT,  string title="", string mode="normal"){ 
	OUTPUT<<title<<endl;
	if (mode=="normal"){											
		for(vector<el_type> v : A){											
			for (auto z : v){
				OUTPUT<<z<<" ";
			}																
			OUTPUT<<endl;												
		}															
	}else if(mode=="easy_read"){
		for(auto v : A){
			for (auto x : v){																								
				complex<double> z=complex<double>(x)*complex<double>(10*A.size());																																
				z={real(z)>=0.?ceil(real(z)):floor(real(z)) , imag(z)>=0.?ceil(imag(z)):floor(imag(z))};
				OUTPUT<<z<<"	";
			}
			OUTPUT<<endl<<endl;;
		}
	}
}

//funzione che stampa a terminale vettore, con arg opzionale la stringa del file .txt dove stamparlo
template<typename el_type>
void print_vec(vector<el_type> v, string wheretxt="terminal"){
	if(wheretxt=="terminal"){
		for (auto z : v){
		cout<<z<<" ";
		}
		cout<<endl;
	}else{
		ofstream OUTPUT(wheretxt);
		for (auto z : v){
		OUTPUT<<z<<" ";
		}
		OUTPUT<<endl;
	}	
}


//funzione che stampa su un ostream dato come reference un vettore
template<typename el_type>
void print_vec(vector<el_type> v, ostream &OUTPUT, string title=""){
	OUTPUT<<title<<endl;
	for (auto z : v){
	OUTPUT<<z<<" ";
	}
	OUTPUT<<endl;
}

// importo una matrice quadrata a valori int da un matrix_file esterno 
vector<vector<int>> import_int_matrix(string matrix_file){
	ifstream INPUT(matrix_file);
	vector<int> prel_vec;
	int z;
	while(INPUT>>z){
		prel_vec.push_back(z);
	}
	int N=sqrt(prel_vec.size());
	vector<vector<int>> M;
	int k=0;	
	for(int i=0; i<N; i++){
		vector<int> riga;
		for (int j=0; j<N; j++){
			riga.push_back(prel_vec[k]);
			k++;	
		}
		M.push_back(riga);
	}
	return M;
}

//funzione che data una reference a una variabile random int, e un valore estremo superiore escluso, se il parametro opzionale mode è posto a "norep" trasforma la variabile in un numero random diverso dal precedente, tranne nel caso in cui eccezione (parametro opzionale di default 1) sia uguale a 0, che rappresenta un(a funzione di) qualche iteratore così che il primo random dell'iterazione non abbia giustamente il vincolo della non ripetizione nemmeno con mode="no_rep" 
void global_random(int &my_random, int modulo, string mode="normal", int eccezione=-1){
	int my_temp=rand()%modulo;
	if(mode=="norep"){
		while(my_temp==my_random && eccezione!=0){
			my_temp=rand()%modulo;
		}
	}
	my_random=my_temp;
}

//classe matrice densità con vari metodi
class density_matrix {
	
	private:
		
		//numero qubits
		int n;
		
		//dimensione spazio di hilbert
		int N;
		
		//array matrice densità
		vector<vector<complex<double>>> M;
			
		//QBITS_LEFT: vettore di interi in (1, n) che tiene memoria delle "identità" dei qubits man mano che alcuni vengono tracciati via
		vector<int> qbits_left;
		
	public:
		
		//costruttore: prende come input il # qubit e la modalità
		density_matrix (int n_=3, string mode="random", vector<vector<complex<double>>> matrix={{}}){
			n=n_;
			N=pow(2, n);
			//inizializzo qbits_left
			for (int i=1; i<=n; i++){
				qbits_left.push_back(i);
			}
			
			//INIZIALIZZAZIONE HA DIVERSE MODALITA': random, random_pure, allbitszero, unisup, manual ecc
			
			//ro inizializzata randomicamente alla piu generica matrice densità (anche mista) tramite lapack (armadillo). Ro random, sommo l'aggiunta, sommo identità per minimo avl, divido per traccia
			if(mode=="random"){
				Mat<complex<double>> R(N, N);
				for(int i=0; i<N;i++){
					for(int j=0; j<N; j++){
						complex<double> z(rand(), rand());
						R(i,j)=z;
					}
				}
				R+=trans(R);
				vec eigvalR;
				eig_sym(eigvalR, R);
				Mat<complex<double>> L(N, N, fill::zeros);
				for(int i=0; i<N;i++){
					L(i,i)=eigvalR[0];
				}
				R+=L;
				R/=arma::trace(R);
				//trasformo la matrice ottenuta da formato "armadillo" a quello utilizzato nel resto del codice
				for(int i=0; i<N;i++){
					vector<complex<double>> riga;
					for(int j=0; j<N; j++){
						riga.push_back(R(i,j));
					}
					M.push_back(riga);
				}	
			}
			
			//ro = stato puro |psiXpsi|, generato randomicamente
			else if (mode=="random_pure"){
				//creo psi casualmente
				vector<complex<double>> psi;
				for (int i=0; i<N; i++){
					complex<double> z(rand(),rand());
					psi.push_back(z);
				}
				double sqnorm=psi_norm(psi);
				//riempio M con |psiXpsi| mentre normalizzo QUESTO è UN PASSAGGIO RIDONDANTE, L'INFO è GIA IN PSI MA è COMODO LAVORARE CON L'ARRAY, DA OTTIMIZZARE
				for (int i=0; i<N; i++){
					vector<complex<double>> riga;
					for(int j=0; j<N; j++){
						riga.push_back(psi[i]*conj(psi[j])/sqnorm);
					}
					M.push_back(riga);
				}
			}
			
			//ro inizializzata allo stato zero di tutti i qubit cioè |m>=1;
			else if (mode=="allbitszero"){
				vector<vector<complex<double>>> array(N, vector<complex<double>>(N, 0));
				array[0][0]=1;
				M=array; //QUESTA è UN OPERAZIONE DI COPIA? CONSUMA RISORSE RISPARMIABILI?
			}
			
			//ro inizializzata alla sovrapposizione uniforme di tutti gli stati della base computazionale TODO
			else if (mode=="unisup"){
				vector<vector<complex<double>>> array(N, vector<complex<double>>(N, 1./N));
				M=array; //QUESTA è UN OPERAZIONE DI COPIA? CONSUMA RISORSE RISPARMIABILI?
			}
			
			//ro inizilizzata manualmente a una certa matrice.
			else if (mode=="manual"){
				M=matrix;
			}
		}
		
		//get n
		int get_n (){
			return n;
		}
		
		//get N
		int get_N (){
			return N;
		}
		
		//get M
		vector<vector<complex<double>>> get_M(){
			return M;
		}
		
		//get qbits_left
		vector<int> get_qbits_left(){
			return qbits_left;
		}
		
		//change M
		void set_M(vector<vector<complex<double>>> my_M){
			M=my_M;
		}
		
		//change Mij
		void set_M(int i, int j, complex<double> val){
			M[i][j]=val;
		}
		
		
		//valore traccia (check che sia uno ecc), la tengo complessa per scovare eventuali errori
		complex<double> trace(){                               
			complex<double> sum(0.0,0.0);
			for (int i=0; i<N; i++){
				sum+=M[i][i];
			}													
			return sum;
		}
		
		//traccia parziale sul qubit k, k intero in (1,n) che indentifica univocamente un dato qubit "fisico"
		void trace_over_qubit(int k){						
			int pos=posizione(k, qbits_left);				
			vector<vector<complex<double>>> newM;				
			for(int i=0; i<N/2;i++){							
				vector<complex<double>> riga;
				for(int j=0; j<N/2;j++){
					vector<int> bits_i=bits(i, n-1);
					vector<int> bits_j=bits(j, n-1);
					complex<double> m0, m1;						
					//trovo la posizione della cifra binaria a cui corrisponde il qubit k e calcolo gli el di matrice m0 ed m1 per i valori 0 e 1 di quella cifra binaria. 
					bits_i.insert(bits_i.begin()+pos, 0);
					bits_j.insert(bits_j.begin()+pos, 0);
					m0=M[comp(bits_i)][comp(bits_j)];			
					bits_i[pos]=1;
					bits_j[pos]=1;
					m1=M[comp(bits_i)][comp(bits_j)];
					//l'elemento ij della nuova matrice è la somma su a=0,1 di ro_ia^ja, dove "ia" va letta come concatenazione binaria per ritrovare l'|m> corrispondente della base computazionale dello spazio di partenza (grande il doppio di quello tracciato)
					riga.push_back(m0+m1);
				}
				newM.push_back(riga);							
			}
			M=newM;												
			n-=1;											
			N/=2;												
			qbits_left.erase(qbits_left.begin()+pos);													
		}
		
		//traccia parziale sul set s di qubits k, k', k''... interi in (1,n) che indentificano univocamente dei qubit "fisici"
		void trace_over_qubit(set<int> s){						
			//numero qubit che voglio tracciare via 
			int o=s.size();												
			vector<int> poss;
			//riempio vettore con le posizioni nell'espansione binaria relativi ai qubit su cui voglio tracciare, in ordine crescente
			for(auto k : s){
				poss.push_back(posizione(k, qbits_left));
			}
			sort(poss.begin(), poss.end());								
			//calcolo										
			vector<vector<complex<double>>> newM;				
			for(int i=0; i<N/pow(2,o);i++){							
				vector<complex<double>> riga;
				for(int j=0; j<N/pow(2,o);j++){
					vector<int> bits_i=bits(i, n-o);							
					vector<int> bits_j=bits(j, n-o);										
					//calcolo l'el come somma su tutte le combinazioni binarie dei qubits degli el della matrice originaria con indici ottenuti per concatenazione
					complex<double> sum(0.,0.);
					//metto degli zeri(temporanei) nelle posizioni relative agli o qubits, inserendoli dalla posizione piu a destra sino alla piu a sinistra, tenendo conto che nel computo della posizione finale ci sono gli spazi occupati dai precedenti inserimenti 
					for(int r=o-1; r>=0; r--){
						bits_i.insert(bits_i.begin()+poss[r]-r, 0);
						bits_j.insert(bits_j.begin()+poss[r]-r, 0);			
					}															
					//comincio a sostituire nelle posizioni appena create tutte le possibili combinazioni binarie 
					sum+=M[comp(bits_i)][comp(bits_j)];
					for(int count=1; count<pow(2,o); count++){
						vector<int> comb=bits(count, o);						
						for(int r=0; r<o; r++){
							bits_i[poss[r]]=comb[r];
							bits_j[poss[r]]=comb[r];													
						}
						sum+=M[comp(bits_i)][comp(bits_j)];							
					}																						
					//l'elemento ij della nuova matrice è la somma su a=0,1 di ro_ia^ja, dove "ia" va letta come concatenazione binaria per ritrovare l'|m> corrispondente della base computazionale dello spazio di partenza (grande il doppio di quello tracciato)
					riga.push_back(sum);
				}
				newM.push_back(riga);							
			}
			M=newM;												
			n-=1;											
			N/=2;
			for(int r=o-1; r>=0; r--){
				qbits_left.erase(qbits_left.begin()+poss[r]);
			}												
																
		}		
		
		//proiezione su un esito in (0,1) di un qubit k. Di default l'outcome dipende dalle probabilità nella ro, ma l'esito può essere forzato.
		int project(int k, int outcome=2){
			complex<double> sum0(0,0);
			int pos=posizione(k, qbits_left);
			//calcolo la traccia della matrice proiettata su 0 per la probabilità e per ri-normalizzare successivamente
			for(int i=0; i<N; i++){ 																		
				if(bits(i, n)[pos]==0){
					sum0+=M[i][i];																																			
				}
			}
			//se outcome è un valore non in {0,1} e dunque non lecito da inserire, oppure se non viene inserito un outcome, questo è generato randomicamente come segue
			if(outcome!=0&&outcome!=1){
				//estraggo l'outcome con una certa probabilità: ma essa è veramente la traccia della matrice proiettata? Anche se la misura è solo in un sottospazio?
				double q_trial=(double)rand()/RAND_MAX;
				if(q_trial<real(sum0)) {
					outcome=0;
				}else{
					outcome=1;
				}
			}																								
																											
			//uccide gli elementi di matrice relativi al sottospazio in cui il qubit misurato ha un outcome diverso da quello imposto
			for(int i=0; i<N; i++){
				for(int j=0; j<N; j++){
					if(bits(i, n)[pos]!=outcome||bits(j, n)[pos]!=outcome){
						M[i][j]=complex<double>(0.0, 0.0);
					}else {
						complex<double> tra=(outcome==0)?sum0:complex<double>(1,0)-sum0; 
						M[i][j]/=tra; //normalizzazione
					}
				}
			}
			return outcome;
		}	
};

//classe ket con vari metodi (DA CONTROLLARE)
class ket {
	
	private:
		
		//numero qubits
		int n;
		
		//dimensione spazio di hilbert
		int N;
		
		//vettore ket
		vector<complex<double>> V;
			
		//QBITS_LEFT: vettore di interi in (1, n) che tiene memoria delle "identità" dei qubits man mano che alcuni vengono tracciati via
		vector<int> qbits_left;
		
	public:
		
		//costruttore: prende come input il # qubit e la modalità
		ket (int n_=3, string mode="random", vector<complex<double>> column={}){
			n=n_;
			N=pow(2, n);
			//inizializzo qbits_left
			for (int i=1; i<=n; i++){
				qbits_left.push_back(i);
			}
			
			//INIZIALIZZAZIONE HA DIVERSE MODALITA': random, allbitszero, unisup, manual ecc
			
			
			//stato puro |psi>, generato randomicamente
			if (mode=="random_pure"){
				//creo psi casualmente
				vector<complex<double>> psi;
				for (int i=0; i<N; i++){
					complex<double> z(rand(),rand());
					psi.push_back(z);
				}
				double sqnorm=psi_norm(psi);
				//riempio M con |psiXpsi| mentre normalizzo QUESTO è UN PASSAGGIO RIDONDANTE, L'INFO è GIA IN PSI MA è COMODO LAVORARE CON L'ARRAY, DA OTTIMIZZARE
				for (int i=0; i<N; i++){
					psi[i]/=sqrt(sqnorm);
				}
			}
			
			//ro inizializzata allo stato zero di tutti i qubit cioè |m>=1;
			else if (mode=="allbitszero"){
				vector<complex<double>> v(N, 0);
				v[0]=1;
				V=v; //QUESTA è UN OPERAZIONE DI COPIA? CONSUMA RISORSE RISPARMIABILI?
			}
			
			//ro inizializzata alla sovrapposizione uniforme di tutti gli stati della base computazionale TODO
			else if (mode=="unisup"){
				vector<complex<double>> v(N, 1./N);
				V=v; //QUESTA è UN OPERAZIONE DI COPIA? CONSUMA RISORSE RISPARMIABILI?
			}
			
			//ro inizilizzata manualmente a una certa matrice.
			else if (mode=="manual"){
				V=column;
			}
		}
		
		//get n
		int get_n (){
			return n;
		}
		
		//get N
		int get_N (){
			return N;
		}
		
		//get M
		vector<complex<double>> get_V(){
			return V;
		}
		
		//get qbits_left
		vector<int> get_qbits_left(){
			return qbits_left;
		}
		
		//change M
		void set_V(vector<complex<double>> my_V){
			V=my_V;
		}
		
		//change Mij
		void set_V(int i, complex<double> val){
			V[i]=val;
		}
		
		//set_qubits_left
		void set_qbits_left(vector<int> my_qbits_left){
			qbits_left=my_qbits_left;
		}
		
		//pensata per funzionare con un ket K  (DA CONTROLLARE)
		int project(int k, int outcome=2){
			double sum0=0;
			int pos=posizione(k, qbits_left);
			//calcolo la traccia della matrice proiettata su 0 per la probabilità e per ri-normalizzare successivamente
			for(int i=0; i<N; i++){ 																		
				if(bits(i, n)[pos]==0){
					sum0+=norm(V[i]);																																			
				}
			}
			//se outcome è un valore non in {0,1} e dunque non lecito da inserire, oppure se non viene inserito un outcome, questo è generato randomicamente come segue
			if(outcome!=0&&outcome!=1){
				//estraggo l'outcome con una certa probabilità: ma essa è veramente la traccia della matrice proiettata? Anche se la misura è solo in un sottospazio?
				double q_trial=(double)rand()/RAND_MAX;
				if(q_trial<real(sum0)) {
					outcome=0;
				}else{
					outcome=1;
				}
			}																								
																											
			//uccide gli elementi di ket relativi al sottospazio in cui il qubit misurato ha un outcome diverso da quello imposto
			for(int i=0; i<N; i++){
				if(bits(i, n)[pos]!=outcome){
					V[i]=complex<double>(0.0, 0.0);
				}else {
					complex<double> sqnorm=(outcome==0)?sum0:complex<double>(1,0)-sum0; 
					V[i]/=sqrt(sqnorm); //normalizzazione
				}
			}
			return outcome;
		}
};


//data una density_matrix by-reference, applico gate di singolo qubit sul qubit k alla matrice densità
void single_qbit_gate(density_matrix &Ro, vector<vector<complex<double>>> gate, int k){
	vector<vector<complex<double>>>my_id_check={{1, 0},{0,1}};
	if(gate!=my_id_check){
		int pos=posizione(k, Ro.get_qbits_left());
		/*
		vector<vector<complex<double>> Id{{1, 0},{0, 1}}; 
		vector<vector<complex<double>> U={{1}};
		for(int g=0; g<qbits_left.size(); g++){
			if(g!=pos){
				U=kron_prod(Id, U);
			} else {
				U=kron_prod(gate, U);
			}
		}
		*/
		//costruisco con kronecker tra Id di dimensioni opportune e gate la matrice U corrispondente
		vector<vector<complex<double>>> U=kron_prod(ID(pow(2, pos)), kron_prod(gate, ID(pow(2, Ro.get_n()-pos-1)))); 	
		//calcolo M'=UMU+ 
		Ro.set_M(mat_prod(U, mat_prod(Ro.get_M(), U, "Badjoint")));
	}
	return;
}

//dato un ket by-reference, applico gate di singolo qubit sul qubit k al ket  (DA CONTROLLARE)
void single_qbit_gate(ket &K, vector<vector<complex<double>>> gate, int k){
	vector<vector<complex<double>>>my_id_check={{1, 0},{0,1}};
	if(gate!=my_id_check){
		int pos=posizione(k, K.get_qbits_left());
		/*
		vector<vector<complex<double>> Id{{1, 0},{0, 1}}; 
		vector<vector<complex<double>> U={{1}};
		for(int g=0; g<qbits_left.size(); g++){
			if(g!=pos){
				U=kron_prod(Id, U);
			} else {
				U=kron_prod(gate, U);
			}
		}
		*/
		//costruisco con kronecker tra Id di dimensioni opportune e gate la matrice U corrispondente
		vector<vector<complex<double>>> U=kron_prod(ID(pow(2, pos)), kron_prod(gate, ID(pow(2, K.get_n()-pos-1)))); 	
		//calcolo K'=UK 
		K.set_V(apply_mat(U, K.get_V()));
	}
	return;
}

//data una density_matrix by-reference, gate di due qubit (da performare tramite riordinamento in modo da avere i due qubit in causa consecutivi, in posizione 1 e 2 , e poi riportare all'ordine iniziale)
void double_qbit_gate(density_matrix &Ro, vector<vector<complex<double>>> gate, int k1, int k2){
	int pos1=posizione(k1, Ro.get_qbits_left());
	int pos2=posizione(k2, Ro.get_qbits_left());
	//calcolo U con l'opportuna trasformazione di Swap in modo che agisca su M pur essendo costruita come gate su qubit 0 e 1
	//lo swap è una composizione di 2 swap, il primo manda il primo qubit in pos 0, il secondo manda il secondo qubit in pos 1
	vector<vector<complex<double>>> mySWAP=mat_prod(SWAP(pos1, 0, Ro.get_n()), SWAP(pos2, 1, Ro.get_n()));								
	//U= SWAP Uswapped SWAP+, con U prod_kron di gate di due qubit e identità
	vector<vector<complex<double>>> U=mat_prod(mySWAP, mat_prod(kron_prod(gate, ID(pow(2, Ro.get_n()-2))), mySWAP, "Badjoint"));
	//calcolo la nuova ro M'=UMU+
	Ro.set_M(mat_prod(U, mat_prod(Ro.get_M(), U, "Badjoint")));
	return;
}	

// (DA CONTROLLARE) data un ket by-reference, gate di due qubit (da performare tramite riordinamento in modo da avere i due qubit in causa consecutivi, in posizione 1 e 2 , e poi riportare all'ordine iniziale)
void double_qbit_gate(ket &K, vector<vector<complex<double>>> gate, int k1, int k2){
	int pos1=posizione(k1, K.get_qbits_left());
	int pos2=posizione(k2, K.get_qbits_left());
	//calcolo U con l'opportuna trasformazione di Swap in modo che agisca su M pur essendo costruita come gate su qubit 0 e 1
	//lo swap è una composizione di 2 swap, il primo manda il primo qubit in pos 0, il secondo manda il secondo qubit in pos 1
	vector<vector<complex<double>>> mySWAP=mat_prod(SWAP(pos1, 0, K.get_n()), SWAP(pos2, 1, K.get_n()));								
	//U= SWAP Uswapped SWAP+, con U prod_kron di gate di due qubit e identità
	vector<vector<complex<double>>> U=mat_prod(mySWAP, mat_prod(kron_prod(gate, ID(pow(2, K.get_n()-2))), mySWAP, "Badjoint"));
	//calcolo il nuovo K'=UK
	K.set_V(apply_mat(U, K.get_V()));
	return;
}

//funzione che, data una density_matrix by-reference, misura in successione tutti i qubits producendo un N-upla come risultato (indicizzata in realtà con il label intero della base canonica relativo);
//aggiungo una condizione per cui se la ro è pura ne ricava il ket e applica la misura più snella con il solo stato puro?
int measure_output(density_matrix &Ro){
	//epsilon di deviazione dalla perfetta conservazione della probabilità che accettiamo, visto l'errore insito nelle operazioni di calcolo
	double eps=1E-06;
	//vettore dei risultati di singolo qubit
	vector<int> my_outcomes;	
	//riempio my_outcomes effettuando in successione le proiezioni dei qubits
	for(auto k : Ro.get_qbits_left()){
		my_outcomes.push_back(Ro.project(k));
	}
	if(real(Ro.trace())<1-eps||imag(Ro.trace())>0+eps) cout<<"errore: probabilità totale= "<<Ro.trace();
	return (comp(my_outcomes));
}

//applica in successione misure proiettive sui qubit sino a ottenere un elemento della base computazionale  (DA CONTROLLARE)
int measure_output(ket &K){
	//epsilon di deviazione dalla perfetta conservazione della probabilità che accettiamo, visto l'errore insito nelle operazioni di calcolo
	double eps=1E-06;
	//vettore dei risultati di singolo qubit
	vector<int> my_outcomes;	
	//riempio my_outcomes effettuando in successione le proiezioni dei qubits
	for(auto k : K.get_qbits_left()){
		my_outcomes.push_back(K.project(k));
	}
	if(psi_norm(K.get_V())<1-eps) cout<<"errore: probabilità totale= "<<psi_norm(K.get_V());
	return (comp(my_outcomes));
}

//funzione che data una (reference a una) matrice densità restituisce il vettore diagonale delle frequnze teoriche degli esiti di proiezione sugli stati di base
vector<double> freqs_theo(density_matrix &Ro){			
	vector<double> my_fs;
	for(int i=0; i<Ro.get_N(); i++){
		my_fs.push_back(real(Ro.get_M()[i][i]));
	}
	return my_fs;
}

//funzione che data una (reference a un) ket restituisce il vettore  delle frequnze teoriche degli esiti di proiezione sugli stati di base (DA CONTROLLARE)
vector<double> freqs_theo(ket &K){			
	vector<double> my_fs;
	for(int i=0; i<K.get_N(); i++){
		my_fs.push_back(norm(K.get_V()[i]));
	}
	return my_fs;
}

//funzione che, data una matrice densità by-reference e un intero n_meas, reitera n_meas volte la funzione measure_output, ottenendo un istogramma di frequenze per i diversi risultati possibili di misura
vector<double> sample(density_matrix &my_Ro, int n_meas){
	//salvo la Ro data in input prima di modificarla
	density_matrix old_Ro=my_Ro;
	//riempio l'istogramma di frequenze
	vector<double> histo(my_Ro.get_N(), 0);	
	for(int i=0; i<n_meas; i++){
		histo[measure_output(my_Ro)]+=1./n_meas;
		my_Ro=old_Ro;
	}
	return histo;		
}

// (DA CONTROLLARE) funzione che, dato un ket by-reference e un intero n_meas, reitera n_meas volte la funzione measure_output, ottenendo un istogramma di frequenze per i diversi risultati possibili di misura
vector<double> sample(ket &my_K, int n_meas){
	//salvo la Ro data in input prima di modificarla
	ket old_K=my_K;
	//riempio l'istogramma di frequenze
	vector<double> histo(my_K.get_N(), 0);	
	for(int i=0; i<n_meas; i++){
		histo[measure_output(my_K)]+=1./n_meas;
		my_K=old_K;
	}
	return histo;		
}

//funzione che manda outcome history-> frequencies histogram, N=2^n
vector<double> ry_to_gram(vector<int> ry, int N){
	vector<double> gram(N, 0);
	for(auto el:ry){
		gram[el]+=1./ry.size();
	}
	return gram;
}
vector<int> history_sample(density_matrix &my_Ro, int n_meas){
	//salvo la Ro data in input prima di modificarla
	density_matrix old_Ro=my_Ro;
	//riempio l'istogramma di frequenze
	vector<int>histo;	
	for(int i=0; i<n_meas; i++){
		histo.push_back(measure_output(my_Ro)); 
		my_Ro=old_Ro;
	}
	return histo;		
}

vector<int> history_sample(ket &my_K, int n_meas){
	//salvo la Ro data in input prima di modificarla
	ket old_K=my_K;
	//riempio l'istogramma di frequenze
	vector<int>histo;	
	for(int i=0; i<n_meas; i++){
		histo.push_back(measure_output(my_K)); 
		my_K=old_K;
	}
	return histo;		
}


//funzione che, data una matrice densità, applica D layers composti ciascuno da un sublayer costituito da gate di singolo qubit scelte randomicamente da un insieme di gates (eventualmente sottoposte a dei vincoli ad es no due gates uguali in successione), e da un sublayer di gates di due qubits, su coppie scelte casualmente tra tutte le possibili o da un sottoinsieme (rappresentato da una matrice di adiacenza), oppure applicate in una sequenza predeterminata. 	
//conviene impostare la funzione con un primo pezzo che genera casualmente o secondo determinati parametri esattamente la matrice/vettore che raccoglie la sequenza di gates da applicare, e poi una seconda parte che li applica alla ro iniziale.

//costruisce circuito
//vettore<layer>, con class layer{vettore<string> N singlequbitgate; matrice<string> NxN doublequbitgates} con gli ingressi a rappresentare a quale qbit applicare, e una funzione costruttore che riempie con stringhe prese da set di stringhe possibili, che sono keys della mappa gate della classe frequent_gates, con opportuni check che evitano 2qbgate di qubits in sè stessi e possibilità di maschere ecc... Tale oggetto salva quindi solo stringhe molto leggere, che poi un'altra funzione, di applicazione, tradurrà in effettive applicazioni lineari agendo sulla matrice 
//calcola output circuito
//struttura costituita da un sublayer di singleqbitgate e un di doublequbit gate

class layer{
	public:
		//larghezza=numero stati base comp dei qubit a cui è applicato
		int W;
		//vettore di single qubit gates
		vector<string> oneql;
		//matrice di double qubit gates
		vector<vector<string>> twoql;
		//costruttore
		layer(int my_W){
			W=my_W;
			oneql=vector<string>(W, "///");
			twoql=vector<vector<string>>(W, vector<string>(W, "///"));
		};		
};

class circuit{
	private:
		//larghezza=numero qubit a cui è applicato
		int W;
		//profondità=numero di layers
		int D;
		//vettore di layers
		vector<layer> V;
		//vettore con il set di (keys di) single qubit gates che il circuito può utilizzare in vari modi
		vector<string> oneqs;
		//vettore con il set di (keys di) double qubit gates che il circuito può utilizzare in vari modi
		vector<string> twoqs;
		//oggetto di classe frequent_gates che porta con sè il significato (la realizzazione matriciale) delle keys presenti nei vettori precedenti
		frequent_gates G;
		//vettore contenente il set di maschere (matrici di 1 e 0, matrici di adiacenza) che il circuito puo utilizzare in vari modi (analogamente ai single qubit gates)
		vector<vector<vector<int>>> masks;	
	public:
		//costruttore che prende numero qubits, numero layers, 2 vettori di gates scelti (di uno e 2 qubits), un oggetto della classe frequent_gates da cui pescare il significato dei nomi, 
		//un vettore di maschere (vettori di vettori di interi) di zeri e uni che indicano su quale coppia di qubit applicare o meno il dato gate, 3 stringhe relative a modalità per single e double qb gate e maschere
		//mode_one= "random": single qubit gates selezionati casualmente dal vettore oneqs;
				//"random_norep": come "random" ma evitando ripetizione successiva dello stesso singlequbitgate sullo stesso qubit;
				//"imported_file.txt": i single qubit gates sono importati da un file i cui elementi sono stringhe di labels di single qubit gates (presenti nella collezione o nella libreria 'frequent_gates.h'), le cui righe rappresentano in ordine i diversi layers e le cui colonne rappresentano in ordine i diversi qubits, tale modalità funziona solo per un preciso n_qubits ed n_layers
		//mode_two= "random_q": twoqubits gate selezionati casualmente in ogni coppia di qubit
				//	"random_l" : twoqubits gate selezionati casualmente ad ogni layer, ma uguali per tutto il layer
				// "random_q_norep": come l'omonimo ma evita la ripetizione successiva dello stesso gate sulla stessa coppia tra due layer adiacenti 
				// "random_l_norep": come l'omonimo ma evita la ripetizione successiva dello stesso gate sulla stessa coppia tra due layer adiacenti 
				//"sequence": twoqubits gate presi in sequenza, ripetendola periodicamente, dal vettore twoqs, layer dopo layer (uniformemente sul layer)
		//mode_mask= "random": seleziona dal vettore masks randomicamente layer dopo layer quale maschera applicare
				//"random_no_rep": come la precedente ma evita ripetizioni successive della stessa maschera
				//"sequence": applica le maschere in sequenza dal vettore masks, periodicamente sui layer
		//OSS: RANDOMIZZARE TOTALMENTE LE MASCHERE, CIOè AVERE UN OPZIONE CHE SCEGLIE CASUALMENTE QUALE COPPIA DI GATE COLLEGARE, E' SUPERFLUO POICHè LO STESSO
		//RISULTATO PUò ESSERE OTTENUTO COL CODICE PRESENTE INTRODUCENDO IN TWOQS LA MATRICE ID2 (identità sullo spazio di due qubit), CON MOLTEPLICITà VARIABILE
		//PER MODULARNE LA FREQUENZA
				
		circuit (int my_W, int my_D, vector<string> my_oneqs, vector<string>my_twoqs, frequent_gates my_G, vector<vector<vector<int>>>my_masks,
					string mode_one="random_norep", string mode_two="sequence", string mode_mask="sequence"){
			//inizializzo i parametri principali
			W=my_W;
			D=my_D;
			oneqs=my_oneqs;
			twoqs=my_twoqs;
			G=my_G;
			masks=my_masks;
			//vettore di profondità D composto di layer, per ora vuoti, di larghezza W;
			for (int i=0; i<D; i++){
				layer L(W);
				V.push_back(L);
			}
			
			//realizzo il circuito assegnando i labels secondo diverse modalità specificate dai parametri opzionali
			//le seguenti variabili serviranno a tenere traccia dei numeri casuali estratti per i diversi oggetti da assegnare, per la modalità "norep". Occupano (molto poca) memoria in più ma rendono il codice molto piu forward compatible.
			vector<int> randoneq(W, 0);
			int randonel=0;
			vector<vector<int>> randtwoq(W, vector<int>(W, 0));
			int randtwol=0;
			int randmask=0;
			//nel caso in cui importo i single qubit gates da un file esterno
			ifstream IMPORTED(mode_one);
			//scorro sui layers
			for(int l=0; l<D; l++){
				//scorro sui qubits
				for(int q1=0; q1<W; q1++){
					
					//assegno i single qubit gates
					if(mode_one=="random"){
						global_random(randoneq[q1], oneqs.size()); 						
						V[l].oneql[q1]=oneqs[randoneq[q1]];	
					}
					else if(mode_one=="random_norep"&&oneqs.size()>1){
						global_random(randoneq[q1], oneqs.size(), "norep", l);			
						V[l].oneql[q1]=oneqs[randoneq[q1]];
					}else{
						string gate_label;
						IMPORTED>>gate_label;
						V[l].oneql[q1]=gate_label;
					}
					
					//assegno i double qubits gates
					for(int q2=q1+1; q2<W; q2++){
						
						//se vettore di maschere è un vettore vuoto, si considerano gates all to all
						if(masks.size()==0){
							goto TWOQG_LAYER_GENERATOR;
						}
						
						//se il vettore di maschere ha almeno un elemento, con diverse opzioni
						else if(mode_mask=="random"){
							if(q1==0) global_random(randmask, masks.size());			
							if(masks[randmask][q1][q2]==1){
								goto TWOQG_LAYER_GENERATOR;
							}	
						}
						else if(mode_mask=="random_norep"){
							if(q1==0) global_random(randmask, masks.size(), "norep", l);     
							if(masks[randmask][q1][q2]==1){
								goto TWOQG_LAYER_GENERATOR;
							}
						}
						else if(mode_mask=="sequence"){
							if(masks[l%masks.size()][q1][q2]==1){
								goto TWOQG_LAYER_GENERATOR;
							}
						} 
						goto MASKED_OUT;
							
						TWOQG_LAYER_GENERATOR:
						if(mode_two=="random_q"){
							global_random(randtwoq[q1][q2], twoqs.size());            
							V[l].twoql[q1][q2]=twoqs[randtwoq[q1][q2]];
						}
						else if(mode_two=="random_q_norep"){
							global_random(randtwoq[q1][q2], twoqs.size(), "norep", l);
							V[l].twoql[q1][q2]=twoqs[randtwoq[q1][q2]];
						}
						
						else if(mode_two=="random_l"){
							if(q1==0) global_random(randtwol, twoqs.size());      
							V[l].twoql[q1][q2]=twoqs[randtwol];
						}
						else if(mode_two=="random_l_norep"){
							if(q1==0) global_random(randtwol, twoqs.size(), "norep", l);
							V[l].twoql[q1][q2]=twoqs[randtwol];
						}
						
						else if(mode_two=="sequence"){
							V[l].twoql[q1][q2]=twoqs[l%twoqs.size()];
						}	
						MASKED_OUT: ;	
					}
				}
			}						
		}
		
		//getter and setter
		
		int get_W(){
			return W;
		}
	
		int get_D(){
			return D;
		}
		
		vector<layer> get_V(){
			return V;
		}
		
	
		vector<string> get_oneqs(){
			return oneqs;
		}
	
		vector<string> get_twoqs(){
			return twoqs;
		}
		
		
		//funzione che applica il dato circuito a una density_matrix (VERSIONE PIU' GENERALE E PIU' LENTA)
		void apply_circuit(density_matrix &Ro){
			if(Ro.get_n()!=W){
				cout<<"matrice densità e circuito hanno dimensioni incompatibili!"<<endl;
			}else{
				//scorro sui layer di gates
				for(int l=0; l<D; l++){
					//applico in successione i single_qubit_gates a tutti i qbits
					for(int q=0; q<W; q++){
						if(V[l].oneql[q]!="ID1"&&V[l].oneql[q]!="///"){
							single_qbit_gate(Ro, G.gate[V[l].oneql[q]], Ro.get_qbits_left()[q]);
						}
					}
					//applico in successione i double_qbit_gates a tutte le coppie
					for(int q1=0; q1<W; q1++){
						for(int q2=q1+1; q2<W; q2++){
							if(V[l].twoql[q1][q2]!="ID2"&&V[l].twoql[q1][q2]!="///"){
								double_qbit_gate(Ro, G.gate[V[l].twoql[q1][q2]], Ro.get_qbits_left()[q1], Ro.get_qbits_left()[q2]);
							}
						}	
					}
				}
			}				
		}
		
		// (DA CONTROLLARE) funzione che applica il dato circuito a un KET (VERSIONE PIU' GENERALE E PIU' LENTA)
		void apply_circuit(ket &K){
			if(K.get_n()!=W){
				cout<<"matrice densità e circuito hanno dimensioni incompatibili!"<<endl;
			}else{
				//scorro sui layer di gates
				for(int l=0; l<D; l++){
					//applico in successione i single_qubit_gates a tutti i qbits
					for(int q=0; q<W; q++){
						if(V[l].oneql[q]!="ID1"&&V[l].oneql[q]!="///"){
							single_qbit_gate(K, G.gate[V[l].oneql[q]], K.get_qbits_left()[q]);
						}
					}
					//applico in successione i double_qbit_gates a tutte le coppie
					for(int q1=0; q1<W; q1++){
						for(int q2=q1+1; q2<W; q2++){
							if(V[l].twoql[q1][q2]!="ID2"&&V[l].twoql[q1][q2]!="///"){
								double_qbit_gate(K, G.gate[V[l].twoql[q1][q2]], K.get_qbits_left()[q1], K.get_qbits_left()[q2]);
							}
						}	
					}
				}
			}				
		}
		
		//funzione che applica il dato circuito a una density_matrix, ma piu veloce del precedente poichè applica i gate (singoli e doppi) in parallelo 
		//risparmiando così prodotti di kroneker inutili con le identità degli altri qubit (funziona dunque solo se non ci sono sovrapposizioni di gates)
		//il meccanismo è che per i single qubit gates faccio il prodotto di kronecker tra i gates e per i double prima li riordino in coppie successive e poi faccio i prodotti 
		void apply_circuit_parallel(density_matrix &Ro){
			if(Ro.get_n()!=W){
				cout<<"matrice densità e circuito hanno dimensioni incompatibili!"<<endl;
			}else{
				//scorro sui layer di gates
				for(int l=0; l<D; l++){
					//applico un unitaria ottenuta dalla produttoria di kronecker dei single qubit gates
					vector<vector<complex<double>>> U1={{1.}};
					for(int q=0; q<W; q++){
						if(V[l].oneql[q]!="///"){
							U1=kron_prod(U1, G.gate[V[l].oneql[q]]);
						}else{
							U1=kron_prod(U1, ID(2));
						}
					}
					Ro.set_M(mat_prod(U1, mat_prod(Ro.get_M(), U1, "Badjoint")));
					
					//applico un unitaria ottenuta dalla produttoria di kronecker dei double qubit gates dopo aver swappato i qubits in modo da averli a coppie di successivi ai quali devo applicare i gates.
					//il numero che mi rappresenta la prima posizione libera nella stringa di qubits, cioè la posizione dove inserire il prossimo che swappo
					vector<int> perm;
					vector<vector<complex<double>>> U2={{1.}};
					for(int q1=0; q1<W; q1++){
						for(int q2=q1+1; q2<W; q2++){
							if(V[l].twoql[q1][q2]!="ID2"&&V[l].twoql[q1][q2]!="///"){
								perm.push_back(q1);
								perm.push_back(q2);
								U2=kron_prod(U2, G.gate[V[l].twoql[q1][q2]]);
								//double_qbit_gate(Ro, G.gate[V[l].twoql[q1][q2]], Ro.get_qbits_left()[q1], Ro.get_qbits_left()[q2]);
							}
						}	
					} 																		
					int n_double_gates=perm.size();
					int U2_size=U2.size();
					U2=kron_prod(U2, ID(pow(2,W-n_double_gates )));							
					while(perm.size()<W){
						perm.push_back(perm.size());
					}																		
					vector<vector<int>> mySWAP=SWAP(perm);															
					//U= SWAP Uswapped SWAP+, con U prod_kron di gate di due qubit e identità
					U2=mat_prod(mySWAP, mat_prod(U2, mySWAP, "Badjoint"));						
					//calcolo la nuova ro M'=UMU+
					Ro.set_M(mat_prod(U2, mat_prod(Ro.get_M(), U2, "Badjoint")));
				}
			}				
		}
		
		// (DA CONTROLLARE) funzione che applica il dato circuito a un ket, ma piu veloce del precedente poichè applica i gate (singoli e doppi) in parallelo 
		//risparmiando così prodotti di kroneker inutili con le identità degli altri qubit (funziona dunque solo se non ci sono sovrapposizioni di gates)
		//il meccanismo è che per i single qubit gates faccio il prodotto di kronecker tra i gates e per i double prima li riordino in coppie successive e poi faccio i prodotti 
		void apply_circuit_parallel(ket &K){
			if(K.get_n()!=W){
				cout<<"matrice densità e circuito hanno dimensioni incompatibili!"<<endl;
			}else{
				//scorro sui layer di gates
				for(int l=0; l<D; l++){
					//applico un unitaria ottenuta dalla produttoria di kronecker dei single qubit gates
					vector<vector<complex<double>>> U1={{1.}};
					for(int q=0; q<W; q++){
						if(V[l].oneql[q]!="///"){
							U1=kron_prod(U1, G.gate[V[l].oneql[q]]);
						}else{
							U1=kron_prod(U1, ID(2));
						}
					}
					K.set_V(apply_mat(U1, K.get_V()));
					
					//applico un unitaria ottenuta dalla produttoria di kronecker dei double qubit gates dopo aver swappato i qubits in modo da averli a coppie di successivi ai quali devo applicare i gates.
					//il numero che mi rappresenta la prima posizione libera nella stringa di qubits, cioè la posizione dove inserire il prossimo che swappo
					vector<int> perm;
					vector<vector<complex<double>>> U2={{1.}};
					for(int q1=0; q1<W; q1++){
						for(int q2=q1+1; q2<W; q2++){
							if(V[l].twoql[q1][q2]!="ID2"&&V[l].twoql[q1][q2]!="///"){
								perm.push_back(q1);
								perm.push_back(q2);
								U2=kron_prod(U2, G.gate[V[l].twoql[q1][q2]]);
								//double_qbit_gate(Ro, G.gate[V[l].twoql[q1][q2]], Ro.get_qbits_left()[q1], Ro.get_qbits_left()[q2]);
							}
						}	
					} 																		
					int n_double_gates=perm.size();
					int U2_size=U2.size();
					U2=kron_prod(U2, ID(pow(2,W-n_double_gates )));							
					while(perm.size()<W){
						perm.push_back(perm.size());
					}																		
					vector<vector<int>> mySWAP=SWAP(perm);															
					//U= SWAP Uswapped SWAP+, con U prod_kron di gate di due qubit e identità
					U2=mat_prod(mySWAP, mat_prod(U2, mySWAP, "Badjoint"));						
					//calcolo la nuova ro M'=UMU+
					K.set_V(apply_mat(U2, K.get_V()));
				}
			}				
		}			
};

//funzione che stampa un circuito in termini di labels
void print_circuit(circuit &C, int &p_c, string title=""){ 
	string output_file=title+"_"+to_string(p_c)+".txt";
	ofstream OUTPUT(output_file);
	p_c++;
	int l=0;
	for(auto L : C.get_V()){
		OUTPUT<<"layer "<<l<<endl;
		l++;
		for(auto el : L.oneql){
			if(el=="///") {
				OUTPUT<<"ID1"<<" ";
			}else OUTPUT<<el<<" ";
		}
		OUTPUT<<endl<<endl;
		for (auto v : L.twoql){
			for (auto el : v){
				OUTPUT<<el<<" ";
			}
			OUTPUT<<endl;
		}
		OUTPUT<<endl<<endl;
	}
}

//funzione che stampa un circuito in termini di labels in ununico ostream dato come reference
void print_circuit(circuit &C, ostream &OUTPUT, string title=""){
	OUTPUT<<title<<endl; 
	int l=0;
	for(auto L : C.get_V()){
		OUTPUT<<"layer "<<l<<endl;
		l++;
		for(auto el : L.oneql){
			if(el=="///") {
				OUTPUT<<"ID1"<<" ";
			}else OUTPUT<<el<<" ";
		}
		OUTPUT<<endl<<endl;
		for (auto v : L.twoql){
			for (auto el : v){
				OUTPUT<<el<<" ";
			}
			OUTPUT<<endl;
		}
		OUTPUT<<endl<<endl;
	}
}

//Fxeb: funzione che dato un vettore di frequenze teoriche e un vettore di frequenze empiriche restituisce la Fxeb definita come 2^n <P(xi)>i-1, n numero di qubits, N=2^n
double Fxeb_lin(vector<double> fs_t, vector<double> fs_e){
	if(fs_t.size()!=fs_e.size()) cout<<"dimensioni vettori di frequenze sono incompatibili!"<<endl;
	double Fxeb=0;
	for(int i=0; i<fs_t.size(); i++){
		Fxeb+=fs_t[i]*fs_e[i];
	}
	Fxeb=fs_t.size()*Fxeb-1;
	return(Fxeb);
}

//Fxeb: funzione che dato un vettore di frequenze teoriche e uno storico dei risultati di misura (elementi di base indicizzati da interi) restituisce la Fxeb definita come 2^n <P(xi)>i-1, n numero di qubits, N=2^n
double Fxeb_lin(vector<double> fs_t, vector<int> histo, bool nozero=false){
	double Fxeb=0;
	int count=0;
	for(auto o : histo){
		if(fs_t[o]>1e-15){
			Fxeb+=fs_t[o];
			count++;
		}
	}
	if(nozero==false){
		count=histo.size();
	}
	Fxeb=fs_t.size()*Fxeb/count-1;
	return(Fxeb);
}



//funzione simile alla precedente Fxeb, cfr supplementary material p.8
double Fxeb_log_google(vector<double> fs_t, vector<double> fs_e){
	if(fs_t.size()!=fs_e.size()) cout<<"dimensioni vettori di frequenze sono incompatibili!"<<endl;
	double Fxeb=0;
	for(int i=0; i<fs_t.size(); i++){
		if(fs_t[i]>1e-15){
			Fxeb+=fs_e[i]*log(fs_t.size()*fs_t[i]);	
		}
		
	}
	Fxeb=Fxeb+0.577;
	return(Fxeb);
}

//Fxeb basica
double Fxeb_log(vector<double> fs_t, vector<double> fs_e){
	if(fs_t.size()!=fs_e.size()) cout<<"dimensioni vettori di frequenze sono incompatibili!"<<endl;
	double Fxeb=0;
	for(int i=0; i<fs_t.size(); i++){
		if(fs_e[i]>1e-15){
			Fxeb+=fs_t[i]*log(fs_e[i]);	
		}
	}
	Fxeb=-Fxeb;
	return(Fxeb);
}

/*
///SBAGLIATO {
//Fxeb: funzione simile alla precedente ma prende la matrice densità teorica e calcola la Fxeb senza approssimazioni, usando F=-Tr(Ro_t log(Ro_e))
double Fxeb_accurate_vers1(density_matrix & Ro, vector<double> fs_e){
	int N=Ro.get_N();
	if(N!=fs_e.size()) cout<<"dimensioni vettori di frequenze sono incompatibili!"<<endl;
	double Fxeb=0;
	//faccio uso di armadillo per il calcolo del logaritmo matriciale (dato che dovrei comunque passare per la diagonalizzazione)
	Mat<complex<double>> R(N,N);
	Mat<complex<double>> E(N,N, fill::zeros);
	for(int i=0;i<N; i++){
		E(i, i)=fs_e[i];
		for(int j=0; j<N; j++){
			R(i, j)=Ro.get_M()[i][j];
		}	
	}
	complex<double> entro=-arma::trace(E*logmat(R));
	if(imag(entro)>1E-06) cout<<"Il risultato dell'entropia di von neumann non è reale!"<<endl;
	Fxeb=log2(exp(1.))*real(entro);
	return(Fxeb);
}

//Fxeb SECONDO TENTATIVO: funzione simile alla precedente ma prende la matrice densità teorica e calcola la Fxeb senza approssimazioni, usando F=-Tr(Ro_t log(Ro_e))
double Fxeb_accurate(density_matrix & Ro, vector<double> fs_e){
	int N=Ro.get_N();
	if(N!=fs_e.size()) cout<<"dimensioni vettori di frequenze sono incompatibili!"<<endl;
	double Fxeb=0;
	//faccio uso di armadillo per il calcolo del logaritmo matriciale (dato che dovrei comunque passare per la diagonalizzazione)
	Mat<complex<double>> R(N,N);
	Mat<complex<double>> E(N,N, fill::zeros);
	for(int i=0;i<N; i++){
		E(i, i)=fs_e[i];
		for(int j=0; j<N; j++){
			R(i, j)=Ro.get_M()[i][j];
		}	
	}
	vec eigvalR;
	cx_mat eigvecR;
	eig_sym(eigvalR, eigvecR,  R);
	for(int i=0; i<N; i++){
		R(i, i)=log2(real(R(i,i)));
	}
	R=eigvecR*R*inv(eigvecR);
	complex<double> entro=-arma::trace(E*R);
	if(imag(entro)>1E-06) cout<<"Il risultato dell'entropia di von neumann non è reale!"<<endl;
	Fxeb=real(entro);
	return(Fxeb);
}

///   } SBAGLIATO 
*/


/*	
//funzione che reitera apply_circuit dato il numero di istanze, il numero di misure per istanza, una referenza alla matrice densità, e i parametri del circuito  per trarre stime di altro tipo eccetera....
void iterate_circuit(int n_istanze, int n_meas, density_matrix &Ro, int my_W, int my_D, vector<string> my_oneqs, vector<string>my_twoqs, frequent_gates my_G, 
						vector<vector<vector<int>>>my_masks, string mode_one="random_no_rep", string mode_two="sequence", string mode_mask="sequence"){
	//vettore che conterrà i sample(vettori di frequenze ottenuti dalla misura sulla matrice risultante dal circuito)
	vector<vector<double>> samples;
	density_matrix temp_Ro;
	for(int it_istanza=0; it_istanza<n_istanze; it_istanza++){
		temp_Ro=Ro;
		//istanza circuito
		circuit my_C(my_W, my_D, my_oneqs, my_twoqs, my_G, my_masks, "random_no_rep", "sequence", "sequence");
		//apply_circuit
		my_C.apply_circuit(temp_Ro);
		//misura cose
		samples.push_back(sample(temp_Ro, n_meas));
	}
}
*/
	
//funzione che calcola l'entropia di von neumann data una referenza a una matrice densità
//passa attraverso la diagonalizzazione

//funzione che produce maschere/matrici di adiacenza utili dato numero di qubit e modalità
vector<vector<int>> frequent_masks(int n_qubits, string mode){
	vector<vector<int>> my_mask(n_qubits, vector<int>(n_qubits, 0));
	if(mode=="even_alternate_A"){
		if(n_qubits%2!=0) cout<<"tipo di maschera incompatibile con il numero di qubit!"<<endl;
		for(int i=0; i<n_qubits-1; i+=2){
			my_mask[i][i+1]=1;	
		}
	}
	else if(mode=="even_alternate_B"){
		if(n_qubits%2!=0) cout<<"tipo di maschera incompatibile con il numero di qubit!"<<endl;
		for(int i=1; i<n_qubits-2; i+=2){
			my_mask[i][i+1]=1;
		}
		my_mask[0][n_qubits-1]=1;
	}
	/*
	else if(){
	}
	...
	*/
	
	return my_mask;
}




