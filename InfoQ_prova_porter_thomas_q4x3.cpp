//main del programma per tesi info, lanciare da terminale tramite
// $ g++ InfoQ_prova.cxx InfoQ.h frequent_gates.h -std=c++11 -larmadillo
#include<iostream>
#include<cstdlib>
#include<fstream>
#include<cmath>
#include<vector>
#include<string>
#include<complex>
#include<random>
#include<ctime>
#include<set>
#include "InfoQ.h"
using namespace std;

//generale (anche stati misti)
void analisi_fxeb_density_matrix(){
	//variabili "globali"
	ofstream MY_OFS1("analisi_fxeb_confronto_back.txt");
	frequent_gates g;
	vector<string> my_oneqs={"SX", "SY", "SW"};
	vector<string> my_twoqs={"GOO"};
	//prima riga file output (labels) ignorata da gnuplot
	MY_OFS1<<"#"<<"qubits "<<"meas "<<"layers "<<"istanze "<<"average_Fxeb "<<"sigma_Fxeb "<<endl;
	//ciclo sui parametri da esplorare
	for(int n_qubits=4; n_qubits<=12; n_qubits+=2){
		//variabili che sono, in generale, funzioni del numero di qubits
		vector<vector<vector<int>>> my_masks={frequent_masks(n_qubits, "even_alternate_A"), frequent_masks(n_qubits, "even_alternate_B")};
		density_matrix Ro0(n_qubits, "allbitszero");
		density_matrix temp_Ro;	
		for(int n_meas=100; n_meas<10002; n_meas*=10){ 
			for(int n_layers=2; n_layers<20; n_layers+=3){	
				vector<double> my_Fxebs;
				double sum_Fxeb=0;
				for(int n_istanza=1; n_istanza<=100; n_istanza++){ 
					temp_Ro=Ro0;
					circuit my_C(n_qubits, n_layers, my_oneqs, my_twoqs, g, my_masks, "random_norep", "sequence", "sequence");
					my_C.apply_circuit_parallel(temp_Ro);
					vector<double> my_fs_e=sample(temp_Ro, n_meas);
					vector<double> my_fs_t=freqs_theo(temp_Ro);
					double my_Fxeb=Fxeb_log(my_fs_t, my_fs_e);
					my_Fxebs.push_back(my_Fxeb);
					sum_Fxeb+=my_Fxeb;
					double average_Fxeb=sum_Fxeb/n_istanza;
					double sigma_Fxeb=0;
					if(n_istanza>1){
						for(auto f : my_Fxebs){
							sigma_Fxeb+=(pow(f-average_Fxeb,2));
						}
						sigma_Fxeb/=(n_istanza-1);
						sigma_Fxeb=sqrt(sigma_Fxeb);
					}
					MY_OFS1<<n_qubits<<"	"<<n_meas<<"	"<<n_layers<<"	"<<n_istanza<<"	"<<average_Fxeb<<"	"<<sigma_Fxeb<<endl;
					cout<<"time: "<<clock()<<endl;	
				}
			}	
		}
	}
}

//come la precedente ma con stati solo puri
void analisi_fxeb_ket(){
	//variabili "globali"
	ofstream MY_OFS1("fs_t_q4x3_l30.txt");
	//ofstream MY_OFS2("fs_e_q4x3_l30.txt");
	
	vector<int> different_n_qubits={12};
	vector<int> different_n_meas={2};
	vector<int> different_n_layers={30};
	int n_istanze=1;
	
	frequent_gates g;
	vector<string> my_oneqs={"SX", "SY", "SW"};
	vector<string> my_twoqs={"GOO"};
	
	//prima riga file output (labels) ignorata da gnuplot
	MY_OFS1<<"#frequenze teoriche circuito 4x3 qubits, 30 layers"<<endl;
	//MY_OFS2<<"#frequenze misurate circuito 4x3 qubits, 30 layers, 10k misure"<<endl;
	//ciclo sui parametri da esplorare
	//for(int n_qubits=4; n_qubits<=8; n_qubits+=2){
	for(int n_qubits : different_n_qubits){ 
		//variabili che sono, in generale, funzioni del numero di qubits
		vector<vector<int>> mask12_A=import_int_matrix("mask12_A.txt");
		vector<vector<int>> mask12_B=import_int_matrix("mask12_B.txt");
		vector<vector<int>> mask12_C=import_int_matrix("mask12_C.txt");
		vector<vector<int>> mask12_D=import_int_matrix("mask12_D.txt");
		vector<vector<vector<int>>> my_masks={mask12_A,mask12_B,mask12_C,mask12_D,mask12_C,mask12_D,mask12_A,mask12_B};
		//vector<vector<vector<int>>> my_masks={mask12_A, mask12_B, mask12_C, mask12_D, mask12_C, mask12_D,mask12_A, mask12_B};
		ket K(n_qubits, "allbitszero");
		ket temp_K;
		//for(int n_meas=30; n_meas<=10000; n_meas*=3){	
		for(int n_meas : different_n_meas){
			//for(int n_layers=1; n_layers<=15; n_layers+=1){
			for(int n_layers : different_n_layers){
				for(int n_istanza=1; n_istanza<=n_istanze; n_istanza++){ 
					temp_K=K;
					circuit my_C(n_qubits, n_layers, my_oneqs, my_twoqs, g, my_masks, "random_norep", "sequence", "sequence");
					my_C.apply_circuit_parallel(temp_K);
					//vector<double> my_fs_e=ry_to_gram(history_sample(temp_K, n_meas), temp_K.get_N()); 
					vector<double> my_fs_t=freqs_theo(temp_K); 
					sort(my_fs_t.rbegin(), my_fs_t.rend());
					for(auto el : my_fs_t) {
						MY_OFS1<<el<<" ";
					}
					MY_OFS1<<endl;
				}
			}	
		}
	}
}


//funzione che cicla su n_qubits e n_layers stampando ogni volta il vettore di frequenze teoriche ordinato in ordine decrescente (verrà poi fittato con esponenziale da altro codice che valuterà il X^2)
void analisi_distr_freqs_theo(){
	//variabili "globali"
	string title="analisi_distr_freqs_theo_X2";
	int n_istanze=1;
	frequent_gates g;
	vector<string> my_oneqs={"SX", "SY", "SW"};
	vector<string> my_twoqs={"GOO"};
	
	for(int n_istanza=0; n_istanza<n_istanze; n_istanza++){
		string filename;
		if(n_istanze>1){
			filename=title+to_string(n_istanza)+".txt";
		}else{
			filename=title+".txt";
		}
		ofstream MY_OFS1(filename);
		//prima riga file output (labels) ignorata da gnuplot
		MY_OFS1<<"#"<<"qubits "<<"layers "<<"freqs_theo "<<endl;
		//ciclo sui parametri da esplorare
		for(int n_qubits=2; n_qubits<=12; n_qubits+=2){
			vector<vector<vector<int>>> my_masks={frequent_masks(n_qubits, "even_alternate_A"), frequent_masks(n_qubits, "even_alternate_B")};
			ket K(n_qubits, "allbitszero");
			ket temp_K;	
			for(int n_layers=2; n_layers<=27; n_layers+=5){
				temp_K=K;
				circuit my_C(n_qubits, n_layers, my_oneqs, my_twoqs, g, my_masks, "random_norep", "sequence", "sequence");
				my_C.apply_circuit_parallel(temp_K);
				vector<double> my_fs_t=freqs_theo(temp_K);
				sort(my_fs_t.rbegin(), my_fs_t.rend());
				MY_OFS1<<n_qubits<<" "<<n_layers<<" ";
				for(auto el : my_fs_t) {
					MY_OFS1<<el<<" ";
				}
				MY_OFS1<<endl;
				cout<<"time: "<<clock()<<endl;
			}
		}
	}	
}

//funzione per testare la Fxeb su distribuzioni Porter-Thomas generate ad hoc invece che dal circuito.
void test_fxeb_montecarlo(){
	ofstream MY_OFS1("test_fxeb_exp_esatta.txt");
	ofstream MY_OFS2("test_fxeb_exp_montecarlo.txt");
	MY_OFS1<<"#"<<"n_qubits "<<"fXEB-lin "<<"fXEB-log "<<endl;
	MY_OFS2<<"#"<<"n_qubits "<<"n_meas "<<"fXEB-lin "<<"fXEB-log "<<endl;
	for(int n_qubits=4; n_qubits<=12; n_qubits+=2){
		//creo istogramma PT_exact che riproduce esattamente la P-T
		int N=pow(2,n_qubits);
		vector<double> PT_exact;
		for(int i=0; i<N;i++){
			PT_exact.push_back(exp_distro(i, N));
		}
		MY_OFS1<<n_qubits<<"	"<<Fxeb_lin(PT_exact, PT_exact)<<"	"<<Fxeb_log(PT_exact, PT_exact)<<"	";
		for(auto el : PT_exact) {
			MY_OFS1<<el<<" ";
		}
		MY_OFS1<<endl;
		//creo istogramma PT_random che approssima la P-T, generato randomicamente tramite metodo accetta-nega, con dimensione del sampling denominata n_meas
		double PT_MAX=1; //(exp(N)-1)/(exp(N)*(1-exp(-N*N))); cout<<"PT_MAX: "<<PT_MAX<<endl;  
		vector<int> different_n_meas={30, 100, 300, 1000, 3000};	
		for(int n_meas : different_n_meas){
			vector<double> PT_random(N, 0);
			for(int m=0; m<n_meas; m++){
				while(true){
					int i_trial=rand()%N;
					double f_trial=(double)rand()/RAND_MAX*PT_MAX;
					if(f_trial<=exp_distro(i_trial, N)){
						PT_random[i_trial]+=1./n_meas;
						break;
					}
				}		
			}
			MY_OFS2<<n_qubits<<"	"<<n_meas<<"	"<<Fxeb_lin(PT_random, PT_random)<<"	"<<Fxeb_log(PT_random, PT_random)<<"	";
			for(auto el : PT_random) {
				MY_OFS2<<el<<" ";
			}
			MY_OFS2<<endl;
		}
	}
}

int main(){	
	clock_t my_t=clock();
	cout<<"start time: "<<clock()<<endl;
	//seed per tutti i metodi pseudocasuali da usare evenualmente nel prosieguo del programma
	srand(time(NULL));
	//titolo generale dei files di output che voglio produre con questa run del codice (i risultati saranno "titolo"+"p_c" con i che aumenta ad ogni stampa di qualcosa)
	string my_of="test_proiezione_senza_traccia";
	//ofstream unico per alcune cose che voglio stampare assieme
	
	//per stampare sul file unico 
	//ofstream MY_OFS1("analisi_fxeb.txt");
	//per stampare sul terminale
	auto& MY_OFS=cout;
	
	
	//variabile "print count" che conta il numero di stampe
	int p_c=1;
	//classe i cui membri sono gate frequentemente utilizzati
	frequent_gates g;
//______________________________________________________________________________________________________________________________________________________________
	//set di qbits che voglio tracciare via
	//set<int> s={1,3};
	//print_mat(g.gate["CNOT"], p_c, "terminal");
	//matrice densità dello stato di bell di 2 qubit
	//vector<vector<complex<double>>> BELL={{0.5, 0,0,0.5},{0,0,0,0},{0,0,0,0},{0.5,0,0,0.5}};
	//importo una matrice arbitraria da un file
	//vector<vector<complex<double>>> my_M=import_complex_matrix("my_matrix.txt");
	//stampo a terminale la matrice appena importata
	//print_mat(my_M, p_c);
	
	//density_matrix Ro1(12, "allbitszero");
	//print_mat(Ro1.get_M(), p_c, my_of);
	//cout<<"traccia: "<<Ro1.trace()<<endl;
	//aggiungo alla libreria il doublequbitgate custom "GOO" che è una combinazione in qualche modo di iSWAP e CZ
	//ma è una somma...così?
	//g.add_gate("WRONGGOO", mat_sum(kron_prod(g.gate["iSWAP"], vector<vector<complex<double>>>(1, vector<complex<double>>(1, 5./6.))), kron_prod(g.gate["CZ"], vector<vector<complex<double>>>(1, vector<complex<double>>(1, 1./6.))) ));
	//print_mat(g.gate["GOO"], p_c, "prova");
	//oppure una composizione.. ma in tal caso in che senso 1/6????
	g.add_gate("WRONGGOO", {{1,0,0,0},{0,0,{0,1},0},{0,{0,1},0,0},{0,0,0,{0,1}}});
	//oppure cos'altro?
	//inizializzo i vettori di gates di singolo e doppio qubit che utilizzerò nel circuito, la mappa nome-matrice è presente nel file frequent_gates.h
	vector<string> my_oneqs={"SX", "SY", "SW"};//{"NOT", "H", "SX", "SY", "SW", "ID1"};
	vector<string> my_twoqs={"GOO"};//{"iSWAP", "iSWAP", "iSWAP", "CZ"};//{"GOO"};//{"CNOT", "CZ", "iSWAP", "ID2"};
	
	/*
	//set di maschere per 12 qubit in grado di realizzare il caso a 1, 2 e 3 colonne con connettività 2, 3 e 4 rispettivamente, vedi appunti
	vector<vector<vector<int>>> my_masks12; //={{{1,1,1},{1,1,1},{1,1,1}}};
	vector<vector<int>> mask12_A=import_int_matrix("mask12_A.txt");
	vector<vector<int>> mask12_B=import_int_matrix("mask12_B.txt");
	vector<vector<int>> mask12_C=import_int_matrix("mask12_C.txt");
	vector<vector<int>> mask12_D=import_int_matrix("mask12_D.txt");
	my_masks12={mask12_A, mask12_B, mask12_C, mask12_D};
	*/
	
		
	/*
	//istanzio un circuito senza maschera (all to all)
	circuit my_C(12, 5, my_oneqs, my_twoqs, g, my_masks, "random_norep", "random_l", "sequence");
	
	//stampo anche la forma del circuito
	print_circuit(my_C, p_c, my_of);
	//applico il circuito alla mia Ro1
	my_C.apply_circuit(Ro1);
	
	//stampo la ro risultante
	print_mat(Ro1.get_M(), p_c, my_of);
	cout<<"traccia: "<<Ro1.trace()<<endl;
	
	
	single_qbit_gate(Ro1, g.gate["H"], 1);
	print_mat(Ro1.get_M(), p_c, my_of);
	cout<<"traccia: "<<Ro1.trace()<<endl;
		
	double_qbit_gate(Ro1, g.gate["CNOT"], 1, 2);
	print_mat(Ro1.get_M(),  p_c, my_of);
	cout<<"traccia: "<<Ro1.trace()<<endl;

	
	//sample restituisce l'istogramma delle frequenze dei vari risultati della tomografia di un oggetto della classe density_matrix, specificato il numero di misure da effettuare
	print_vec(sample(Ro1, 3), "my_histo.txt");
	*/
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/*
	
	//FOR LOOP PER REITERARE n_istanze VOLTE LE OPERAZIONI PRECEDENTI, ESEGUENDO NEL CONTEMPO ALTRE MISURE
	//numero di istanze di circuiti con la stessa struttura ma generati diversamente ogni volta se nelle varie modalità random
	int n_istanze=1;
	//numero qubits (larghezza circuito)
	int n_qubits=8;
	//numero di layers (profonditò circuito)
	int n_layers=27;
	//numero di misure da eseguire su ogni matrice densità quando si applica sample()
	int n_meas=1000;
	//...............................................................................
	//creo vettore di maschere 
	vector<vector<vector<int>>> my_masks={frequent_masks(n_qubits, "even_alternate_A"), frequent_masks(n_qubits, "even_alternate_B")};
	//vector<vector<vector<int>>> my_masks={{{0, 0, 0},{0, 0, 1},{0, 0, 0}}};	
	//istanzio la matrice densità di partenza
	ket K (n_qubits, "allbitszero");
	//reitero l'istanza del circuito
	ket temp_K;
	vector<double> my_Fxebs;
	double average_Fxeb=0;
	double sigma_Fxeb=0;
	for(int i=1; i<=n_istanze; i++){
	//while(true){
		temp_K=K;
		circuit my_C(n_qubits, n_layers, my_oneqs, my_twoqs, g, my_masks, "random", "sequence", "sequence");
		my_C.apply_circuit_parallel(temp_K);		
		print_circuit(my_C, MY_OFS, "istanza "+to_string(i));
		//print_mat(temp_Ro.get_M(), MY_OFS);
		//calcolo Fxeb usando i valori teorici di frequenze nella matrice densità prodotta e quelli "euristici" ricavati dal sampling... ora è solo una valutazione delle librerie pseudocasuali con cui effettuo le misure forse, ma poi aggiungendo modelli di rumore eccetera sarà più interessante, e comunque non ho dati sperimentali
		
		vector<double> my_fs_e=sample(temp_K, n_meas);
		vector<double> my_fs_t=freqs_theo(temp_K);
		double my_Fxeb=Fxeb_log(my_fs_t, my_fs_e);
		my_Fxebs.push_back(my_Fxeb);
		MY_OFS<<"Fxeb: "<<my_Fxeb<<endl;
		average_Fxeb+=my_Fxeb;
		print_vec(my_fs_t, MY_OFS, "freq teoriche"); 
		print_vec(my_fs_e, MY_OFS, "freq euristiche");
		
		//vector<string> my_check1={"H", "ID1","ID1"};
		//if(my_C.get_V()[0].oneql==my_check1) break;
	}
	average_Fxeb/=n_istanze;
	if(n_istanze>1){
		for(auto f : my_Fxebs){
			sigma_Fxeb+=(pow(f-average_Fxeb,2));
		}
		sigma_Fxeb/=(n_istanze-1);
		sigma_Fxeb=sqrt(sigma_Fxeb);
	}
	MY_OFS<<"average_Fxeb: "<<average_Fxeb<<", sigma: "<<sigma_Fxeb<<endl;
	
	*/
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//CHIAMO LA FUNZIONE CHE ANALIZZA LA FXEB SONDANDO LO SPAZIO DEI PARAMETRI:
	analisi_fxeb_ket();
	
	
	//CHIAMO LA FUNZIONE CHE STAMPA LA DISTRIBUZIONE DI FREQUENZE TEORICHE
	//analisi_distr_freqs_theo();
	
	//CHIAMO LA FUNZIONE CHE TESTA LA FXEB GENERANDO DISTRIBUZIONI ESPONENZIALI RANDOM SENZA BISOGNO DEL CIRCUITO
	//test_fxeb_montecarlo();
	
	//print_mat(g.gate["SX"], p_c, "terminal");
	
	/*
	Ro1.trace_over_qubit(s);
	print_mat(Ro1.get_M(), p_c, my_of);
	cout<<"traccia: "<<Ro1.trace()<<endl;
	
	density_matrix Ro2(2, "manual", BELL);
	print_mat(Ro2.get_M(), p_c, my_of);
	cout<<"traccia: "<<Ro2.trace()<<endl;
	
	Ro2.project(1);
	print_mat(Ro2.get_M(),  p_c, my_of);
	cout<<"traccia: "<<Ro2.trace()<<endl;
	
	Ro2.trace_over_qubit(1);
	print_mat(Ro2.get_M(),  p_c, my_of);
	cout<<"traccia: "<<Ro2.trace()<<endl;
	*/
	my_t=clock()-my_t;
	MY_OFS<<"# clicks: "<<my_t<<" ("<<(double)my_t/CLOCKS_PER_SEC<<" secondi)"<<endl;
	return 0;
}
