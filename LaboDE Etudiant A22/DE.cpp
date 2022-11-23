#include "Entete.h"
#pragma comment (lib,"EvoDiffDLL.lib")  
//%%%%%%%%%%%%%%%%%%%%%%%%% IMPORTANT: %%%%%%%%%%%%%%%%%%%%%%%%% 
//Les fichiers de la DLL (DiffEvolutionDLL.dll et DiffEvolutionDLL.lib) doivent se trouver dans le même répertoire que l'exécutable (.exe) et 
//dans le répertoire courant du projet pour une exécution à l'aide du compilateur.
//Indiquer les arguments du programme dans les propriétés du projet - débogage - arguements.
//Sinon, utiliser le répertoire execution.

//*****************************************************************************************
// Prototype des fonctions se trouvant dans la DLL 
//*****************************************************************************************
//DESCRIPTION: Dimension des vecteurs de la population & initialisation des solutions avec des valeurs entre Xmin et Xmax.
extern "C" _declspec(dllimport) void InitialisationPopulation(std::vector<tSolution> &unePop, tProblem unProb, tAlgoDE &unDE);

//DESCRIPTION: SÉLECTION: On retient le meilleur vecteur entre le vecteur trial (essai) et le vecteur target (cible) pour prendre la place de ce dernier dans la population
//PARAMETRES: unTarget-Vecteur cible qui pourra être remplacé par le vecteur trial si ce dernier présente une meilleure fonction objectif, iTarget-Indice dans la population du vecteur cible 
//(nécessaire à la MAJ de la Best), unTrial-Vecteur d'essai, uneBest-Meilleure solution depuis le début de l'algorithme, iBest-Indice dans la population de la meilleure solution qui pourra être modifié
//unProb-Définition du probleme, unDE-Définition de l'algorithme
extern "C" _declspec(dllimport) void Selection(tSolution &unTarget, int iTarget, tSolution unTrial, tSolution uneBest, int &iBest, tProblem unProb, tAlgoDE unDE);

//DESCRIPTION: Renvoie un double aléatoire entre a et b
extern "C" _declspec(dllimport) double AleaDouble(double a, double b);

//DESCRIPTION: Fonction qui affiche le détail des solutions (de Debut jusqu'a Fin-1) de la population
extern "C" _declspec(dllimport) void AfficherSolutions(std::vector<tSolution> unePop, int Debut, int Fin, int Iter, tProblem unProb);

//DESCRITPION: Fonction qui affiche une solution (bool Detail: avec ou sans détail)
extern "C" _declspec(dllimport) void AfficherUneSolution(tSolution P, int Iter, tProblem unProb, bool Detail);
//DESCRIPTION: Fonction qui affiche une solution dans un fichier de sortie (bool Detail: avec ou sans détail)
extern "C" _declspec(dllimport) void AfficherUneSolutionFichier(tSolution Sol, int Iter, tProblem unProb, bool Detail, ofstream &Fichier);

//DESCRIPTION: Fonction affichant les résultats de l'algorithme
extern "C" _declspec(dllimport) void AfficherResultats(tSolution uneBest, tProblem unProb, tAlgoDE unDE);
//DESCRIPTION: Fonction affichant les résultats de l'algorithme dans un fichier texte
extern "C" _declspec(dllimport) void AfficherResultatsFichier(tSolution uneBest, tProblem unProb, tAlgoDE unDE, std::string FileName);

//**Liberation de la mémoire allouée dynamiquement
extern "C" _declspec(dllimport) void LibererMemoireFinPgm(std::vector<tSolution> &unePop, tProblem unProb, tAlgoDE unDE);

//*****************************************************************************************
// Prototype des fonctions locales
//*****************************************************************************************
void InitialisationDomaineVariable(tProblem &unProb);
void EvaluationSolution(tSolution &Sol, tProblem unProb, tAlgoDE &unDE);
void Mutation(std::vector<tSolution> unePop, int iTarget, int iBest, tSolution &unMutant, tProblem unProb, tAlgoDE unDE);
void Croisement(tSolution unTarget, tSolution unMutant, tSolution &unTrial, tProblem unProb, tAlgoDE &unDE);

//******************************************************************************************
// Fonction main
//******************************************************************************************
int main(int NbParam, char *Param[])
{
	tProblem LeProb;					//**Définition de l'instance de problème
	tAlgoDE LeDE;						//**Définition des paramètres de l'algorithme
	std::vector<tSolution> Pop;			//**Ensemble de solutions 
	tSolution Mutant, Trial;			//**Vecteurs: mutant(donneur) et essai
	int NoBest;							//**Indice dans la Pop de la Meilleure solution depuis le début de l'algorithme
	int i;
	
		//**Lecture des paramètres
	LeDE.NP				= atoi(Param[1]);
	LeDE.F				= atof(Param[2]);
	LeDE.CR				= atof(Param[3]);
	LeDE.NB_EVAL_MAX	= atoi(Param[4]);
	LeDE.Iter			= 0;
	LeDE.CptEval		= 0;
	
	srand((unsigned)time(NULL));							//**Precise un germe pour le generateur aleatoire
	cout.setf(ios::fixed | ios::showpoint);

	//**Choix de la stratégie de mutation/croisement
	LeDE.TypeMut = RAND1;									//**Spécifie le type de Mutation (sélection solutions + #perturbations) - Voir Entete.h
	LeDE.TypeCr = EXP;										//**Spécifie le type de croisement  - Voir Entete.h
	
	//**Spécifications du problème à traiter
	LeProb.Fonction = BOOTH;								//**Spécifie le problème traité  - Voir Entete.h
	InitialisationDomaineVariable(LeProb);

	//**Dimension du vecteur de la population, initialisation des solutions avec des valeurs entre Xmin et Xmax
	InitialisationPopulation(Pop, LeProb, LeDE);			//**NB: Pop est un vecteur de 0 à NP - 1
	//**Evaluation des solutions et conservation de la meilleure solution
	for (i = 0; i < LeDE.NP; i++)
	{
		EvaluationSolution(Pop[i], LeProb, LeDE);
		//Conservation de la meilleure solution initiale
		if (i == 0)
			NoBest = i;
		else
			if (Pop[i].FctObj <= Pop[NoBest].FctObj)
				NoBest = i;
	}

	//Dimension des vecteurs: mutant et essai
	Mutant.X.resize(LeProb.D);
	Trial.X.resize(LeProb.D);								
	
	//AfficherSolutions(Pop, 0, LeDE.NP, LeDE.Iter, LeProb);
	AfficherUneSolution(Pop[NoBest], LeDE.Iter, LeProb, false);

	//**Boucle principale de l'algorithme
	while (LeDE.CptEval < LeDE.NB_EVAL_MAX) 	//**NE PAS ENLEVER/MODIFIER LA CONDITION SUR LE NOMBRE D'ÉVALUATION
	{
		LeDE.Iter++;
		
		for (i = 0; i < LeDE.NP; i++)
		{
			//**MUTATION: Création du vecteur mutant
			Mutation(Pop, i, NoBest, Mutant, LeProb, LeDE);

			//**CROISEMENT: Création du vecteur trial (essai)
			Croisement(Pop[i], Mutant, Trial, LeProb, LeDE);

			//**SELECTION: entre le vecteur target(cible) et le vecteur trial (essai)
			Selection(Pop[i], i, Trial, Pop[NoBest], NoBest, LeProb, LeDE);
		}

		//AfficherSolutions(Pop, 0, LeDE.NP, LeDE.Iter, LeProb);
		AfficherUneSolution(Pop[NoBest], LeDE.Iter, LeProb, false);
	};

	AfficherResultats(Pop[NoBest], LeProb, LeDE);		//**NE PAS ENLEVER
	AfficherResultatsFichier(Pop[NoBest], LeProb, LeDE, "Resultats.txt");
	
	LibererMemoireFinPgm(Pop, LeProb, LeDE);

	system("PAUSE");
	return 0;
}

//**-----------------------------------------------------------------------
//**Détermine l'intervalle de recherche selon la fonction choisie
void InitialisationDomaineVariable(tProblem &unProb)
{
	switch(unProb.Fonction)
	{
		case BOOTH:		unProb.Xmin = -10.0;	unProb.Xmax = 10.0; unProb.D = 2; break;
		case SPHERE:	unProb.Xmin = -5.12;	unProb.Xmax = 5.12; unProb.D = 10; break;
		default:		unProb.Xmin = 0.0;		unProb.Xmax = 0.0;	unProb.D = 0; break;
	}
}

//**-----------------------------------------------------------------------
//**Calcul de la fonction objectif d'une solution selon la fonction continue (problème) sélectionnée
void EvaluationSolution(tSolution &Sol, tProblem unProb, tAlgoDE &unDE)
{
	double valeur = 0.0;
	int d;

	switch (unProb.Fonction)
	{
		case BOOTH: //Rastringin: Min 0 en (0, 0 ... 0)
			valeur = pow(Sol.X[0] + 2 * Sol.X[1] - 7, 2);
			valeur += pow(2 * Sol.X[0] + Sol.X[1] - 5, 2);
			break;
	
		case SPHERE: //Sphere: Min 0 en (0,0 ... 0)
			for (d = 0; d < unProb.D; d++)
			{
				valeur +=  pow(Sol.X[d], 2);
			}
			break;

		default: valeur = FLT_MAX;
	}
	//Ne pas enlever/modifier
	Sol.FctObj = valeur;
	unDE.CptEval++;
}

//**-----------------------------------------------------------------------
//MUTATION: Creation du vecteur mutant à l'aide d'autres vecteurs de la population (ces vecteurs doivent être différents).
//PARAMETRES: unePop-Ensemble des solutions, iTarget-Indice dans la population du vecteur cible impliqué dans la mutation, iBest-Indice dans la population de la meilleure solution
//unMutant-Vecteur mutant qui sera produit et retourné, unProb-Définition du probleme, unDE-Définition de l'algorithme
void Mutation(std::vector<tSolution> unePop, int iTarget, int iBest, tSolution &unMutant, tProblem unProb, tAlgoDE unDE)
{
	int R1, R2, R3,			//indices des solutions choisies aléatoirement
		d;				

	/**********************************************************************************************************************************************/
	/*NB: Pour simplification: ne pas faire la vérification que la meilleure solution et le vecteur Cible (Target) sont des solutions différentes.*/
	/*    Faire les vérifications pour les autres vecteurs                                                                                        */
	/**********************************************************************************************************************************************/
	switch (unDE.TypeMut)
	{
		case RAND1:	//le vecteur mutant est crée en ajoutant 1 perturbation à l'aide de 3 solutions choisies aléatoirement (1 différence pondérée)
			do R1 = rand() % unDE.NP; while (R1 == iTarget);
			do R2 = rand() % unDE.NP; while (R2 == iTarget || R2 == R1);
			do R3 = rand() % unDE.NP; while (R3 == iTarget || R3 == R1 || R3 == R2);
			for (d = 0; d < unProb.D; d++)
				unMutant.X[d] = unePop[R1].X[d] + unDE.F * (unePop[R2].X[d] - unePop[R3].X[d]);
			break;
		case BEST2:	//le vecteur mutant est créé en ajoutant une perturbation à Best à travers 2 différences pondérées de solutions sélectionnées aléatoirement
			break;
		case CURRENTtoBEST1: //Le vecteur mutant est créé à l’aide de deux vecteurs choisis au hasard, ainsi que du meilleur vecteur (Best)
			break;
	}
	//Confinement d'intervalle pour chaque dimension de la solution
	for (int d = 0; d < unProb.D; d++)
	{
		if (unMutant.X[d] < unProb.Xmin || unMutant.X[d] > unProb.Xmax) unMutant.X[d] = AleaDouble(unProb.Xmin, unProb.Xmax);
	}
}

//**-----------------------------------------------------------------------
//CROISEMENT: Échange de genes entre le vecteur Target (cible) et le vecteur Mutant afin de créer le vecteur Trial (essai-progéniture)
//PARAMETRES: unTarget-Vecteur cible, unMutant-Vecteur mutant, unTrial-Vecteur d'essai qui sera produit et retourné, 
//unProb-Définition du probleme, unDE-Définition de l'algorithme
void Croisement(tSolution unTarget, tSolution unMutant, tSolution &unTrial, tProblem unProb, tAlgoDE &unDE)
{
	int n,			//croisement exponentionel: Point de départ (point coupure) [0,D-1]
		L;			//croisement exponentionel: Nombre de dimensions à copier
	double Alea;

	switch (unDE.TypeCr)
	{
			case BIN:	//Croisement de type binomial (ou uniform)
				break;
			case EXP:	//Croisement de type exponentiel (ou 2-point)
				unTrial = unTarget;  //Copie du vecteur cible dans le vecteur essai
				n = rand() % unProb.D;
				L = 0;
				do
				{
					unTrial.X[n] = unMutant.X[n];
					n = (n+1) % unProb.D; 
					L++;
					Alea = AleaDouble(0, 1);
				} while (Alea<unDE.CR && L < unProb.D);
				break;
	}

	//Evaluation de la fonction objectif - Ne pas enlever
	EvaluationSolution(unTrial, unProb, unDE);
}