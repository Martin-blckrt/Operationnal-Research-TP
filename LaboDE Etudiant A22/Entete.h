#ifndef __ENTETE_H_
#define __ENTETE_H_

#include <cstdio>
#include <cstdlib> 
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <ctime>  
#include <cmath>
#include <vector>
#include <limits>
using namespace std;

enum eProb { BOOTH, SPHERE, ALPINE, RASTRIN };			//**Ne pas modifer
enum eDefinitionXY { RAND1, BEST2, CURRENTtoBEST1 };	//**Ne pas modifer
enum eDefinitionZ { BIN, EXP };							//**Ne pas modifer


struct tProblem									//**Définition du problème
{
	eProb	Fonction;							//**Nom de la fonction ou du problème à traiter
	int		D;									//**Nbre de Variables (dimensions)
	double	Xmin;								//**Domaine des variables: valeur minimale 
	double	Xmax;								//**Domaine des variables: valeur maximale
};

struct tSolution								//**Définition d'une solution: 
{
	std::vector<double>		X;					//**Valeur réelle pour chacune des dimensions de la solution. NB: Vecteur de 0 à D-1.
	double					FctObj;				//**Valeur de la fonction objectif
};

struct tAlgoDE
{
	int		NP;							//**Taille de la population (nombre de solutions)
	double	CR;							//**Probabilité de croisement [0,1]
	double	F;							//**Facteur de mutation [0,1]
	eDefinitionXY TypeMut;				//**Choix de la stratégie de croisement DE/x/y/z avec x: mode de sélection, y: nombre de différentiations pour la perturbation
	eDefinitionZ TypeCr;				//**Choix de la stratégie de croisement DE/x/y/z avec z: schéma de croissement (exp ou bin) 
	int		Iter;						//**Compteur du nombre de générations
	int		CptEval;					//**COMPTEUR DU NOMBRE DE SOLUTIONS EVALUEES. SERT POUR CRITERE D'ARRET.
	int		NB_EVAL_MAX;				//**CRITERE D'ARRET: MAXIMUM "NB_EVAL_MAX" EVALUATIONS.
};

#endif