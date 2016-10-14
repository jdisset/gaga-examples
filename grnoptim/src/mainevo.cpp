#include "external/gaga/gaga.hpp"
#include "external/grgen/classic.hpp"
#include "external/grgen/grn.hpp"
#include "external/grgen/real.hpp"

//////////////////////////////////////////////
// Code exemple montrant l'utilisation de gaga avec une classe ADN personnalisée qui est
// en fait un simple wrapper de GRN.
// Ici, l'ADN d'un individu est constitué d'un unique GRN à 2 entrées et 1 sortie. Le but
// est que le grn apprenne à fixer son output à la valeur moyenne de ses deux inputs. La
// fonction de fitness à minimiser est l'erreur moyenne.

struct MyDNA {  // ma classe ADN. Ici c'est un simple wrapper de GRN, mais on peut
	// imaginer y rajouter pleins d'autres attributs (il faudra alors penser à
	// modifier serialize, mutate et crossover).
	using GRN_t = GRN<RealCoords>;
	GRN_t grn;
	MyDNA() {
		// utile pour générer les adn random de la première génération. C'est aussi ici que
		// l'on va pouvoir définir les entrées / sorties du GRN
		grn.randomParams();   // attribue des valeurs aléatoires à Beta et Delta
		grn.randomReguls(1);  // on commence avec une régulatrice aléatoire

		// là on décrit les input et outputs de notre grn.
		grn.addRandomProtein(ProteinType::input, "in0");  // 2 inputs
		grn.addRandomProtein(ProteinType::input, "in1");
		grn.addRandomProtein(ProteinType::output, "out");  // 1 output
	}

	MyDNA(GRN_t g) : grn(g) {}  // constructeur via un grn

	/***** REQUIS PAR GAGA : ********/
	// il faut un constructeur qui puisse désérialiser un string
	MyDNA(const string& s) : grn(s) {}  // on utilise juste le constructeur de GRN

	// il faut un méthode de sérialisation qui produise un string (on utilise celle du GRN)
	std::string serialize() const { return grn.serialize(); }

	// il faut une méthode reset() qui sera appelée entre chaque génération pour
	// réinitialiser un individu
	void reset() { grn.reset(); }

	// mutation & crossover : on réutilise les méthodes du GRN
	void mutate() { grn.mutate(); }
	MyDNA crossover(const MyDNA& other) { return MyDNA(grn.crossover(other.grn)); }
};

int main(int argc, char** argv) {
	GAGA::GA<MyDNA> ga(argc, argv);  // déclaration de l'algo G; le paramètre template
	// définit quelle classe on utilisera comme ADN à faire
	// évoluer (ici la classe MyDNA)

	ga.setEvaluator(  // on passe en paramètre la fonction d'évaluation. Celle-ci peut être
	                  // une fonction définie ailleurs ou une lambda (ici c'est une lambda)
	    [](auto& individu) {  // l'évaluateur prend une instance d'Individual<DNA> en
		                        // argument.
		    double error = 0;
		    const double nbSteps = 400;
		    for (double i = 0; i < nbSteps; ++i) {
			    double v0 = sin(i * 0.05) * 0.5 + 0.5;
			    double v1 = cos(i * 0.03) * 0.5 + 0.5;
			    double target = ((v0 + v1) * 0.5);
			    individu.dna.grn.setProteinConcentration("in0", ProteinType::input, v0);
			    individu.dna.grn.setProteinConcentration("in1", ProteinType::input, v1);
			    individu.dna.grn.step(10);
			    error += std::abs(
			        individu.dna.grn.getProteinConcentration("out", ProteinType::output) -
			        target);
		    }
		    // on peut maintenant fixer la fitness de cet individu. Ici nous n'avons qu'un
		    // seul objectif : "Erreur moyenne" mais on pourrait truès bien en rajouter
		    // d'autre. On passe l'erreur en négatif car Gaga va par
		    // défaut chercher à maximiser la fitness. On pourrait aussi choisir de
		    // fonctionner en minimisation via la méthod setIsBetterMethod qui nous permet de
		    // choisir la fonction de comparaison entre deux fitness
		    individu.fitnesses["Erreur moyenne"] = -error / nbSteps;
		  },
	    "erreurMoyenne");  // le deuxième argument est le nom de l'évaluateur (facultatif)

	ga.setPopSize(200);  // population de 200 individus
	ga.setMutationProba(0.8);
	ga.setCrossoverProba(0.2);

	ga.setVerbosity(1);  // afficher les stats par génération uniquement
	ga.initPopulation(
	    []() { return MyDNA(); });  // on initialise la population avec des individus random

	ga.step(400);  // on lance l'algo G pour 400 générations
	return 0;
}
