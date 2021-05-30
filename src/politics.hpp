#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <iostream>
#include <random>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <opencv2/opencv.hpp>
#include <opencv2/core/eigen.hpp>
#include <math.h>       /* exp */
#include <unordered_map>

using namespace Eigen;

static cv::Vec3b landcol(168, 187, 191);
static cv::Vec3b watercol(170, 155, 39);


static Matrix<int,8,2> DIRECTIONS((
	Matrix<int,8,2>() << 0,1, 
			     1,1,
			     1,0,
			     1,-1,
			     0,-1,
			     -1,-1,
			     -1,0,
			     -1,1).finished());



class City{
	private:
		Vector2i position;
		std::string name;
		float pop;
		float old_pop;
		float gini;
		float range;
		float defence;

	public:
		void update_pop(float inc_res);
		void update_def(float coeff);
		void update_range(float rge);

		float Pop(){return pop;}
		float OldPop(){return old_pop;}
		float Gini(){return gini;}
		float Range(){return range;}
		float Def(){return defence;}
		Vector2i Pos(){return position;}

		City(const Vector2i &pos, 
		     std::string nm,
		     float init_pop=1., float init_gini=0., float init_def=0.);
};

class CityGraph{
	private:
		int size;
		int n_cities;
		std::vector<City> city_nodes;
		MatrixXi Roads; //0: no road, 1: road, 2: road_possible

	public:
		int nCities() {return n_cities;}
		void explore_step(MatrixXi &prov, MatrixXi &terr, MatrixXi &exp, std::unordered_map<int, int> &cts, int EXPLORE=50);
		bool add_city(City cit);
		void update_pop(VectorXf inc_res);
		void update_road(int id1, int id2, int value);
		CityGraph(int sz);
		VectorXf Pops(){
			VectorXf out = VectorXf::Zero(size);
			for(int i=0; i<n_cities; i++)
				out(i) = city_nodes[i].Pop();
			return out;
		}

		Vector2i possible_war();
		City get_city(int i){return city_nodes[i];}
		int get_road(int i, int j){return Roads(i,j);}

};



class State{
	private:
		std::string name;	
		float taxe_rate;
		float military_wealth;
		float diplomacy_wealth;
		float infra_wealth;

		//Budget
		float infrastructure; //Will reduce inequalities
		float military; //Will increase attack and defense
		float diplomacy; //Will increase diplomatic power

		Vector3f old_pol;

		float power;

	public:
		State(std::string nm,
		      float init_taxe=0.1, float init_wealth=2, 
		      float init_inf=0.33, float init_mil=0.33, float init_dip=0.33);

		float apply_taxes(float inc_money);
		void change_policy(float popdiff);
		float MilWealth() {return military_wealth;}
		void update_power(float new_pow){power=new_pow;}
		float Power(){return power;}

		Vector3f policy() const{return Vector3f(infrastructure, military, diplomacy);}
};

float affinity(const State& s1, const State& s2);

class StateGraph{
	private:
		int size;
		int n_states;
		std::vector<State> state_nodes;
		MatrixXi Relations; //0: unknown, 1:agressive, 2:neutral, 3:ally
	public:
		int nStates(){return n_states;}
		StateGraph(int sz);


		MatrixXf get_money();
		float apply_taxes(int i, float inc_money);
		void update_power(MatrixXi &states, MatrixXi &canexp);
		void add_state(State stat);
		void update_dip(int s1, int s2, int value);
		int war(int s1, int s2); //return winner's id
		State getState(int i) {return state_nodes[i];}
};

class World{
	private:
		MatrixXi Map; //0: sea 1: land
		MatrixXf Topography;
		MatrixXi Provinces; //0: unknown, i:i-1
		MatrixXi States; 
		MatrixXi CanExpand; 
		MatrixXf Resources;

		float scaling;

		int width;
		int height;

		int n_cities;
		int n_states;

		std::unordered_map<int, int> city_T_state;

		CityGraph *CitiesG;
		StateGraph *StateG;
		
		std::vector<cv::Vec3b> province_palette;
		std::vector<cv::Vec3b> state_palette;

	public:
		int cols() {return Map.cols();}
		int rows() {return Map.rows();}
		bool isLand(int i, int j){
			if(Map(i,j)==1)
				return true;
			return false;
		}

		void genStateMap();

		int nCities(){return CitiesG->nCities();}
		int nStates(){return StateG->nStates();}
		void add_city_state(City cit, State stat);
		void update_resources();
		void update_pop();
		void update_policies();
		void update_diplomacy(int new_wars_per_turn );
		VectorXf resources_per_city();
		VectorXf popdiff_per_state();
		void expand_provinces(int thresh, int range_px, int new_cities_per_turn);
		World(std::string path, int cit, int stat, float coef=1);

		cv::Mat get_provinces_image();
		cv::Mat get_pop_image();
		cv::Mat get_res_image();
		cv::Mat get_states_image(bool show_boundaries=false, bool show_cities=true, bool show_roads=true);
		
};

