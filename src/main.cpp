#include "politics.hpp"
#include <iomanip>
#include <stdlib.h>
#include <time.h>




using namespace std;


void init_cities(int n, World& world){

	for(int i=0; i<n; i++){

		bool keep_searching = true;

		while(keep_searching){
			int x = rand() % world.rows();
			int y = rand() % world.cols();
			
			if(world.isLand(x,y)){
				keep_searching = false;
				City cit(Vector2i(x,y), "City"+to_string(i));
				State stat("State"+to_string(i));
				world.add_city_state(cit, stat);

			}
		}
	}
}




int main( int argc, char** argv  )
{
	bool SAVE_FIG = false;
	string FIG_PATH = "figs/";
	int step = 0;
	int new_cities_per_turn = 10;
	int new_wars_per_turn = 1000;
	float coef = 0.5;
	int range_px = 50*coef;

	World world("../maps/mini_europe/", 400, 400, coef);
	
	init_cities(50, world);

	char keyboard = ' ';

	while (keyboard != 'q') {
		cout<<"---------"<<endl;
		cout<<"step: "<< step<< " Cities "<<world.nCities()<<" States "<<world.nStates()<<endl;
		world.update_resources();
		world.update_pop();
		world.expand_provinces(5, range_px, new_cities_per_turn);
		
		if(step > 1)
			world.update_policies();
		world.update_diplomacy(new_wars_per_turn);
		cv::Mat out = world.get_states_image(false, true, false);
		cv::Mat out_pop = world.get_pop_image();
		cv::imshow("map", out);
		cv::imshow("pop", out_pop);

		std::stringstream stream;
		stream << std::setw(10) << std::setfill('0') << step;
		std::string step_string = stream.str();
		
		if(SAVE_FIG)
			cv::imwrite(FIG_PATH+"out_"+step_string+".png", out);
		step++;
		keyboard = cv::waitKey(1);

	}

	cv::destroyWindow("map");

}
