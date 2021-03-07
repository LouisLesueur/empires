#include "politics.hpp"

using namespace std;



int main( int argc, char** argv  )
{

	World world("../maps/greece/", 100, 10);

	City cit(Vector2i(107,150), "City");
	State stat("State");
	City cit2(Vector2i(50,50), "City2");
	State stat2("State2");

	world.add_city_state(cit, stat);
	world.add_city_state(cit2, stat2);
	char keyboard = ' ';

	while (keyboard != 'q') {
		cout<<"---------"<<endl;
		cout<<"Cities "<<world.nCities()<<" States "<<world.nStates()<<endl;
		world.update_resources();
		world.update_pop();
		world.expand_provinces(10);
		cv::Mat out = world.get_states_image();
		cv::imshow("map", out);
		keyboard = cv::waitKey(0);
	}

	cv::destroyWindow("map");

}
