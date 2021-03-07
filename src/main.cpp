#include "politics.hpp"
#include <iomanip>

using namespace std;


int main( int argc, char** argv  )
{
	bool SAVE_FIG = true;
	string FIG_PATH = "figs/";
	int step = 0;

	World world("../maps/greece/", 100, 10);

	City cit(Vector2i(107,150), "City");
	State stat("State");
	City cit2(Vector2i(50,50), "City2");
	State stat2("State2");
	City cit3(Vector2i(150,150), "City3");
	State stat3("State3");

	world.add_city_state(cit, stat);
	world.add_city_state(cit2, stat2);
	world.add_city_state(cit3, stat3);
	char keyboard = ' ';

	while (keyboard != 'q') {
		cout<<"---------"<<endl;
		cout<<"Cities "<<world.nCities()<<" States "<<world.nStates()<<endl;
		world.update_resources();
		world.update_pop();
		world.expand_provinces(10);
		cv::Mat out = world.get_states_image();
		cv::imshow("map", out);

		std::stringstream stream;
		stream << std::setw(10) << std::setfill('0') << step;
		std::string step_string = stream.str();
		
		if(SAVE_FIG)
			cv::imwrite(FIG_PATH+"out_"+step_string+".png", out);
		step++;
		keyboard = cv::waitKey(0);
	}

	cv::destroyWindow("map");

}
