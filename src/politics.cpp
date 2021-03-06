#include "politics.hpp"


float cumsum(int start, int end, const std::vector<float> &vec){
	float cum = 0;
	for(int i=start; i<end; i++)
		cum += vec[i];
	return cum;
}


cv::Vec3b rdmColor(){
	return cv::Vec3b(rand()%255, rand()%255, rand()%255);
}


std::vector<Vector2i> gen_path(Vector2i init_pos, MatrixXi &terr, int size, const std::vector<float> &prb){
	// 0: N, 1:NE, ...
	std::vector<Vector2i> out = {};
	int count = 0;

	Vector2i pose = init_pos;
	out.push_back(init_pos);

	while(count < size){

		float num = float(rand() % 1000000 + 1)/1000000;
		int out_id = 0;

		count ++;

		if(0 < num < prb[0])
			out_id = 0;
		else if(prb[0]<num<cumsum(0,1,prb))
			out_id = 1;
		else if(cumsum(0,1,prb)<num<cumsum(0,2,prb))
			out_id = 2;
		else if(cumsum(0,2,prb)<num<cumsum(0,3,prb))
			out_id = 3;
		else if(cumsum(0,3,prb)<num<cumsum(0,4,prb))
			out_id = 4;
		else if(cumsum(0,4,prb)<num<cumsum(0,5,prb))
			out_id = 5;
		else if(cumsum(0,5,prb)<num<cumsum(0,6,prb))
			out_id = 6;
		else if(cumsum(0,6,prb)<num<1)
			out_id = 7;

		if(0<=pose(0)+DIRECTIONS.row(out_id)(0) && pose(0)+DIRECTIONS.row(out_id)(0)<terr.rows() && 
		   0<=pose(1)+DIRECTIONS.row(out_id)(1) && pose(1)+DIRECTIONS.row(out_id)(1)<terr.cols()){
			pose += DIRECTIONS.row(out_id);
			out.push_back(pose);

		}
		else
			out.push_back(pose);

	}

	return out;
}

//------------------------------------------------------------------------------------

void City::update_pop(float inc_res){
	float death = 0.01;
	float surplus = inc_res-pop;
        float natality = 0;
       	if (-0.5<inc_res<0.5) 
		natality = 0.02*(surplus/inc_res);

	pop += (natality-death)*pop;
	gini += 0.1*natality;
}

void City::update_def(float coeff){
	defence += coeff*defence;
}

void City::update_range(float rge){
	range = rge;
}
City::City(const Vector2i &pos, 
		std::string nm, float init_pop, float init_gini, float init_def){
	name = nm;
	position = pos;
	pop = init_pop;
	gini = init_gini;
	defence = init_def;
	range=0;
}

//-------------------------------------------------------------

void CityGraph::explore_step(MatrixXi &prov, MatrixXi &terr, MatrixXi &exp, MatrixXi &cts){
	int state;
	bool search=true;

	int origin = rand() % n_cities;
	
	std::vector<float> prb = {0,0,0,0,0,0,0,0};
	int direction = rand() % 8;
	int dirb, dira;
	if(direction==0)
		dirb = 7;
	else
		dirb = direction-1;
	if(direction==7)
		dira = 0;
	else
		dira = direction+1;

	prb[direction] = 0.7;
	prb[dirb] = 0.15;
	prb[dira] = 0.15;

	int EXPLORE = 100;
	std::vector<Vector2i> path = gen_path(city_nodes[origin].Pos(), terr, EXPLORE, prb);

	bool isOk = true;
	if(terr(path[EXPLORE-1](0), path[EXPLORE-1](1))==0)
		isOk = false;
	
	else{
		for(int i=0; i<EXPLORE; i++){
			if(prov(path[i](0), path[i](1))>0 && prov(path[i](0), path[i](1)) != origin+1){
				isOk = false;
				this->update_road(origin, prov(path[i](0), path[i](1))-1, 2);
			}
		}
	}

	if(isOk){
		for(int i=0; i<EXPLORE; i++){
			if(terr(path[i](0), path[i](1))==1){
				if(prov(path[i](0), path[i](1))==0){
					if(i<int(EXPLORE/2))
						prov(path[i](0), path[i](1))=origin+1;
					else
						prov(path[i](0), path[i](1))=n_cities+1;
					terr(path[i](0), path[i](1))=2;
					exp(path[i](0), path[i](1))=1;
				}
			} else if(terr(path[i](0), path[i](1))==0){
				terr(path[i](0), path[i](1))=2;

			}
		}
	
		this->update_road(origin, n_cities, 1);
		state=0;
		search=true;
		while(search){
			if(cts(origin,state)==1)
				search=false;
			else
				state++;
		}
		
		if(this->add_city(City(path[EXPLORE-1], "test"))){
			std::cout<<state<<std::endl;
			cts(n_cities-1, state) = 1;
			exp(path[EXPLORE-1](0), path[EXPLORE-1](1)) = 1;
		}
	}

}

bool CityGraph::add_city(City cit){

	if(n_cities < size-1){
		city_nodes.push_back(cit);
		n_cities ++;
		return true;
		}
	return false;
}

void CityGraph::update_road(int id1, int id2, int value){
	Roads(id1, id2) = value;
}


CityGraph::CityGraph(int sz){
	size = sz;
	n_cities = 0;
	city_nodes = {};
	Roads = MatrixXi::Zero(sz, sz);
}

void CityGraph::update_pop(VectorXf inc_res){
	for(int i=0; i<n_cities; i++)
		city_nodes[i].update_pop(inc_res(i));
}


//-------------------------------------------------------------


State::State(std::string nm,
	     float init_taxe, float init_wealth, 
	     float init_inf, float init_mil, float init_dip){
	
	name = nm;
	taxe_rate = init_taxe;
	wealth = init_wealth;
	infrastructure = init_inf;
	military = init_mil;
	diplomacy = init_dip;
}


Vector3f State::get_money(){
	Vector3f out(infrastructure*wealth, military*wealth, diplomacy*wealth);
	wealth -= out.sum();
	return out;
}

float State::apply_taxes(float money){
	wealth = taxe_rate*money;
	return money*=(1-taxe_rate);
}

void State::change_policy(float new_taxes, 
		float new_inf, float new_mil, float new_dip){
	taxe_rate = new_taxes;
	infrastructure = new_inf;
	military = new_mil;
	diplomacy = new_dip;
}

//-----------------------------------------------------------


StateGraph::StateGraph(int sz){
	size = sz;
	n_states=0;
	state_nodes = {};
	Relations = MatrixXi::Zero(sz, sz);
}

void StateGraph::add_state(State stat){
	if(n_states < size){
		state_nodes.push_back(stat);
		n_states ++;
	}
}

MatrixXf StateGraph::get_money(){
	MatrixXf out(3, n_states);
	for(int i=0; i<n_states; i++)
		out.col(i) = state_nodes[i].get_money();
	return out;
} 


void StateGraph::apply_taxes(VectorXf &money){
	for(int i=0; i<n_states; i++)
		money(i) = state_nodes[i].apply_taxes(money(i));
}

//--------------------------------------------------------------


World::World(std::string path, int cit, int stat){
	cv::Mat terr = cv::imread(path+"bound.png", cv::IMREAD_GRAYSCALE);
	cv::cv2eigen(terr/255, Map);

	height = terr.rows;
	width = terr.cols;

	cv::Mat topo = cv::imread(path+"topo.png", cv::IMREAD_GRAYSCALE);
	cv::cv2eigen(topo/255, Topography);

	Provinces = MatrixXi::Zero(height, width);
	CanExpand = MatrixXi::Zero(height, width);
	Resources = MatrixXf::Constant(height, width, 1);

	n_cities = cit;
	n_states = stat;
	
	state_city = MatrixXi::Zero(n_cities, n_states); 
	
	CitiesG = new CityGraph(cit);
	StateG = new StateGraph(stat);
	
	province_palette = {};
	for(int i=0; i<cit; i++)
		province_palette.push_back(rdmColor());
	
	state_palette = {};
	for(int i=0; i<cit; i++)
		state_palette.push_back(rdmColor());
}

void World::add_city_state(City cit, State stat){
	CitiesG->add_city(cit);
	StateG->add_state(stat);
	state_city(CitiesG->nCities()-1, StateG->nStates()-1) = 1;
	CanExpand(cit.Pos()(0), cit.Pos()(1)) = 1;
}

void World::update_resources(){
	Resources += (MatrixXf::Constant(height, width, 1)-Topography);
}
	
VectorXf World::resources_per_city(){

	VectorXf city_res = VectorXf::Zero(n_cities);
	int res=0;
	for(int i=0; i<height; i++){
		for(int j=0; j<width; j++){
			if(Map(i,j)>=1)
				city_res(Map(i,j)-1) += Resources(i,j);
		}
	}


	return city_res;
}

void World::update_pop(){
	VectorXf res = this->resources_per_city();
	CitiesG->update_pop(res);
}

void World::expand_provinces(int thresh){
	CitiesG->explore_step(Provinces, Map, CanExpand, state_city);
	int roll;

	std::vector<Vector2i> expidx = {};
	for(int i=1; i<height-1; i++){
		for(int j=1; j<width-1; j++){
			if(CanExpand(i,j)==1){
				expidx.push_back(Vector2i(i,j));
			}
		}
	}
	
	for(int k=0; k<expidx.size(); k++){
		roll = rand() % 100;
		if(roll<thresh){
			//int dir_id = rand() % 8;
			//Vector2i dir = DIRECTIONS.row(dir_id);
			int i = expidx[k](0);
			int j = expidx[k](1);

			for(int k=0; k<8; k++){
				Vector2i dir = DIRECTIONS.row(k);

				if(Map(i+dir(0), j+dir(1))==1){
					if(Provinces(i+dir(0), j+dir(1))==0){
						Provinces(i+dir(0), j+dir(1))=Provinces(i,j);
						CanExpand(i+dir(0), j+dir(1))=1;
					}
				}
			}
		CanExpand(i,j) = 0;
			
		}
	}
}


cv::Mat World::get_provinces_image(){

	//BGR !!!!!!!!!!!
	cv::Mat img(height, width, CV_8UC3);

	for(int i=0; i<height; i++){
		for(int j=0; j<width; j++){
			if(Map(i,j) == 0){
				img.at<cv::Vec3b>(i,j) = watercol;
			} else if(Map(i,j) == 1){
				img.at<cv::Vec3b>(i,j) = landcol;
				if(Provinces(i,j)>0){
					img.at<cv::Vec3b>(i,j) = province_palette[Provinces(i,j)-1];
				}
			}
			
			if(Map(i,j) == 2){
				img.at<cv::Vec3b>(i,j) = roadcol;
			}

		}
	}

	return img;
}

cv::Mat World::get_pop_image(){
	VectorXf pops = CitiesG->Pops();
	cv::Mat img(height, width, CV_8UC3);

	for(int i=0; i<height; i++){
		for(int j=0; j<width; j++){
			if(Map(i,j) == 0){
				img.at<cv::Vec3b>(i,j) = watercol;
			} else {
				img.at<cv::Vec3b>(i,j) = landcol;
				if(Provinces(i,j)>0){
					uint8_t col = uint8_t((pops(Provinces(i,j)-1)/pops(pops.maxCoeff()))*255);
					img.at<cv::Vec3b>(i,j) = cv::Vec3b(col,col,col);
				}
			}

		}
	}

	return img;
}


cv::Mat World::get_states_image(){

	//BGR !!!!!!!!!!!
	cv::Mat img(height, width, CV_8UC3);

	for(int i=0; i<height; i++){
		for(int j=0; j<width; j++){
			if(Map(i,j) == 0){
				img.at<cv::Vec3b>(i,j) = watercol;
			} else if(Map(i,j) == 1){
				img.at<cv::Vec3b>(i,j) = landcol;
				if(Provinces(i,j)>0){
					int prov = Provinces(i,j)-1;
					int state = 0;
					for(int k=0; k<n_states; k++){
						if(state_city(prov, k)==1)
							state=k;
					}

					img.at<cv::Vec3b>(i,j) = state_palette[state];
				}
			}
			
			if(Map(i,j) == 2){
				img.at<cv::Vec3b>(i,j) = roadcol;
			}

		}
	}

	return img;
}
