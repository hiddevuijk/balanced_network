#ifndef GUARD_read_input_h
#define GUARD_read_input_h

#include "network_struct.h"
#include <fstream>
#include <string>


std::string remove_text(std::string str)
{
	std::size_t pos = str.find('=');
	return str.substr(pos+1);
}

void read_input(Network& NW,int& N, int& K, double& theta_e, double& theta_i, double& D, double& mo,double& mf, int& ti, int& tf,int & tt, int& seed, std::string infile_name)
{
	ifstream input(infile_name);
	std::string temp;

	std::getline(input,temp);
	temp = remove_text(temp);
	N = std::stoi(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	K = std::stoi(temp);	

	std::getline(input,temp);
	temp = remove_text(temp);
	theta_e = std::stod(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	theta_i = std::stod(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	D = std::stod(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	NW.Jeo = std::stod(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	NW.Jee = std::stod(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	NW.Jie = std::stod(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	NW.Jio = std::stod(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	NW.Jei = std::stod(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	NW.Jii = std::stod(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	NW.tau = std::stod(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	mo = std::stod(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	mf = std::stod(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	ti = std::stoi(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	tf = std::stoi(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	tt = std::stoi(temp);

	std::getline(input,temp);
	temp = remove_text(temp);
	seed = std::stoi(temp);
	

}



#endif


