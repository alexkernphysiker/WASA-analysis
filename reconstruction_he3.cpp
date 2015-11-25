// this file is distributed under 
// GPL v 3.0 license
#include <unistd.h>
#include <reconstruct.h>
Genetic::RANDOM engine;
int main(int,char**){
	Plotter::Instance().SetOutput(SimulationDataPath(),"he3reconstruction");
	printf("begin\n");
	ProcessReconstruction("He3.E.FRH1",0.00,0.60,60,[](double&,double&){return true;},engine);
	ProcessReconstruction("He3.E.FRH2",0.25,0.60,60,[](double&x,double&y){return (x>0.25)&&(y>0.45);},engine);
	ProcessReconstruction("He3.th",0,0.2,60,[](double&x,double&y){return pow(x-y-0.01,2)<0.0003;},engine);
	ProcessReconstruction("He3.phi",0,6.3,100,[](double&x,double&y){return pow(x-y,2)<0.5;},engine);
}