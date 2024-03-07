#pragma once
class GlobalData
{
public:
	int SimulationTime;
	int SimulationStepTime;
	int Conductivity;
	int Alfa;
	int Tot;
	int InitialTemp;
	int Density;
	int SpecificHeat;
	int Nodes_number;
	int Elements_number;

	GlobalData(int st,int sst, int con, int a,int t,int it,int d,int sh,int nn,int en):
		SimulationTime(st),SimulationStepTime(sst),Conductivity(con),Alfa(a),Tot(t),
		InitialTemp(it),Density(d),SpecificHeat(sh),Nodes_number(nn),Elements_number(en)
	{}
};

