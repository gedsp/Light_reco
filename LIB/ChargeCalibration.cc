double WorkFunctioneVPerElectron = 23.6; // eV/e^-
double ElectronTofC = 1.602e-4; // fC/e^-
double WorkFunctionMeVperfC = WorkFunctioneVPerElectron*1e-6 / (ElectronTofC);

double BirksA = 0.8;
double BirksK = 0.0486; // kV*g/(cm^3*MeV)
double EField = 0.5; // kV/cm
double LArTemperature = 87.; // K
double LArDensity = -0.00615*LArTemperature + 1.928; // g/cm^3

double BirksKDividedByLArDensity = BirksK/LArDensity; // kV/MeV 
  
// returns effective gain for a given run 
// from table by Philippe Cotte 
float getEffectiveGain(int run)
{
	switch(run)
	{
		case 989: return 0.647;
		case 987: return 0.672;
		case 988: return 0.676;
		case 986: return 0.682;
		case 993: return 0.726;
		case 991: return 0.735;
		case 994: return 0.735;
		case 990: return 0.742;
		case 996: return 0.804;
		case 995: return 0.81;
		case 999: return 0.881;
		case 998: return 0.907;
		case 997: return 0.928;
		case 1004: return 0.997;
		case 1000: return 1.03;
		case 1002: return 1.03;
		case 1006: return 1.19;
		case 1035: return 1.31;
		case 1036: return 1.31;
		case 1038: return 1.32;
		case 1012: return 1.36;
		case 1016: return 1.36;
		case 1037: return 1.36;
		case 1013: return 1.37;
		case 1014: return 1.38;
		case 786: return 1.49;
		case 785: return 1.53;
		case 787: return 1.53;
		case 783: return 1.56;
		case 788: return 1.56;
		case 784: return 1.57;
		case 796: return 1.66;
		case 789: return 1.67;
		case 790: return 1.67;
		case 791: return 1.68;
		case 793: return 1.72;
		case 792: return 1.73;
		case 795: return 1.73;
		case 794: return 1.76;
		case 842: return 2.07;
		case 840: return 2.15;
		default: return -1;
	}
}

