
/*

Macro written in root 6.04. If you are running in lxplus, just run:
	source /cvmfs/sft.cern.ch/lcg/contrib/gcc/4.8/x86_64-slc6/setup.sh
	source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.04.18/x86_64-slc6-gcc48-opt/root/bin/thisroot.sh

Then just do:
	root -l attach_all();


*/

//const string indirL= "/eos/experiment/wa105/data/311_PMT/data/root/"; // input folder for raw light-crt data files.
const string indirL= "/eos/user/c/chlastor/ReprocessedData/21Feb/"; // input folder for raw light-crt data files
const string indirC= "/eos/experiment/wa105/offline/LArSoft/Data/Reco/2018_August_24/ROOT/recofast/";// input folder for raw charge data files.

//const string outfolder = "/eos/user/j/jsotooto/root_files/CRT/"; //output folder where you want to store the matched data.
const string tempfolder = "/eos/user/l/leyton/2018Aug24/temporary/"; //output folder where you want to store the matched data.
const string outfolder = "/eos/user/l/leyton/2018Aug24/"; //output folder where you want to store the matched data.

#include "attach_lib_Aug2018.h"

void attach_all_Aug2018(int option)
{


	int light_runs[100]=
{
//1667,1668,1669,1670,1671
//1225	,
//1228	,
//1229	,
//1235	,
//1237	,
//1241	,
//1242	,
//1295	,
/*
1320	,
1321	,
1322	,
1399	,
1400	,
1406	,
1407	,
1408	,
1409	,
1666	,
1670	,
1671	,
*/
1672	
/*
1682	,
1683	,
1684	,
1685	,
1686	,
1687	,
1688	,
1689	,
1690	,
1718	,
1719	,
1720	,
1721	,
1722	,
1723	,
1724	,
1725	,
1726	,
1727	,
1728	,
1729	,
1730	,
1731	,
1732	,
1733	,
1735	,
1736	,
1737	,
1738	,
1739	,
1740	,
1741	,
1742	,
1743	,
1744	,
1745	
*/
}; //light run numbers to be assigned to the charge runs.
	//int n_lightruns=49; // size of the vector above.
	int n_lightruns=1;

	int Lstarttime[100], Lendtime[100];
	for (int i=0; i<n_lightruns;i++){ cout << i << " " << light_runs[i] << endl; timing_light(light_runs[i],Lstarttime[i],Lendtime[i],"V");}

	//charge runs to process
	int charge_run_number;

switch(option)
{
case 1:
charge_run_number=	1000	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1002	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1003	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime); // peta siempre!
charge_run_number=	1004	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1005	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1006	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1007	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1008	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1009	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime); // peta siempre!
charge_run_number=	1010	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1011	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1012	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1013	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1014	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1016	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
break;
case 2:
charge_run_number=	1035	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1036	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1037	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1038	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1039	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1165	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1166	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1167	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1172	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1173	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1174	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1175	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1176	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1177	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1178	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1180	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1181	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1182	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1183	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1184	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1185	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1186	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
break;
case 3://done!
charge_run_number=	1187	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1188	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1189	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1190	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1191	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1192	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1193	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1194	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1195	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
break;
case 4://done
charge_run_number=	1196	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1197	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1198	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	1199	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	771	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	783	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	784	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	785	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	786	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	787	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	788	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	789	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	790	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
break;
case 5: //done
charge_run_number=	791	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	792	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	793	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	794	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	795	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	796	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	798	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	799	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	801	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	833	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	834	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	835	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	836	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
break;
case 6://done!
charge_run_number=	837	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	838	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	840	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	841	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	842	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	843	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	986	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	987	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	988	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	989	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
break;
case 7:
  //charge_run_number=	990	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
  //charge_run_number=	991	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
  //charge_run_number=	992	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
  //charge_run_number=	993	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
  //charge_run_number=	994	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
charge_run_number=	995	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
//charge_run_number=	996	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
//charge_run_number=	997	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
//charge_run_number=	998	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);
//charge_run_number=	999	;attach_charge_run(charge_run_number,n_lightruns,light_runs,Lstarttime,Lendtime);

}

}
