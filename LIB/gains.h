void get_gains(double gains[5], double voltage[5])
{
	std::vector<TF1*> _myfit;

	for (int i=0; i<5; i++)  _myfit.push_back(new TF1(Form("%s %i","myfit_",i),"pow(10,[0])*pow(x,[1])"));  //<-- RANGE!
	_myfit[0]->SetParameter(0,-27.59);
	_myfit[0]->SetParameter(1,10.88);
	_myfit[1]->SetParameter(0,-26.57);
	_myfit[1]->SetParameter(1,10.59);
	_myfit[2]->SetParameter(0,-27.11);
	_myfit[2]->SetParameter(1,10.9);	
	_myfit[3]->SetParameter(0,-25.86);
	_myfit[3]->SetParameter(1,10.53);
	_myfit[4]->SetParameter(0,-25.1);
	_myfit[4]->SetParameter(1,10.17);

	for (int i=0; i<5; i++) gains[i]=_myfit[i]->Eval(voltage[i]);
	for (int i=0; i<5; i++)  cout << " Gain of channel "<< i << " at "<< voltage[i] << "V: "<< gains[i] << endl;

}
