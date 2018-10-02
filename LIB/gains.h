void get_gains(double gains[5], double voltage[5], bool lemsON=true)
{
	std::vector<TF1*> _myfit;

	for (int i=0; i<5; i++)  _myfit.push_back(new TF1(Form("%s %i","myfit_",i),"pow(10,[0])*pow(x,[1])"));  //<-- RANGE!
	if (lemsON)
	{
		// October calibration
		cout << "Using October 2017 calibration" << endl;
		_myfit[0]->SetParameter(0,-26.49);
		_myfit[0]->SetParameter(1,10.54);
		_myfit[1]->SetParameter(0,-27.96);
		_myfit[1]->SetParameter(1,11.03);
		_myfit[2]->SetParameter(0,-28.51);
		_myfit[2]->SetParameter(1,11.34);	
		_myfit[3]->SetParameter(0,-27.72);
		_myfit[3]->SetParameter(1,11.12);
		_myfit[4]->SetParameter(0,-27.41);
		_myfit[4]->SetParameter(1,10.89);
	}
	else
	{
		// July calibration
		cout << "Using July 2017 calibration" << endl;
		_myfit[0]->SetParameter(0,-27.69);
		_myfit[0]->SetParameter(1,10.93);
		_myfit[1]->SetParameter(0,-26.42);
		_myfit[1]->SetParameter(1,10.57);
		_myfit[2]->SetParameter(0,-27.12);
		_myfit[2]->SetParameter(1,10.92);	
		_myfit[3]->SetParameter(0,-27.21);
		_myfit[3]->SetParameter(1,10.99);
		_myfit[4]->SetParameter(0,-26.75);
		_myfit[4]->SetParameter(1,10.70);
	}
	
	// previous calibration data 
	/*
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
	*/		

	for (int i=0; i<5; i++) gains[i]=_myfit[i]->Eval(voltage[i]);
	for (int i=0; i<5; i++)  cout << " Gain of channel "<< i << " at "<< voltage[i] << "V: "<< gains[i] << endl;

}
