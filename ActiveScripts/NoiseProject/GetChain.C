



TChain *GetChain(TString flist, TString ddir="./",TString tree="TreeGReco", TString suffix="", int max=9999999){

  TString path = ddir+flist+suffix;

  /// if is a root file
  bool temp = false;
  TFile *f = new TFile(path.Data());
  if(f && !f->IsZombie() ){
    cout << endl << "ROOT file "<<flist;
    ofstream outt;
    outt.open("temp.txt", ios::out | ios::trunc); //open input file list
    outt << flist << "   1";
    outt << endl;
    outt.close();
    flist="temp.txt";
    temp=true;
  }else{
    cout << endl << " NOT a ROOT file...";
  }

  cout << endl << "Tree --> "<<tree.Data();
  
  TChain *ch = new TChain(tree.Data());
  /// if is a list o files
  if(!temp)cout << endl << "Opening file list: "<<flist<<endl;
  ifstream innn;
  innn.open(flist.Data(), ios::in); //open input file list
  int nfile=0;
  while (1) {
    TString file;
    int flag = 1;
    innn >> file;
    ///    innn >> flag;
    if (!innn.good()) break;
    TString path = ddir+file+suffix;
    //         cout <<endl<< nfile <<" - "<< path;
    TFile *f = new TFile(path.Data());
    if(!f) continue;
    if(f->IsZombie() )continue;
    delete f;

    if( flag ){
      ch->Add(path.Data());
      // cout <<endl<< nfile <<" - "<< path;
      // cout << "  ++ "<<ch->GetEntries();
    }
    nfile++;
    if(nfile==max)break;
  };
  innn.close();
  cout << endl << endl;
  cout <<ch->GetEntries()<<" entries "<<endl;;
  return ch;
  
}

