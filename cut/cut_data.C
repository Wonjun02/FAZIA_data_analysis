void cut_data()
{
  gStyle->SetOptStat(0);
  int En = 70; 

  TFile* fFit = new TFile(Form("data_%dMeV.root", En), "READ");
  TF1* fTG[4][4];
  for(int i=0; i<4; i++) {
    for(int j=0; j<4; j++) {
      fTG[i][j] = (TF1*)fFit->Get(Form("fTG_q%d_t%d", i+1, j+1));
    }
  }

  TChain chain("faziatree");
  //int startRun = 43, endRun = 381; // 50MeV_1 /day2
  //int startRun = 382, endRun = 545; // 50MeV_2 /day2
  //int startRun = 546, endRun = 589; // 50MeV_3 /day2
  //int startRun = 43, endRun = 589; // total 50MeV /day2
  //int startRun = 590, endRun = 1095; // 40MeV /day2
  int startRun = 340, endRun = 1032; // 70MeV /day1
  //int startRun = 1033, endRun = 1699; // 60MeV /day1

  for (int run=startRun; run<=endRun; run++) 
    chain.Add(Form("/media/wonjun/T7/converted_data/run_%d.root", run));

  chain.SetBranchStatus("*", 0);
  chain.SetBranchStatus("Mtot",   1);
  chain.SetBranchStatus("fQua",   1);
  chain.SetBranchStatus("fTel",   1);
  chain.SetBranchStatus("sQ3max", 1);
  chain.SetBranchStatus("sQ2max", 1);

  UShort_t mtot, fQua[100], fTel[100];
  float sQ3max[100], sQ2max[100];
  chain.SetBranchAddress("Mtot",   &mtot);
  chain.SetBranchAddress("fQua",   fQua);
  chain.SetBranchAddress("fTel",   fTel);
  chain.SetBranchAddress("sQ3max", sQ3max);
  chain.SetBranchAddress("sQ2max", sQ2max);

  TH1D* histDiff_y[4][4];
  for(int i=0; i<4; i++) {
    for(int j=0; j<4; j++) {
      
      histDiff_y[i][j] = new TH1D(Form("histDiff_y_q%d_t%d",i+1,j+1), 
                                  Form("histDiff_y_q%d_t%d",i+1,j+1), 
                                  200, -50, 50);
    }
  }

  Long64_t nEntries = chain.GetEntries();
  for (Long64_t iEntry = 0; iEntry < nEntries; iEntry++) {
    chain.GetEntry(iEntry);
    int nHit = (mtot > 100) ? 100 : (int)mtot;

    for (int i = 0; i < nHit; i++) {
      int q = fQua[i] - 1;
      int t = fTel[i] - 1;
      if (q < 0 || q >= 4 || t < 0 || t >= 4) continue;
      if (!fTG[q][t]) continue; 

      float x_data = sQ3max[i]; 
      float y_data = sQ2max[i]; 

      double fit_y = fTG[q][t]->Eval(x_data);
      double dy = y_data - fit_y;
      histDiff_y[q][t]->Fill(dy);
    }
    if(iEntry % 100000 == 0) cout << "Processing... " << (100.0*iEntry/nEntries) << "%" << endl;
  }

  auto outputFile = new TFile(Form("Diff_y_%dMeV.root", En), "RECREATE");
  for(int i=0; i<4; i++) {
    for(int j=0; j<4; j++) {
      if(fTG[i][j]) fTG[i][j]->Write();
      if(histDiff_y[i][j]) histDiff_y[i][j]->Write();
    }
  }
  outputFile->Close(); 
  cout << " Results saved in cut_data_" << En << "MeV.root" << endl;
}
