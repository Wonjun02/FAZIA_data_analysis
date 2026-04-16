void draw_PID()
{
  gStyle->SetOptStat(0);

  int En = 60;
  //int Energies[] = {70, 60, 50, 40};
  int qSel = 1;
  int tSel = 2;

  TF1* fTG[4][4];
  TH1D* histProj_x[4][4];
  TH1D* histProj_y[4][4];
  TH1D* histDiff[4][4];
  auto file = new TFile( Form("histogram_corr_%dMeV.root", En), "READ" );

  for( int i=0; i<4; i++ )
    for( int j=0; j<4; j++ )
    {
      qSel = i+1;
      tSel = j+1;
      auto h2 = (TH2D*) file->Get( Form("histCorr_q%d_t%d", qSel, tSel) );

      // 6) Draw 2D and pick ridge points
      TCanvas* c = new TCanvas("c2D", Form("Q_%d, T_%d", qSel, tSel), 800, 600);
      gPad->SetRightMargin(0.12);
      gPad->SetLogz();
      h2->Draw("colz");
      c->Update();

      const int nPick = 12;
      float vx[nPick], vy[nPick];

      printf("\n>>> Click %d points on the proton ridge (LEFT mouse)\n", nPick);
      printf(">>> Release button to register. A red marker will appear.\n");

      for (int i = 0; i < nPick; ) {
        gSystem->ProcessEvents();
        gPad->WaitPrimitive();

        int ev = gPad->GetEvent();
        if (ev != kButton1Up) continue;

        int px = gPad->GetEventX();
        int py = gPad->GetEventY();

        double xx = gPad->AbsPixeltoX(px);
        double yy = gPad->AbsPixeltoY(py);
        xx = gPad->PadtoX(xx);
        yy = gPad->PadtoY(yy);

        vx[i] = (float)xx;
        vy[i] = (float)yy;

        TMarker* m = new TMarker(vx[i], vy[i], 29);
        m->SetMarkerColor(kRed);
        m->SetMarkerSize(1.4);
        m->Draw("SAME");

        printf("  point %2d : x = %.2f , y = %.2f\n", i, vx[i], vy[i]);

        gPad->Modified();
        gPad->Update();
        i++;
      }



      double xmin = 0;
      double xmax = 1000;
      fTG[i][j] = new TF1( Form("fTG_q%d_t%d", qSel, tSel),
          "pow( pow([0]*x,[1]+1) + pow([2],[1]+1)*pow([3],2)*pow([4],[1]) , 1./([1]+1) ) - [0]*x + [5]",
          xmin, xmax);

      fTG[i][j]->SetParNames("g","mu","lambda","Z","A","pdy");
      fTG[i][j]->SetParameters(1.0, 0.6, 160.0, 1.0, 1.0, 0.0);  // lambda는 ΔE(0) 스케일이라 40 근방 추천

      // Considering only proton
      fTG[i][j]->FixParameter(3, 1.0); // Z
      fTG[i][j]->FixParameter(4, 1.0); // A

      TGraph* grFit = new TGraph(nPick, vx, vy);
      grFit->Fit(fTG[i][j], "R");

      float g      = fTG[i][j]->GetParameter(0);
      float mu     = fTG[i][j]->GetParameter(1);
      float lambda = fTG[i][j]->GetParameter(2);
      float pdy    = fTG[i][j]->GetParameter(5);

      fTG[i][j]->SetLineColor(kRed);
      fTG[i][j]->SetLineWidth(3);
      fTG[i][j]->Draw("SAME");
      c->Update();

      /*
         std::cout << "\n=== Fit result ===\n";
         std::cout << "g      = " << g << "\n";
         std::cout << "mu     = " << mu << "\n";
         std::cout << "lambda = " << lambda << "\n";
         std::cout << "pdy    = " << pdy << "\n";
       */



      histProj_x[i][j] = new TH1D(Form("histProj_x_q%d_t%d", qSel, tSel),Form("histProj_x_q%d_t%d", qSel, tSel), 250, 50, 300);
      histProj_y[i][j] = new TH1D(Form("histProj_y_q%d_t%d", qSel, tSel),Form("histProj_y_q%d_t%d", qSel, tSel), 100, 0, 100);
      histDiff[i][j] = new TH1D(Form("histDiff%d_t%d", qSel, tSel),Form("histDiff%d_t%d", qSel, tSel), 150, -20, 50);

      // misc
      //auto h2Clone = new TH2D( "h2Clone", "", h2->GetNbinsX(), h2->GetXaxis()->GetXmin(), h2->GetXaxis()->GetXmax(), 140, -20, 50 );

      int nx = h2->GetNbinsX();
      int ny = h2->GetNbinsY();
      double compareVal = 1.;

      for( int ix=0; ix<nx; ix++ )
      {
        double xx = h2->GetXaxis()->GetBinCenter( ix+1 );
        double fitVal = fTG[i][j]->Eval( xx );
        for( int iy=0; iy<ny; iy++ )
        {
          double yy = h2->GetYaxis()->GetBinCenter( iy+1 );
          double counts = h2->GetBinContent( ix+1, iy+1 );

          double diff = yy - fitVal;

          // misc
          //h2Clone->Fill( xx, diff, counts );

          histDiff[i][j]->Fill( diff, counts );

          if( std::abs(diff) < compareVal )
          {
            // Finding x binning of histProj
            int idx_x = histProj_x[i][j]->FindBin( xx );

            // Then accumulate the histogram
            //x-axis
	    double val_x = histProj_x[i][j]->GetBinContent( idx_x ) + counts;
            histProj_x[i][j]->SetBinContent( idx_x, val_x );
            //y-axis
	    int idx_y = histProj_y[i][j]->FindBin(yy);
	    double val_y = histProj_y[i][j] -> GetBinContent( idx_y ) + counts;
	    histProj_y[i][j] -> SetBinContent( idx_y, val_y);

	  }
        }
      }

      //TCanvas* can = new TCanvas("c","c",700,500);
      //histProj->Draw();

      // misc
      //new TCanvas();
      //h2Clone->Draw( "colz" );
    }


    auto outputFile = new TFile( Form("data_%dMeV.root", En), "RECREATE" );
  for( int i=0; i<4; i++ )
    for( int j=0; j<4; j++ )
    {
      qSel = i+1;
      tSel = j+1;
      fTG[i][j]->Write();
      auto h2 = (TH2D*) file->Get( Form("histCorr_q%d_t%d", qSel, tSel) );
      h2->Write();
      histProj_x[i][j]->Write();
      histProj_y[i][j]->Write();
      histDiff[i][j]->Write();
    }
}
