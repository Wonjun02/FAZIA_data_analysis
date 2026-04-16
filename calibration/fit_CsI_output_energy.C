// fit_qxtx_CsI_output_energy.C
#include "TMath.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TString.h"
#include <iostream>

//output-energy function 정의 (x[0] = E)
double output_energy_function(double *x, double *par) {
    double E = x[0];
    if (E <= 0) return 0.0; 

    double a1 = par[0]; // Gain (Variable)
    double a2 = par[1]; // Quenching parameter
    double a3 = par[2]; // Low energy integration constant
    double a4 = par[3]; // Delta-ray contribution
    double A  = par[4]; // Mass number
    double Z  = par[5]; // Atomic number

    double term1 = a2 * A * Z * Z;
    if (term1 <= 0) return 0.0;

    double partA = E * (1.0 - (term1 / E) * log(1.0 + E / term1));
    double partB = a4 * term1 * log((E + term1) / (a3 + term1));

    return a1 * (partA + partB);
}

void fit_CsI_output_energy()
{
    gStyle->SetOptStat(0);

    const int n = 4;
    const int nGraph = 8;

    // X values (Energy -> MeV)
    double E[n] = {65.98, 60.62, 55.8, 50.82};
    double sigma_E[n] = {0.1964, 1.323, 0.1876, 1.437};

    // Y values (Channel -> adc)
    double Q1T2[n] = {172.5, 155.0, 143.3, 127.6};
    double Q4T4[n] = {208.3, 190.0, 176.2, 158.8};
    double Q1T1[n] = {168.3, 150.0, 144.5, 132.7};
    double Q1T3[n] = {187.8, 169.8, 152.4, 136.2};
    double Q2T3[n] = {173.2, 156.7, 148.6, 132.6};
    double Q2T4[n] = {137.4, 120.2, 115.0, 105.0};
    double Q3T2[n] = {140.6, 129.5, 119.1, 105.0};
    double Q4T1[n] = {162.4, 150.3, 142.7, 129.2};

    double sigma_Q1T2[n] = {0., 0., 0., 0.};
    double sigma_Q4T4[n] = {0., 0., 0., 0.};
    double sigma_Q1T1[n] = {0., 0., 0., 0.};
    double sigma_Q1T3[n] = {0., 0., 0., 0.};
    double sigma_Q2T3[n] = {0., 0., 0., 0.};
    double sigma_Q2T4[n] = {0., 0., 0., 0.};
    double sigma_Q3T2[n] = {0., 0., 0., 0.};
    double sigma_Q4T1[n] = {0., 0., 0., 0.};

    const char* names[nGraph] = {
        "Q1T2", "Q4T4", "Q1T1", "Q1T3",
        "Q2T3", "Q2T4", "Q3T2", "Q4T1"
    };

    double* adc_vals[nGraph] = {
        Q1T2, Q4T4, Q1T1, Q1T3,
        Q2T3, Q2T4, Q3T2, Q4T1
    };

    double* sigma_adc_vals[nGraph] = {
        sigma_Q1T2, sigma_Q4T4, sigma_Q1T1, sigma_Q1T3,
        sigma_Q2T3, sigma_Q2T4, sigma_Q3T2, sigma_Q4T1
    };

    int colors[nGraph] = {
        kBlue+1, kRed+1, kGreen+2, kMagenta+1,
        kOrange+7, kCyan+2, kBlack, kViolet
    };

    for (int i = 0; i < nGraph; i++) {
        TCanvas* c = new TCanvas(Form("c_%s", names[i]), names[i], 800, 650);
        c->SetGrid();
        c->SetTicks();

        TGraphErrors* gr = new TGraphErrors(n, E, adc_vals[i], sigma_E, sigma_adc_vals[i]);
        gr->SetName(Form("gr_%s", names[i]));
        gr->SetTitle("");

        gr->SetMarkerStyle(20);
        gr->SetMarkerColor(colors[i]);
        gr->SetLineColor(colors[i]);
        gr->SetMarkerSize(1.5);
        gr->SetLineWidth(2);

        gr->Draw("AP");

        gr->GetXaxis()->SetTitle("Energy (MeV)");
        gr->GetYaxis()->SetTitle("Output (adc)");

        double min_E = TMath::MinElement(n, E);
        double max_E = TMath::MaxElement(n, E);
        double fit_min = min_E - 5.0;
        double fit_max = max_E + 5.0;

        //gr->GetXaxis()->SetLimits(fit_min, fit_max);
        gr->GetXaxis()->SetLimits(0., 70.);
        
        double min_adc = TMath::MinElement(n, adc_vals[i]);
        double max_adc = TMath::MaxElement(n, adc_vals[i]);
        gr->GetYaxis()->SetRangeUser(min_adc - 20.0, max_adc + 20.0);
        gr->GetYaxis()->SetRangeUser(0., max_adc + 20.0);
        gr->GetXaxis()->SetRangeUser(0., 70.);
        gr->GetXaxis()->SetTitleSize(0.045);
        gr->GetYaxis()->SetTitleSize(0.045);
        gr->GetXaxis()->SetLabelSize(0.035);
        gr->GetYaxis()->SetLabelSize(0.035);
        gr->GetXaxis()->SetTitleOffset(1.0);
        gr->GetYaxis()->SetTitleOffset(1.15);

        TF1* fFit = new TF1(Form("fFit_%s", names[i]), output_energy_function, fit_min, fit_max, 6);
        //TF1* fFit = new TF1(Form("fFit_%s", names[i]), output_energy_function, 0., 70., 6);
        fFit->SetParNames("a1", "a2", "a3", "a4", "A", "Z");
        fFit->SetLineColor(colors[i]);
        fFit->SetLineWidth(2);

        // 파라미터 초기값 설정 및 고정
        fFit->SetParameter(0, 3.0);  // a1 (Gain - 피팅할 변수)
        fFit->FixParameter(1, 0.25); // a2
        fFit->FixParameter(2, 3.1); // a3
        fFit->FixParameter(3, 0.27); // a4
        fFit->FixParameter(4, 1.0);  // A (Proton)
        fFit->FixParameter(5, 1.0);  // Z (Proton)

        gr->Fit(fFit, "RMI");
        gr->Draw("P SAME");
        fFit->SetRange(0.,70.);
	fFit->Draw("SAME");
        // 모든 파라미터 값 가져오기
        double a1_val = fFit->GetParameter(0);
        double a2_val = fFit->GetParameter(1);
        double a3_val = fFit->GetParameter(2);
        double a4_val = fFit->GetParameter(3);
        double A_val  = fFit->GetParameter(4);
        double Z_val  = fFit->GetParameter(5);
        double chi2   = fFit->GetChisquare();
        int ndf       = fFit->GetNDF();

        TLegend* leg = new TLegend(0.12, 0.80, 0.35, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.04);
        leg->AddEntry(gr, Form("%s (Proton)", names[i]), "lp");
        leg->Draw();

        
        TLatex* tex = new TLatex();
        tex->SetNDC();
        tex->SetTextSize(0.035); 

        //tex->DrawLatex(0.12, 0.75, TString::Format("Fit: a_{1} = %.3f", a1_val).Data()); 
        tex->DrawLatex(0.12, 0.75, TString::Format("a_{1} = %.3f", a1_val).Data()); 
        //tex->DrawLatex(0.12, 0.70, TString::Format("Fix: a_{2}=%.2f, a_{3}=%.2f, a_{4}=%.2f", a2_val, a3_val, a4_val).Data()); 
        tex->DrawLatex(0.12, 0.70, TString::Format("a_{2}=%.2f, a_{3}=%.2f, a_{4}=%.2f", a2_val, a3_val, a4_val).Data()); 
        tex->DrawLatex(0.12, 0.65, TString::Format("A=%.0f, Z=%.0f", A_val, Z_val).Data()); 
        //tex->DrawLatex(0.12, 0.65, TString::Format("Particle: A=%.0f, Z=%.0f", A_val, Z_val).Data()); 
        tex->DrawLatex(0.12, 0.60, TString::Format("#chi^{2}/ndf = %.2f / %d", chi2, ndf).Data()); 

        c->Update();
        c->SaveAs(Form("fit_CsI_output_energy_%s.pdf", names[i]));


        std::cout << "=====================================" << std::endl;
        std::cout << "Channel: " << names[i] << std::endl;
        std::cout << "  a1 (Gain) = " << a1_val << std::endl;
        std::cout << "  a2 = " << a2_val << ", a3 = " << a3_val << ", a4 = " << a4_val << std::endl;
        std::cout << "  A = " << A_val << ", Z = " << Z_val << std::endl;
        std::cout << "  Chi2/ndf  = " << chi2 << " / " << ndf << std::endl;
        std::cout << "=====================================\n" << std::endl;
    }
}
