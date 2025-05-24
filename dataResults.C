#include <iostream>
#include <chrono>
#include <string>
#include "dataHelperFunctions.h"
using namespace std;
using namespace std::chrono;

template<typename T>
void PrintTime(T Start, std::string String){
  auto Stop = chrono::high_resolution_clock::now();
  auto Duration = duration_cast<microseconds>(Stop-Start);
  cout<<String<<float(Duration.count())/float(1000000)<<" seconds"<<endl;
}

int GetBinFromTHnSparse(THnSparse* sparse, int AxisNo, double val){return sparse->GetAxis(AxisNo)->FindBin(val);}
double GetBinWidthOfAxis(THnSparse* sparse, int AxisNo, int BinNo){return sparse->GetAxis(AxisNo)->GetBinWidth(BinNo);}
double GetBinCenterOfAxis(THnSparse* sparse,int AxisNo, int BinNo){ return sparse->GetAxis(AxisNo)->GetBinCenter(BinNo);}
double GetBinLowEdgeOfAxis(THnSparse* sparse,int AxisNo, int BinNo){ return sparse->GetAxis(AxisNo)->GetBinLowEdge(BinNo);}
double GetBinUpEdgeOfAxis (THnSparse* sparse,int AxisNo, int BinNo){ return (sparse->GetAxis(AxisNo)->GetBinLowEdge(BinNo)+sparse->GetAxis(AxisNo)->GetBinWidth(BinNo));}


int GetRndmRandomNumber(TRandom3 &randomGenerator, int Length){
  int randomValue = randomGenerator.Rndm()*Length; // Uniform number between 0 and 1
  if (randomValue == Length){ randomValue = Length -1;}
  return randomValue;
}
// Required Results printing functions***************************************
template<typename T, typename U, typename V>
void Get_NM_FromPowerSum(T PowerSum, U iEvtEntries, V &Moment, int nClass, int nOrder){
  cout<<"Obtaining Normal Moments ******"<<endl;
  for(int iClass=0; iClass<nClass; iClass++){
    for(int iOrder=0; iOrder<=nOrder; iOrder++){
      if(iEvtEntries[iClass] != 0){Moment[iClass][iOrder]=double(PowerSum[iClass][iOrder])/static_cast<double>(iEvtEntries[iClass]);}
      else{Moment[iClass][iOrder]=0;
      }
    }
  }
}

template<typename T, typename V>
void Get_C_From_NM(T NM_nClass, V &C_NM_nClass, int nClass){
  cout<<"Finding Cumulants from Normal Moments**********"<<endl;
  double X0, X1, X2, X3, X4, X5, X6;
  for(int iClass=0;iClass<nClass; iClass++){
    X0 =NM_nClass[iClass][0];
    X1 =NM_nClass[iClass][1];
    X2 =NM_nClass[iClass][2];
    X3 =NM_nClass[iClass][3];
    X4 =NM_nClass[iClass][4];
    X5 =NM_nClass[iClass][5];
    X6 =NM_nClass[iClass][6];

    C_NM_nClass[iClass][0] = X0;
    C_NM_nClass[iClass][1] = X1;
    C_NM_nClass[iClass][2] = X2 - TMath::Power(X1,2);
    C_NM_nClass[iClass][3] = X3 - 3.0*X2*X1 + 2.0*TMath::Power(X1,3);
    C_NM_nClass[iClass][4] = X4 - 4.0*X3*X1 - 3.0*TMath::Power(X2,2) + 12.0*X2*TMath::Power(X1,2) - 6.0*TMath::Power(X1,4);
    C_NM_nClass[iClass][5] = X5 - 5.0*X4*X1 -10.0*X3*X2 + 20*X3*TMath::Power(X1,2) + 30.0*TMath::Power(X2,2)*X1 - 60.0*X2*TMath::Power(X1,3) + 24.0*TMath::Power(X1,5);
    C_NM_nClass[iClass][6] = X6 - 6.0*X5*X1 -15.0*X4*X2 + 30*X4*TMath::Power(X1,2) - 10.0*TMath::Power(X3,2)    +120.0*X3*X2*X1 -120.0*X3*TMath::Power(X1,3) + 30.0*TMath::Power(X2,3) -270.0*TMath::Power(X2,2)*TMath::Power(X1,2)+360*X2*TMath::Power(X1,4)-120.0*TMath::Power(X1,6);
  }
}

template<typename T, typename U>
void Get_FC_From_C(T C_CM_nClass, U &FC_C_nClass, int nClass ){
  cout << "Finding Factorial Cumulants from Cumulants"<<endl;
  double C0, C1, C2, C3, C4, C5, C6;
  for(int iClass=0; iClass< nClass; iClass++){
      C0 = C_CM_nClass[iClass][0];
      C1 = C_CM_nClass[iClass][1];
      C2 = C_CM_nClass[iClass][2];
      C3 = C_CM_nClass[iClass][3];
      C4 = C_CM_nClass[iClass][4];
      C5 = C_CM_nClass[iClass][5];
      C6 = C_CM_nClass[iClass][6];

      FC_C_nClass[iClass][0] = C0;
      FC_C_nClass[iClass][1] = C1;
      FC_C_nClass[iClass][2] = C2 - C1;
      FC_C_nClass[iClass][3] = C3 -  3*C2 +  2*C1;
      FC_C_nClass[iClass][4] = C4 -  6*C3 + 11*C2 -   6*C1;
      FC_C_nClass[iClass][5] = C5 - 10*C4 + 35*C3 -  50*C2 +  24*C1;
      FC_C_nClass[iClass][6] = C6 - 15*C5 + 85*C4 - 225*C3 + 274*C2 - 120*C1;
  }
}
//other functions*******************************

void dataResults(){
    TFile *infile = new TFile("AnalysisResults.root");
    infile->ls();

    THnSparseD *sparse1 = (THnSparseD*)infile->Get("nch-cumulants-id/sparse1");
    THnSparseD *sparse2 = (THnSparseD*)infile->Get("nch-cumulants-id/sparse2");

    sparse1->Print();
    sparse2->Print();

    TCanvas *c = new TCanvas();
    auto hist = sparse1->Projection(8);
    hist->SetTitle("centFT0M");
    hist->Draw();

    const int nClass = 5;
    const int nSubSample = 50;
    const int nOrder = 6;

    vector<double> FT0M_classedges = {5,20,40,60,80,95};
    vector<double> FT0M_classLow = {5,20,40,60,80};
    vector<double> FT0M_classUp = {20,40,60,80,95};
    int modeType  = kFloatTypeAxis;
    // cout<<"fourth class edge is :"<<FT0M_classedges[3]<<endl;

    if(FT0M_classedges.size() != (nClass+1) ){ cout<<"ERROR::ERROR:: class and classedge mismatch at "<<(nClass+1)<<endl;return;}

    //treat each bin separately
    //Find the number of small bins available.
    cout<<"HistName = "<<hist->GetName()<<endl;
    cout<<"nBins    = "<<hist->GetNbinsX()<<endl;
    cout<<"xAxisLow = "<<hist->GetBinLowEdge(1)<<endl;
    cout<<"xAxisUp  = "<<hist->GetBinLowEdge(hist->GetNbinsX())+hist->GetBinWidth(hist->GetNbinsX())<<endl;
   
    //To Reduce the memory size find the first and the last filled bin 
    int firstFilledBin = getFirstFilledBin(hist);
    int lastFilledBin  = getLastFilledBin(hist); 
  
    cout<<"firstFilledBin = "<<firstFilledBin<<endl;
    cout<<"lastFilledBin  = "<<lastFilledBin<<endl;
    cout<<"lowEdgeOfLowestBin  = "<<hist->GetBinLowEdge(firstFilledBin)<<endl;
    cout<<"lowEdgeOfHighestBin = "<<hist->GetBinLowEdge(lastFilledBin)<<endl;
    const int nSmallClass =hist->GetBinLowEdge(lastFilledBin) - hist->GetBinLowEdge(firstFilledBin) + 1;//
    //const int nSmallClass = 10;
    cout<<"nSmallClass = "<<nSmallClass<<endl;

    vector<double> smallClassLow(nSmallClass);
    vector<double> smallClassUp(nSmallClass);
    vector<int> nClassIDX_smallClass(nSmallClass);
    vector<int> iBinLow_smallClass(nSmallClass);
    vector<int> iBinUp_smallClass(nSmallClass);

    cout<<"Int Mode"<<endl;
    getSmallBinProperEdges(hist, nSmallClass, smallClassLow, smallClassUp, iBinLow_smallClass, iBinUp_smallClass, nClassIDX_smallClass, nClass, FT0M_classLow,FT0M_classUp, firstFilledBin, lastFilledBin, kIntTypeAxis, "centFT0M");
    modeType  = kFloatTypeAxis;
    cout<<"float Mode"<<endl;
    getSmallBinProperEdges(hist, nSmallClass, smallClassLow, smallClassUp, iBinLow_smallClass, iBinUp_smallClass, nClassIDX_smallClass, nClass, FT0M_classLow,FT0M_classUp, firstFilledBin, lastFilledBin, modeType, "centFT0M");
  
    //Print integral of each big class:
    for(int i=0; i< nClass; i++){
    cout<<"class:"<<i<<" :: bins = ["<<hist->GetXaxis()->FindBin(FT0M_classLow[i])<<","<<hist->GetXaxis()->FindBin(FT0M_classUp[i])-1<<"]"<<" :: integral = "<<hist->Integral(hist->FindBin(FT0M_classLow[i]),hist->FindBin(FT0M_classUp[i])-1)<<endl;}

    cout<<"hist total Entries = "<<hist->GetEntries()<<endl;
    //Find Entries Per SubSample
    int64_t totalSubSampleEntries = 0;
    int64_t EntriesPerSubSample = FindEntriesPerSubSample(hist, nSubSample);
    cout<<"Entries Per SubSample = "<<EntriesPerSubSample<<endl;
  
    int64_t SubSampleEntryCounter[nSubSample]; 
    for (int i = 0 ; i < nSubSample; i++){ SubSampleEntryCounter[i] = 0;}

    //How to assign random SubSample number to one event
    unsigned int seed = 123456789; // Seed value; you can change it to any fixed integer
    TRandom3 randomGenerator(seed);

    // INITIALIZING TO ZERO USING VECTOR INSTEAD OF MANNUAL GROUPING.....................
    //Using Vectors
    //WITHOUT CBWC****************************************************
    std::vector<std::vector<double>> dN_PowerSum_nClass_NoCBWC_vec(nClass+2, std::vector<double>(nOrder + 1,0));
    std::vector<int64_t>              iEvtEntries_nClass_NoCBWC_vec(nClass+2,0);

    //For CBWC *******************************************************
    std::vector<std::vector<double>> dN_PowerSum_nSmallClass_vec(nSmallClass, std::vector<double>(nOrder + 1,0)); //nth power of req variable; dN can be Np, Nk etc
    std::vector<int64_t> iEvtEntries_nSmallClass_vec(nSmallClass,0); //no of events of that class 

    //FOR STATISTICAL UNCERTAINITY FOR CBWC AND WITHOUT CBWC**********
    std::vector<std::vector<std::vector<double>>> dN_PowerSum_nSubSample_nClass_NoCBWC_vec(nSubSample,std::vector<std::vector<double>>(nClass+2, std::vector<double>(nOrder + 1,0)));
    std::vector<std::vector<int64_t>>             iEvtEntries_nSubSample_nClass_NoCBWC_vec(nSubSample, std::vector<int64_t>(nClass+2,0));
  
    std::vector<std::vector<std::vector<double>>> dN_PowerSum_nSubSample_nSmallClass_CBWC_vec(nSubSample,std::vector<std::vector<double>>(nSmallClass, std::vector<double>(nOrder + 1,0)));
    std::vector<std::vector<int64_t>>             iEvtEntries_nSubSample_nSmallClass_CBWC_vec(nSubSample, std::vector<int64_t>(nSmallClass,0));
    //ZERO INITIALIZING ENDS********************************************

    double pow_dN_iOrder[nOrder + 1] ;

    int BinCountCheck0 = 200000;
    int BinCountCheck1 = 20000;
    int BinCountCheck2 = 1000;
    int BinCountCheck3 = 100;
    int BinCountCheck4 = 10;
    const auto& sparse = sparse1; // change sprse here******the goto Xvar[]**** then dN
    const int nDim = sparse->GetNdimensions();
    int THnBinValue[nDim];
    Double_t Xvar[nDim];
    Long64_t nBinsSparse = sparse->GetNbins();
    Long64_t iBinEntries = 0;
    int64_t eventsCounted = 0;
    cout<<" nDim = "<<nDim<<" :: "<<"nBinsSparse = "<<nBinsSparse<<endl;

    std::vector<TH1D*> h(nClass+2, nullptr);
    for (int i = 0 ; i < nClass+2 ; i++){
      h[i] = new TH1D(Form("h_%d",i), Form("h_%d",i), 103, -1, 102);
    }
    vector<TH1D*> h1D(nDim, nullptr);
    vector<vector<TH1D*>> h1D_nClass(nDim, vector<TH1D*>(nClass+2, nullptr));

    auto Start1 = chrono::high_resolution_clock::now();
    auto Start2 = chrono::high_resolution_clock::now();

    int    dN       ;
    int    nPi      ;
    int    nAPi     ;
    int    Np       ; 
    int    Nm       ; 
    int    nPr      ; 
    int    nAPr     ; 
    int    nKa      ; 
    int    nAKa     ; 
    int    Nch      ; 
    double centFT0M ;

    int iClassIDX    = -1;
    int iSmallBinIDX = -1;
    int iSubSample   = -1;
    //BIN by BIN LOOP STARTS HERE
    for(int iBin = 0; iBin < nBinsSparse; iBin++){
      if(iBin < BinCountCheck0){
        if(iBin < BinCountCheck1){
          if(iBin < BinCountCheck2){
            if(iBin < BinCountCheck3){
              if(iBin % BinCountCheck4 == 0){
                cout<<endl;
                PrintTime(Start2, Form("    %s :: BinReading :: BinTime     :: iBin = %d :: eventsCounted = (%f)%ld :: ",sparse->GetName(),iBin, static_cast<float>(100. * double(eventsCounted)/(1.94E+10)),eventsCounted));
                PrintTime(Start1, Form("    %s :: BinReading :: ElapsedTime :: iBin = %d :: eventsCounted = (%f)%ld :: ",sparse->GetName(),iBin, static_cast<float>(100. * double(eventsCounted)/(1.94E+10)),eventsCounted));
                Start2 = chrono::high_resolution_clock::now();
              }
            } else if (iBin % BinCountCheck3 == 0){
                cout<<endl;
                PrintTime(Start2, Form("    %s :: BinReading :: BinTime     :: iBin = %d :: eventsCounted = (%f)%ld :: ",sparse->GetName(),iBin, static_cast<float>(100. * double(eventsCounted)/(1.94E+10)),eventsCounted));
                PrintTime(Start1, Form("    %s :: BinReading :: ElapsedTime :: iBin = %d :: eventsCounted = (%f)%ld :: ",sparse->GetName(),iBin, static_cast<float>(100. * double(eventsCounted)/(1.94E+10)),eventsCounted));
                Start2 = chrono::high_resolution_clock::now();
            }
          }
          if(iBin % BinCountCheck2 == 0){
          cout<<endl;
          PrintTime(Start2, Form("    %s :: BinReading :: BinTime     :: iBin = %d :: eventsCounted = (%f)%ld :: ",sparse->GetName(),iBin, static_cast<float>(100. * double(eventsCounted)/(1.94E+10)),eventsCounted));
          PrintTime(Start1, Form("    %s :: BinReading :: ElapsedTime :: iBin = %d :: eventsCounted = (%f)%ld :: ",sparse->GetName(),iBin, static_cast<float>(100. * double(eventsCounted)/(1.94E+10)),eventsCounted));
          Start2 = chrono::high_resolution_clock::now();
          }
        }
        else if(iBin % BinCountCheck1 == 0){
          cout<<endl;
          PrintTime(Start2, Form("    %s :: BinReading :: BinTime     :: iBin = %d :: eventsCounted = (%f)%ld :: ",sparse->GetName(),iBin, static_cast<float>(100. * double(eventsCounted)/(1.94E+10)),eventsCounted));
          PrintTime(Start1, Form("    %s :: BinReading :: ElapsedTime :: iBin = %d :: eventsCounted = (%f)%ld :: ",sparse->GetName(),iBin, static_cast<float>(100. * double(eventsCounted)/(1.94E+10)),eventsCounted));
          Start2 = chrono::high_resolution_clock::now();
        }
      }
      else if(iBin % BinCountCheck0 == 0){
        cout<<endl;
        PrintTime(Start2, Form("    %s :: BinReading :: BinTime     :: iBin = %d :: eventsCounted = (%f)%ld :: ",sparse->GetName(),iBin, static_cast<float>(100. * double(eventsCounted)/(1.94E+10)),eventsCounted));
        PrintTime(Start1, Form("    %s :: BinReading :: ElapsedTime :: iBin = %d :: eventsCounted = (%f)%ld :: ",sparse->GetName(),iBin, static_cast<float>(100. * double(eventsCounted)/(1.94E+10)),eventsCounted));
        Start2 = chrono::high_resolution_clock::now();
      }

      sparse->GetBinContent(static_cast<int>(iBin),THnBinValue);
      for(int iAxis = 0 ; iAxis < nDim ; iAxis++){
        Xvar[iAxis] = GetBinCenterOfAxis(sparse,iAxis,THnBinValue[iAxis]);
      }
      iBinEntries = sparse->GetBinContent(static_cast<int>(iBin));
      eventsCounted += iBinEntries ;
      if(iBin < 10) {cout<<"    iBin = "<<iBin<<" :: Entries = "<<iBinEntries<<" :: eventsCounted = ("<<static_cast<float>(100. * double(eventsCounted)/(1.94E+10))<<")"<<eventsCounted<<endl;}
      // Direct Weight Filling
      // MasterHistogram->Fill(Xvar,iBinEntries);

      dN       = Xvar[0];
      Np       = Xvar[1];
      Nm       = Xvar[2];
      nPr      = Xvar[3];
      nAPr     = Xvar[4];
      //nPi      = Xvar[3];
      //nAPi     = Xvar[4];
      nKa      = Xvar[5];
      nAKa     = Xvar[6];
      Nch      = Xvar[7];
      centFT0M = Xvar[8];

      // cout<<"dN = "<<dN<<endl;
      iClassIDX    = -1;
      iSmallBinIDX = -1;
      dN  = Np-Nm;
      int dNP = nPr-nAPr;
      int dNK = nKa-nAKa;
      // int dNPI= nPi-nAPi;

      for(int i = 0 ; i < nClass ; i++){
        if( FT0M_classLow[i] <= centFT0M &&  centFT0M < FT0M_classUp[i]) {
          iClassIDX = i+1 ; break;
        }
        if(centFT0M < FT0M_classLow[0]      ) {iClassIDX = 0; break;}
        if(centFT0M > FT0M_classUp[nClass-1]) {iClassIDX = nClass+1; break;}
      }

      if(iClassIDX == -1) {cout<<"ERROR iClassIDX == -1 :: nClass = "<<nClass<<" :: centFT0M = "<<centFT0M<<endl;}

      h[iClassIDX]->Fill(centFT0M,iBinEntries);

      for(int i = 0 ; i < nSmallClass ; i++){
        if( smallClassLow[i] <= centFT0M &&  centFT0M < smallClassUp[i]) { iSmallBinIDX = i ;  break;}
      }
      if( centFT0M < 0 ){ cout<<"ERROR centFT0M < 0 :: iBin = "<<iBin<<endl; }

      if(iSmallBinIDX == -1) {cout<<"ERROR iSmallBinIDX == -1 :: nClass = "<<nClass<<" :: centFT0M = "<<centFT0M<<" :: iBin = "<<iBin<<endl;}

      if(nClassIDX_smallClass[iSmallBinIDX] != iClassIDX) { cout<<"ERROR :: nClassIDX_smallClass["<<iSmallBinIDX<<"] != iClassIDX :: "<<nClassIDX_smallClass[iSmallBinIDX]<<" != "<<iClassIDX<<
        " :: "<<centFT0M<<
        
        endl;}

       for(int iOrder = 0; iOrder<= nOrder; iOrder++){
        //No CBWC PowerSum
        dN_PowerSum_nClass_NoCBWC_vec[iClassIDX][iOrder] += TMath::Power(dN,iOrder)*iBinEntries;
      }
       //No CBWC PowerSum
      iEvtEntries_nClass_NoCBWC_vec[iClassIDX] += iBinEntries;
      
      for(int iOrder = 0; iOrder<= nOrder; iOrder++){
        //for CBWC PowerSum
        dN_PowerSum_nSmallClass_vec[iSmallBinIDX][iOrder] += TMath::Power(dN,iOrder)*iBinEntries;        
      }
       //for CBWC PowerSum
      iEvtEntries_nSmallClass_vec[iSmallBinIDX] += iBinEntries;

      for (int iOrder = 0; iOrder <= nOrder; iOrder++) {
        pow_dN_iOrder[iOrder] = TMath::Power(dN, iOrder); // debug the case of power(0,0), it gives the value = 1.0 by default
      }

      //  Get Momentum Sum in each Small Bin of each SubSample
      for (int iEntry = 0 ; iEntry < iBinEntries; iEntry++){
        iSubSample = GetRndmRandomNumber(randomGenerator, nSubSample);
        while ( SubSampleEntryCounter[iSubSample] >= EntriesPerSubSample ){ 
          iSubSample = GetRndmRandomNumber(randomGenerator, nSubSample); 
        }
        SubSampleEntryCounter[iSubSample]++;
        if(SubSampleEntryCounter[iSubSample] > EntriesPerSubSample){ cout<<"Problem :"<<iSubSample<<" :: ";}

        for (int iOrder = 0; iOrder <= nOrder; iOrder++) {
          //__Subsampling_NoMBWC - Start __________________//Fill large class power sums and entries for subsample____________________________________
          dN_PowerSum_nSubSample_nClass_NoCBWC_vec[iSubSample][iClassIDX][iOrder] += pow_dN_iOrder[iOrder];
          //__Subsampling_NoMBWC - End  ______________________________________________________
          //__Subsampling_MBWC - Start _____________________// // Fill small class power sums and entries for subsampleS
          dN_PowerSum_nSubSample_nSmallClass_CBWC_vec[iSubSample][iSmallBinIDX][iOrder] += pow_dN_iOrder[iOrder];
          //__Subsampling_MBWC - End ______________________________________________________
        }
        iEvtEntries_nSubSample_nClass_NoCBWC_vec[iSubSample][iClassIDX]++;
        iEvtEntries_nSubSample_nSmallClass_CBWC_vec[iSubSample][iSmallBinIDX]++;
        //___________________Debug Histogram______________
        //Histogram filling consumes a lot of time 
        // int iAxis = 8;
        //   h1D[iAxis]                  ->Fill(Xvar[iAxis]);
        //   h1D_nClass[iAxis][classIDX] ->Fill(Xvar[iAxis]);
      }
    }//BIN BY BIN LOOP ENDS ************************* 
    cout<<"CHECK1 :upto here passed::"<<endl;
    
    TCanvas *c10 = new TCanvas("c10", "c10");
    c10->Divide(3,3);
    for(int i = 0 ; i < nClass+2; i++){
      c10->cd(i+1);
      h[i]->Draw("hist");
    }
/*
    //PrintVariables(dN_PowerSum_nClass_NoCBWC_vec, iEvtEntries_nClass_NoCBWC_vec   , nClass, nOrder, "Printing PowerSum n class no CBWC", "iClass");
    //PrintVariables(dN_PowerSum_nSmallClass_vec, iEvtEntries_nSmallClass_vec  , nSmallClass, nOrder, "Printing PowerSum nsmallClass", "iSmallClass");

    //Calculations WITHOUT CBWC ****************************************

    std::vector<std::vector<double>> dN_NM_nClass_NoCBWC_vec(nClass, std::vector<double>(nOrder + 1,0));
    Get_NM_FromPowerSum(dN_PowerSum_nClass_NoCBWC_vec, iEvtEntries_nClass_NoCBWC_vec, dN_NM_nClass_NoCBWC_vec, nClass, nOrder);
    //PrintVariables(dN_NM_nClass_NoCBWC_vec    , iEvtEntries_nClass_NoCBWC_vec    , nClass, nOrder, "Printing PowerSum by vec", "iClass");

    std::vector<std::vector<double>> dN_C_NM_nClass_NoCBWC_vec(nClass, std::vector<double>(nOrder+1,0));
    Get_C_From_NM(dN_NM_nClass_NoCBWC_vec, dN_C_NM_nClass_NoCBWC_vec, nClass);
    //PrintVariables(dN_C_NM_nClass_NoCBWC_vec, iEvtEntries_nClass_NoCBWC_vec, nClass, nOrder, "Printing cumulants NM by vec", "iClass");

    std::vector<std::vector<double>> dN_FC_C_nClass_NoCBWC_vec(nClass, std::vector<double>(nOrder+1,0));
    Get_FC_From_C(dN_C_NM_nClass_NoCBWC_vec, dN_FC_C_nClass_NoCBWC_vec, nClass);
    //PrintVariables(dN_FC_C_nClass_NoCBWC_vec, iEvtEntries_nClass_NoCBWC_vec, nClass, nOrder, "Printing factorial cumulants from c by vec", "iClass");

    //calculations WITH CBWC *********************************************
    std::vector<std::vector<double>> dN_NM_nSmallClass_vec(nSmallClass, std::vector<double>(nOrder + 1,0));
    Get_NM_FromPowerSum(dN_PowerSum_nSmallClass_vec, iEvtEntries_nSmallClass_vec, dN_NM_nSmallClass_vec, nSmallClass, nOrder);

    std::vector<std::vector<double>> dN_C_NM_nSmallClass_vec(nSmallClass, std::vector<double>(nOrder+1,0));
    Get_C_From_NM(dN_NM_nSmallClass_vec, dN_C_NM_nSmallClass_vec, nSmallClass);

    std::vector<std::vector<double>> dN_FC_C_nSmallClass_vec(nSmallClass, std::vector<double>(nOrder+1,0));
    Get_FC_From_C(dN_C_NM_nSmallClass_vec, dN_FC_C_nSmallClass_vec, nSmallClass);

    std::vector<std::vector<double>> dN_CBWC_FC_vec(nClass, std::vector<double>(nOrder+1,0));
    Get_CBWC_C(dN_CBWC_FC_vec, dN_FC_C_nSmallClass_vec, nClass, nOrder, nSmallClass,FT0M_classLow,FT0M_classUp,smallClassLow,smallClassUp, iEvtEntries_nSmallClass_vec);
    //PrintVariables(dN_CBWC_FC_vec, iEvtEntries_nSmallClass_vec, nClass, nOrder, "Printing cbwc factorial cumulants C by vec", "iClass");
    

    TGraph *gr[7];
    double Xval[nClass] = {10.0,30.0,50.0,70.0,90.0};
    double Yval[nClass];
    for(int i=0; i<7; i++){
      for(int j=0; j < nClass; j++){
        Yval[j] = dN_CBWC_FC_vec[j][i];
      }
      gr[i] = new TGraph(nClass,Xval,Yval);
      gr[i]->SetName(Form("netCharge%d",i));
      gr[i]->SetTitle("#DeltaN");
      gr[i]->GetXaxis()->SetTitle("percentile%");
      gr[i]->GetYaxis()->SetTitle(Form("CBWC_FC_{%d}",i));
      gr[i]->SetMarkerStyle(33);
      gr[i]->SetMarkerSize(2);
      gr[i]->SetMarkerColor(kRed);
      gr[i]->SetLineColor(kGreen);

    }
    
    //for(int i=0; i<7; i++){
    //  TCanvas* c = new TCanvas(Form("c2_%d",i),Form("c2_%d",i),1600,1200);
    //  gr[i]->Draw("ALP");
    //  c->SaveAs(Form("cbwc_dN_%d.png",i));
    //}
*/

    
}
