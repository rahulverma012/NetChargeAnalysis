template<typename T>
int findDigitCount(T number){
  // Its a digits counter function that counts no of didgits in that long int no. its 64 bit
  int digitCount = 0;
  if (number > 0) { digitCount +=1;}
  else            { digitCount +=2; number = -number;}
  // cout<<"number = "<<std::setw(15)<<static_cast<double>(number)<<endl;
  // cout<<"number = 1234567890123456789";
  //int64_t will work upto 19 digits. 19 digits are good
  while (static_cast<int64_t>(number)/10 > 0){
    //while loop will go unless the integer val or number is > 0
   number = static_cast<int64_t>(number)/10;
   digitCount++;
  }
  return digitCount;
}

template<typename T, typename U>
void PrintVariables(T PowerSum, U iEvtEntries, int nClass, int nOrder, std::string line1=" ", std::string line2 = " "){
  cout<<line1<<endl;
  int width_nClass = 0 ;
  for(int iClass = 0; iClass < nClass ; iClass++){
    if(width_nClass < findDigitCount(iEvtEntries[iClass]) ) { width_nClass = findDigitCount(iEvtEntries[iClass]);}
  }
  
  cout<<"width_nClass = "<<width_nClass<<endl;

  int width_nOrder[nOrder+1];
  for(int iOrder = 0; iOrder <= nOrder ; iOrder++){
    int digitCount = 0;
    for(int iClass = 0; iClass < nClass ; iClass++){
      if( digitCount < findDigitCount(PowerSum[iClass][iOrder]) ){ digitCount = findDigitCount(PowerSum[iClass][iOrder]); };
    }
    width_nOrder[iOrder] = digitCount;
    if(digitCount > 10) width_nOrder[iOrder] = 10+5;
  }

  for(int iClass = 0; iClass < nClass ; iClass++){
    cout<<line2<<" = "<<std::setw(findDigitCount(nClass))<<iClass<<" :: nEntries = "<<std::setw(width_nClass)<<iEvtEntries[iClass]<<" :: ";
    for(int iOrder = 0; iOrder <= nOrder ; iOrder++){
      cout<<"["<<iOrder<<"] = "<<std::setw(width_nOrder[iOrder])<<PowerSum[iClass][iOrder]<<" : ";
    }
    cout<<endl;
  }
  cout<<endl;
}

template<typename T, typename U>
void PrintVariables(T PowerSum, U iEvtEntries, int nSubSample, int nClass, int nOrder, std::string line1=" ", std::string line2 = " ", std::string line3= " "){
  cout<<line1<<endl;
  for(int iSubSample = 0; iSubSample < nSubSample; iSubSample++){
    PrintVariables(PowerSum[iSubSample], iEvtEntries[iSubSample], nClass, nOrder, line1+"["+std::to_string(iSubSample)+"]", line3);
  }
  cout<<endl;
}

template<typename T>
  int getFirstFilledBin(const T& histAxis){
    int firstFilledBin = -1;
    for(int iBin = 0 ; iBin <= histAxis->GetNbinsX()+1; iBin++){
      if(histAxis->GetBinContent(iBin) != 0) { firstFilledBin = iBin; break;}
    }
    if(firstFilledBin == -1){ cout<<"Histogram is empty"<<endl;}
    else if (firstFilledBin == 0) { cout<<"LastFilledBin is Underflow"<<endl;}
    return firstFilledBin;
  }

  template<typename T>
  int getLastFilledBin(const T& histAxis){
    int lastFilledBin = -1;
    for(int iBin = 0 ; iBin <= histAxis->GetNbinsX()+1; iBin++){
      if(histAxis->GetBinContent(iBin) != 0) { lastFilledBin = iBin;}
    }
    if(lastFilledBin == histAxis->GetNbinsX()+1){ cout<<"LastFilledBin is Overflow"<<endl;}
    else if (lastFilledBin == -1) { cout<<"Histogram is empty"<<endl;}
    return lastFilledBin;
  }

  enum MulitplicityAxisType{
    kIntTypeAxis=0,
    kFloatTypeAxis
  };


  template<typename T> 
  void getSmallBinProperEdges(const T& histAxis, const int& nSmallClass, vector<double>& smallClassLow, vector<double>& smallClassUp, 
                              vector<int>& iBinLow_smallClass, vector<int>& iBinUp_smallClass, vector<int>& nClassIDX_smallClass, 
                              const int &nClass, const vector<double>& classLow, const vector<double>& classUp,
                              const int& firstFilledBin, const int& lastFilledBin, const int mode, std::string var = "mult"){
    cout<<"Processing smallClass information for Multiplicity/Centraility bin width correction"<<endl;
    int binPositionsLow  = firstFilledBin;
    // int binPositionsHigh = lastFilledBin;
    for(int i = 0; i < nSmallClass; i++) {
      smallClassLow[i] = histAxis->GetBinLowEdge(binPositionsLow + i) ; 
      //For integer Bins 
      if      (mode == kIntTypeAxis){ smallClassUp[i]  = smallClassLow[i] ;}//+ histAxis->GetBinWidth(i);
      else if (mode == kFloatTypeAxis){ smallClassUp[i]  = histAxis->GetBinLowEdge(binPositionsLow + i) + histAxis->GetBinWidth(binPositionsLow + i);}
      for(int iClass = 0 ; iClass < nClass+2; iClass++) {
        if( smallClassUp[i] < classLow[0]) { nClassIDX_smallClass[i] = 0; break; } 
        if      (mode == kIntTypeAxis  ){if( classLow[iClass] <= smallClassLow[i]  && smallClassUp[i]  <= classUp[iClass]) { nClassIDX_smallClass[i] = iClass+1; break; }} 
        else if (mode == kFloatTypeAxis){if( classLow[iClass] <= smallClassLow[i]  && smallClassUp[i]  <= classUp[iClass]) { nClassIDX_smallClass[i] = iClass+1; break; }}
        if( classUp[nClass-1] <= smallClassLow[i]) { nClassIDX_smallClass[i] = nClass+1; break; } 
      }
    }
    // smallClassLow[0] = -1.0; smallClassUp[nSmallClass-1] = 101.0;

    int64_t totEntryCounter = 0;
    for(int i=0 ; i < nSmallClass ; i++){
      iBinLow_smallClass[i] = histAxis->FindBin(smallClassLow[i]);
      if     (mode == kIntTypeAxis  ){ iBinUp_smallClass[i]  = histAxis->FindBin(smallClassUp[i]);}  //Here smallClassLow[i] and smallClassUp[i] are same
      else if(mode == kFloatTypeAxis){ iBinUp_smallClass[i]  = histAxis->FindBin(smallClassLow[i]);} //Here smallClassLow[i] and smallClassUp[i] are differernt
      cout<<"i = "<<std::setw(3)<<i<<" :: "<<var<<" \u2208 "<<"["<<std::setw(5)<<smallClassLow[i]<<","<<std::setw(5)<<smallClassUp[i];
      if     (mode == kIntTypeAxis   ){ cout<<"] ";}
      else if(mode == kFloatTypeAxis ){ cout<<") ";}
      cout<<" :: binPos \u2208 ["<<std::setw(5)<<iBinLow_smallClass[i]<<","<<std::setw(5)<<iBinUp_smallClass[i]<<"] "
      <<" :: nClassIDX = "<<nClassIDX_smallClass[i]
      <<" :: \u222B = "<<std::setw(11)<<histAxis->Integral(iBinLow_smallClass[i], iBinUp_smallClass[i])
      <<" :: \u0025 = "<<100.0*static_cast<double>(histAxis->Integral(iBinLow_smallClass[i], iBinUp_smallClass[i]))/static_cast<double>(histAxis->GetEntries())
      <<endl;
      if     (mode == kIntTypeAxis  ){totEntryCounter += histAxis->Integral(histAxis->FindBin(smallClassLow[i]), histAxis->FindBin(smallClassUp[i]));}  //Here smallClassLow[i] and smallClassUp[i] are same
      else if(mode == kFloatTypeAxis){totEntryCounter += histAxis->Integral(histAxis->FindBin(smallClassLow[i]), histAxis->FindBin(smallClassLow[i]));} //Here smallClassLow[i] and smallClassUp[i] are differernt
    }
    if(totEntryCounter != histAxis->GetEntries()){ cout<<"ERROR :: total entries in classes are not matching"<<endl;}
  }

  template<typename T>
  void checkClassificationHist(const std::vector<T>& h1D_nClass, const int& axisCl, const int& nClass){
    for(int iClass = 0 ; iClass < nClass-1; iClass++){
      int binDiff = getFirstFilledBin(h1D_nClass[axisCl][iClass+1]) - getLastFilledBin(h1D_nClass[axisCl][iClass]);
      if(binDiff != 1){ 
        cout<<"ERROR :: Severe :: bin classification Improper :: "
        <<" :: lastBin  = "<<getLastFilledBin(h1D_nClass[axisCl][iClass])
        <<" :: firstBin = "<<getFirstFilledBin(h1D_nClass[axisCl][iClass+1])
        <<endl;
      }
    }
  }

  //  template<typename T> 
  // void getSmallBinProperEdges(const T& histAxis, const int& nSmallClass, vector<double>& smallClassLow, vector<double>& smallClassUp, 
  //                             vector<int>& iBinLow_smallClass, vector<int>& iBinUp_smallClass, vector<int>& nClassIDX_smallClass, 
  //                             const int &nClass, const vector<double>& FT0M_classLow,const vector<double>& FT0M_classUp,
  //                             const int& firstFilledBin, const int& lastFilledBin){
  //   cout<<"Processing smallClass information for Multiplicity/Centraility bin width correction"<<endl;
  //   int binPositionLow = firstFilledBin;
  //   for(int i=0; i < nSmallClass; i++){     
  //     smallClassLow[i] = histAxis->GetBinLowEdge(binPositionLow + i);
  //     smallClassUp[i]  = histAxis->GetBinLowEdge(binPositionLow + i) + histAxis->GetBinWidth(binPositionLow + i);

  //     nClassIDX_smallClass[i] = -1;// default initialization

  //     for(int iClass=0; iClass < nClass; iClass++){
  //       if(FT0M_classLow[iClass]<= smallClassLow[i] && smallClassUp[i] <= FT0M_classUp[iClass]){
  //         nClassIDX_smallClass[i] = iClass;
  //         break; 
  //       }
  //     }
  //   }

  //   int64_t totEntryCounter = 0;
  //   double totalPercentge =0 ;
  //   for(int i=0 ; i < nSmallClass ; i++){
  //     iBinLow_smallClass[i] = histAxis->FindBin(smallClassLow[i]);
  //     iBinUp_smallClass[i]  = histAxis->FindBin(smallClassUp[i]) -1;
  //     cout<<"i = "<<std::setw(3)<<i<<" :: cent \u2208 "<<"["<<std::setw(5)<<smallClassLow[i]<<","<<std::setw(5)<<smallClassUp[i]<<"] "
  //     <<" :: binPos \u2208 ["<<std::setw(5)<<iBinLow_smallClass[i]<<","<<std::setw(5)<<iBinUp_smallClass[i]<<"] "
  //     <<" :: nClassIDX = "<<nClassIDX_smallClass[i]
  //     <<" :: \u222B = "<<std::setw(11)<<histAxis->Integral(iBinLow_smallClass[i], iBinUp_smallClass[i])
  //     <<" :: \u0025 = "<<100.0*static_cast<double>(histAxis->Integral(iBinLow_smallClass[i], iBinUp_smallClass[i]))/static_cast<double>(histAxis->GetEntries())
  //     <<endl;
  //     totEntryCounter += histAxis->Integral(iBinLow_smallClass[i], iBinUp_smallClass[i]);
  //   }
  //   if(totEntryCounter != histAxis->GetEntries()){ cout<<"ERROR :: total entries in classes are not matching"<<endl;}
  // }

  template<typename T> 
  int64_t FindEntriesPerSubSample(const T& hist, const int &nSubSample){
    cout<<endl<<"Finding Entries Per SubSample"<<endl;
    cout<<hist->GetName()<<" Total Entries = "<<std::setprecision(20)<<hist->GetEntries()<<endl;
    int64_t totalEntries = hist->GetEntries();
    int64_t Entries = totalEntries/nSubSample;
    int64_t Remainder  = totalEntries%nSubSample;
    cout<<"Entries              = "<<Entries<<endl; //9995026
    cout<<"Remainder            = "<<Remainder<<endl;
    if( Remainder != 0){Entries++;}
    cout<<endl;
    // cout<<"Entries Per SubSample = "<<Entries<<endl;
    return Entries;
  }

  template<typename T> 
  int64_t FindEntriesPerSubSample(const std::vector<T> &histList, int nSubSample, int64_t &totalEntries){
    cout<<endl<<"Finding Entries Per SubSample"<<endl;
    totalEntries = 0;
    int nHist = histList.size();
    for(int iHist = 0 ; iHist < nHist; iHist++){
      totalEntries += histList[iHist]->GetEntries();
    }
    cout<<histList[0]->GetName()<<" :: total Entries = "<<std::setprecision(20)<<totalEntries<<endl;
    int64_t Entries = totalEntries/nSubSample;
    int64_t Remainder  = totalEntries%nSubSample;
    cout<<"Entries              = "<<Entries<<endl; //9995026
    cout<<"Remainder            = "<<Remainder<<endl;
    if( Remainder != 0){Entries++;}
    cout<<endl;
    cout<<"Entries Per SubSample = "<<Entries<<endl;
    return Entries;
  }
void Get_CBWC_C(vector<vector<double>> &dN_CBWC_C_vec,  vector<vector<double>> &dN_C_NM_nSmallClass_vec, const int& nClass, const int& nOrder, const int& nSmallClass, const vector<double> classLow,const vector<double> classUp,const vector<double> smallClassLow,const vector<double> smallClassUp, const vector<int64_t> iEvtEntries_nSmallClass_vec ){
  vector<vector<double>>SumForCBWC_C(nClass,vector<double>(nOrder+1,0));
  if(dN_C_NM_nSmallClass_vec.size()!= nSmallClass){cout<<"Error in cbwc"<<endl;}
  for(int iClass=0; iClass<nClass; iClass++){
    double countSum =0;
    for(int iOrder=0; iOrder< nOrder+1; iOrder++){
    SumForCBWC_C[iClass][iOrder]=0;}
   for(int iSmallClass=0; iSmallClass< nSmallClass; iSmallClass++){
    if(classLow[iClass]<= smallClassLow[iSmallClass] && smallClassUp[iSmallClass] <= classUp[iClass]){
      countSum += iEvtEntries_nSmallClass_vec[iSmallClass];
      for(int iOrder=0; iOrder<=nOrder; iOrder++){
      SumForCBWC_C[iClass][iOrder] += dN_C_NM_nSmallClass_vec[iSmallClass][iOrder]*double(iEvtEntries_nSmallClass_vec[iSmallClass]);
      }
     }
    }
    for(int iOrder=0; iOrder<=nOrder; iOrder++){
      dN_CBWC_C_vec[iClass][iOrder] = SumForCBWC_C[iClass][iOrder]/countSum;
    }
   }
  }

  void Get_NM_FromPowerSum_Err(vector<vector<vector<double>>> &dN_PowerSum_nSubSample_nSmallClass_CBWC_vec,const vector<vector<int64_t>> iEvtEntries_nSubSample_nSmallClass_CBWC_vec, vector<vector<vector<double>>> &dN_NM_nSubSample_nSmallClass_CBWC, const int& nSubSample, const int& nSmallClass,const int& nOrder){
  cout<<endl<<"obtaining Normal Moments for Error Calcucation"<<endl;
  for(int iSubSample = 0; iSubSample < nSubSample ; iSubSample++){
    for(int iSmallClass = 0; iSmallClass < nSmallClass ; iSmallClass++){
      for(int iOrder = 0; iOrder <= nOrder ; iOrder++){
        if(iEvtEntries_nSubSample_nSmallClass_CBWC_vec[iSubSample][iSmallClass]!= 0){ //Event Entries of that SubSample must not be zero.
          dN_NM_nSubSample_nSmallClass_CBWC[iSubSample][iSmallClass][iOrder] = double(dN_PowerSum_nSubSample_nSmallClass_CBWC_vec[iSubSample][iSmallClass][iOrder])/double(iEvtEntries_nSubSample_nSmallClass_CBWC_vec[iSubSample][iSmallClass]);
        }else{
          dN_NM_nSubSample_nSmallClass_CBWC[iSubSample][iSmallClass][iOrder] = 0; 
        }
      }
      // cout<<endl;
    }
  }

}

void Get_C_From_NM_Err(vector<vector<vector<double>>> &dN_NM_nSubSample_nSmallClass_CBWC, vector<vector<vector<double>>> &dN_C_NM_nSubSample_nSmallClass_CBWC, const int& nSubSample, const int& nSmallClass){
  cout<<endl<<"obtaining Cumulants from normal moments for Error Calcucation"<<endl;
  double X0, X1, X2, X3, X4, X5, X6;
  for(int iSubSample = 0; iSubSample < nSubSample ; iSubSample++){
    for(int iSmallClass = 0; iSmallClass < nSmallClass ; iSmallClass++){
      X0 = dN_NM_nSubSample_nSmallClass_CBWC[iSubSample][iSmallClass][0];
      X1 = dN_NM_nSubSample_nSmallClass_CBWC[iSubSample][iSmallClass][1];
      X2 = dN_NM_nSubSample_nSmallClass_CBWC[iSubSample][iSmallClass][2];
      X3 = dN_NM_nSubSample_nSmallClass_CBWC[iSubSample][iSmallClass][3];
      X4 = dN_NM_nSubSample_nSmallClass_CBWC[iSubSample][iSmallClass][4];
      X5 = dN_NM_nSubSample_nSmallClass_CBWC[iSubSample][iSmallClass][5];
      X6 = dN_NM_nSubSample_nSmallClass_CBWC[iSubSample][iSmallClass][6];

      dN_C_NM_nSubSample_nSmallClass_CBWC[iSubSample][iSmallClass][0] = X0;
      dN_C_NM_nSubSample_nSmallClass_CBWC[iSubSample][iSmallClass][1] = X1;
      dN_C_NM_nSubSample_nSmallClass_CBWC[iSubSample][iSmallClass][2] = X2 - TMath::Power(X1,2);
      dN_C_NM_nSubSample_nSmallClass_CBWC[iSubSample][iSmallClass][3] = X3 - 3.0*X2*X1 + 2.0*TMath::Power(X1,3);
      dN_C_NM_nSubSample_nSmallClass_CBWC[iSubSample][iSmallClass][4] = X4 - 4.0*X3*X1 - 3.0*TMath::Power(X2,2) + 12.0*X2*TMath::Power(X1,2) - 6.0*TMath::Power(X1,4);
      dN_C_NM_nSubSample_nSmallClass_CBWC[iSubSample][iSmallClass][5] = X5 - 5.0*X4*X1 -10.0*X3*X2 + 20*X3*TMath::Power(X1,2) + 30.0*TMath::Power(X2,2)*X1 - 60.0*X2*TMath::Power(X1,3) + 24.0*TMath::Power(X1,5);
      dN_C_NM_nSubSample_nSmallClass_CBWC[iSubSample][iSmallClass][6] = X6 - 6.0*X5*X1 -15.0*X4*X2 + 30*X4*TMath::Power(X1,2) - 10.0*TMath::Power(X3,2)    +120.0*X3*X2*X1 -120.0*X3*TMath::Power(X1,3) + 30.0*TMath::Power(X2,3) -270.0*TMath::Power(X2,2)*TMath::Power(X1,2)+360*X2*TMath::Power(X1,4)-120.0*TMath::Power(X1,6);
    }
  }

}

void Get_FC_From_C(vector<vector<vector<double>>> &dN_C_NM_nSubSample_nSmallClass_CBWC, vector<vector<vector<double>>> &dN_FC_C_nSubSample_nSmallClass_CBWC, const int& nSubSample, const int& nSmallClass){
  cout<<"obtaining factorial cumulants from cumulants for error calculation"<<endl;
  double C0, C1, C2, C3, C4, C5, C6;
  for(int iSubSample = 0; iSubSample<nSubSample; iSubSample++){
    for(int iSmallClass =0; iSmallClass<nSmallClass; iSmallClass++){
      C0 = dN_C_NM_nSubSample_nSmallClass_CBWC[iSubSample][iSmallClass][0];
      C1 = dN_C_NM_nSubSample_nSmallClass_CBWC[iSubSample][iSmallClass][1];
      C2 = dN_C_NM_nSubSample_nSmallClass_CBWC[iSubSample][iSmallClass][2];
      C3 = dN_C_NM_nSubSample_nSmallClass_CBWC[iSubSample][iSmallClass][3];
      C4 = dN_C_NM_nSubSample_nSmallClass_CBWC[iSubSample][iSmallClass][4];
      C5 = dN_C_NM_nSubSample_nSmallClass_CBWC[iSubSample][iSmallClass][5];
      C6 = dN_C_NM_nSubSample_nSmallClass_CBWC[iSubSample][iSmallClass][6];

      dN_FC_C_nSubSample_nSmallClass_CBWC[iSubSample][iSmallClass][0] = C0;
      dN_FC_C_nSubSample_nSmallClass_CBWC[iSubSample][iSmallClass][1] = C1;
      dN_FC_C_nSubSample_nSmallClass_CBWC[iSubSample][iSmallClass][2] = C2 - C1;
      dN_FC_C_nSubSample_nSmallClass_CBWC[iSubSample][iSmallClass][3] = C3 -  3*C2 +  2*C1;
      dN_FC_C_nSubSample_nSmallClass_CBWC[iSubSample][iSmallClass][4] = C4 -  6*C3 + 11*C2 -   6*C1;
      dN_FC_C_nSubSample_nSmallClass_CBWC[iSubSample][iSmallClass][5] = C5 - 10*C4 + 35*C3 -  50*C2 +  24*C1;
      dN_FC_C_nSubSample_nSmallClass_CBWC[iSubSample][iSmallClass][6] = C6 - 15*C5 + 85*C4 - 225*C3 + 274*C2 - 120*C1;
    }
  }
}

void Get_delta_of_C(vector<vector<vector<double>>> &dN_C_NM_nSubSample_nSmallClass_CBWC, vector<vector<double>> &dN_delta_C_NM_nSubSample_nSmallClass, const int& nSubSample, const int& nSmallClass, const int& nOrder)
{
    cout << "getting delta for error" << endl;

    for (int iSmallClass = 0; iSmallClass < nSmallClass; iSmallClass++) {
        for (int iOrder = 0; iOrder <= nOrder; iOrder++) {
            double cumulantSum = 0.0;

            // First: compute average cumulant over subsamples
            for (int iSubSample = 0; iSubSample < nSubSample; iSubSample++) {
                cumulantSum += dN_C_NM_nSubSample_nSmallClass_CBWC[iSubSample][iSmallClass][iOrder];
            }

            double averageCumulant = cumulantSum / static_cast<double>(nSubSample);

            // Then: compute standard deviation (sample variance)
            double cumulantErrorSum = 0.0;
            for (int iSubSample = 0; iSubSample < nSubSample; iSubSample++) {
                double diff = dN_C_NM_nSubSample_nSmallClass_CBWC[iSubSample][iSmallClass][iOrder] - averageCumulant;
                cumulantErrorSum += diff * diff;
            }

            double sigma = TMath::Sqrt(cumulantErrorSum / static_cast<double>(nSubSample - 1));
            dN_delta_C_NM_nSubSample_nSmallClass[iSmallClass][iOrder] = sigma / TMath::Sqrt(nSubSample);
        }
    }
}

void Get_Error(vector<vector<double>> &dN_delta_C_NM_nSubSample_nSmallClass, vector<vector<double>> &dN_Error_C_NM, const int& nClass, const int& nOrder, const int& nSmallClass, const vector<double> classLow,const vector<double> classUp,const vector<double> smallClassLow,const vector<double> smallClassUp, const vector<int64_t> iEvtEntries_nSmallClass_vec){
  cout<<"getting error ******"<<endl;
  for(int iClass=0; iClass<nClass; iClass++){
    double totalEvents = 0.0;
    vector<double> errorSquaredSum(nOrder+1,0);
    for(int iSmallClass=0; iSmallClass<nSmallClass; iSmallClass++ ){
      if(classLow[iClass]<= smallClassLow[iSmallClass] && smallClassUp[iSmallClass] <= classUp[iClass]){
        totalEvents += iEvtEntries_nSmallClass_vec[iSmallClass];
      }
    }//smallcls1
     if (totalEvents == 0.0) {
      cout << "Warning: totalEvents = 0 in wide class " << iClass << endl;
      continue;
    }

   for(int iSmallClass=0; iSmallClass<nSmallClass; iSmallClass++){
    if(classLow[iClass]<= smallClassLow[iSmallClass] && smallClassUp[iSmallClass] <= classUp[iClass]){
      double weight = static_cast<double>(iEvtEntries_nSmallClass_vec[iSmallClass])/totalEvents;
      for(int iOrder=0; iOrder<=nOrder; iOrder++){
        double delta = dN_delta_C_NM_nSubSample_nSmallClass[iSmallClass][iOrder];
        errorSquaredSum[iOrder] += weight*weight*delta*delta;
      }
    }
   }//smllcls2
   for(int iOrder=0; iOrder<=nOrder; iOrder++){
    dN_Error_C_NM[iClass][iOrder] = TMath::Sqrt(errorSquaredSum[iOrder]);

   }
    
  }//nclass
}//void

void Plot_CBWC_With_Error(const vector<vector<double>>& dN_CBWC_C_vec,           // [nClass][nOrder+1]
    const vector<vector<double>>& dN_Error_C_NM,           // [nClass][nOrder+1]
    const vector<double>& classLow,                        // class lower edge
    const vector<double>& classUp,                         // class upper edge
    int nClass,
    int orderToPlot,                                       // which order to plot
    const char* yLabel = "Cumulant",
    const char* title = "CBWC Cumulant with Error"
){
  std::vector<double> x(nClass), ex(nClass), y(nClass), ey(nClass);
  for(int iClass=0; iClass<nClass; iClass++){
    x[iClass] = 0.5*(classLow[iClass]+classUp[iClass]);
    ex[iClass]= 0;
    y[iClass]= dN_CBWC_C_vec[iClass][orderToPlot];
    ey[iClass]= dN_Error_C_NM[iClass][orderToPlot];
  }
  TGraphErrors* g = new TGraphErrors(nClass, &x[0], &y[0], &ex[0], &ey[0]);
    g->SetTitle(title);
    g->SetMarkerStyle(20);
    g->SetMarkerColor(kBlue+2);
    g->SetLineColor(kBlue+2);
    g->SetLineWidth(1);
    g->GetXaxis()->SetTitle("NchV0M(%)");
    g->GetYaxis()->SetTitle(yLabel);

    g->Draw("AP");
}
template <typename T>
void Plot_Cumulant_By_ClassIndex(
    T &gr,
    const vector<vector<double>>& dN_CBWC_C_vec,     // [nClass][nOrder+1]
    const vector<vector<double>>& dN_Error_C_NM,     // [nClass][nOrder+1]
    int nClass,
    int orderToPlot,                                 // which cumulant order
    const char* yLabel = "Cumulant",
    const char* title = "Cumulant vs Class Index"
) {
    std::vector<double> x(nClass), ex(nClass, 0.0), y(nClass), ey(nClass);

    for (int i = 0; i < nClass; ++i) {
        x[i] = (i + 0.5)*20;                               // class center (e.g., 0.5, 1.5, ...)
        y[i] = dN_CBWC_C_vec[i][orderToPlot];         // cumulant value
        ey[i] = dN_Error_C_NM[i][orderToPlot];        // statistical error
    }

    TGraphErrors* g = new TGraphErrors(nClass, &x[0], &y[0], &ex[0], &ey[0]);
    g->SetTitle(title);
    g->SetMarkerStyle(20);
    g->SetMarkerColor(kBlue+1);
    g->SetLineColor(kBlue+1);
    g->SetLineWidth(2);
    g->GetXaxis()->SetTitle("Multiplicity Percentile(%)");
    g->GetYaxis()->SetTitle(yLabel);

    g->Draw("AP");
    
    g->SetName(Form("netchargeFC_%d",orderToPlot));
    gr[orderToPlot] = g;

    // Optional: Improve axis style
    //gPad->SetGrid();
}
