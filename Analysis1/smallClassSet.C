void smallClassSet(){
  const int nSmallClass = 100;
  double smallClassLow[nSmallClass];
  double smallClassUp[nSmallClass];
  int classNumber = 0;
  for(int i = 0; i < nSmallClass; i++) { smallClassLow[i] = i * 100.0/double(nSmallClass) ; smallClassUp[i] = (i+1) * 100.0/double(nSmallClass) ;}
  smallClassLow[0] = -1.0; smallClassUp[nSmallClass-1] = 101.0;

  for(int i = 0; i < nSmallClass; i++) { 
    cout<<"i =  "<<i<<" :: "<<smallClassLow[i]<<" :: "<<smallClassUp[i]<<endl;
 }

 const int nClass = 5;
 std::vector<double> classLow(nClass);
 std::vector<double> classUp(nClass);
 // double classLow[nClass];
 // double classUp[nClass];
 for(int i = 0; i < nClass; i++) { classLow[i] = i * 100.0/double(nClass) ; classUp[i] = (i+1) * 100.0/double(nClass) ;}
 classLow[0] = -1.0; classUp[nClass-1] = 101.0;

 for(int iClass = 0; iClass < nClass; iClass++) { //upper bound was included
    cout<<"iClass =  "<<iClass<<" :: "<<classLow[iClass]<<" :: "<<classUp[iClass]<<endl;
  for(int iSmallClass = 0 ; iSmallClass < nSmallClass ; iSmallClass++){
    if( classLow[iClass] < smallClassLow[iSmallClass] && smallClassLow[iSmallClass] <= classUp[iClass] ){
      cout<<"iSmallClass = "<<iSmallClass<<"  :: its Class = "<<iClass<<" :: "<<smallClassLow[iSmallClass]<<" : "<<smallClassUp[iSmallClass]<<endl;
    }
    else {
        if(iSmallClass == 0 && iClass == 0 ) { 
          cout<<"iSmallClass = "<<iSmallClass<<"  :: its Class = "<<iClass<<" :: "<<smallClassLow[iSmallClass]<<" : "<<smallClassUp[iSmallClass]<<endl;
        }else {
          cout<<"iSmallClass = "<<iSmallClass <<" ERROR NO BIN "<<endl;
        }
    }
  }
 }




//  double sumForCBWC_Raa_nClass[nClass];
//  double sumForCBWC_Rbb_nClass[nClass];
//  double sumForCBWC_Rab_nClass[nClass];
//  double sumForCBWC_NDy_nClass[nClass];
 
//  for(int iClass = 0; iClass < nClass; iClass++) { //upper bound was included
//    double countSum = 0;
//    sumForCBWC_Raa_nClass[iClass] = 0;
//    sumForCBWC_Rbb_nClass[iClass] = 0;
//    sumForCBWC_Rab_nClass[iClass] = 0;
//    sumForCBWC_NDy_nClass[iClass] = 0;
//    // cout<<"iClass =  "<<iClass<<" :: "<<classLow[iClass]<<" :: "<<classUp[iClass]<<endl;
//    for(int iSmallClass = 0 ; iSmallClass < nSmallClass ; iSmallClass++){
//      if( classLow[iClass] < smallClassLow[iSmallClass] && smallClassLow[iSmallClass] <= classUp[iClass] ){
//        // cout<<"iSmallClass = "<<iSmallClass<<"  :: its Class = "<<iClass<<" :: "<<smallClassLow[iSmallClass]<<" : "<<smallClassUp[iSmallClass]<<endl;
//        countSum += nEntries_nSmallClass[iSmallClass];
//        sumForCBWC_Raa_nClass[iClass] +=  nEntries_nSmallClass[iSmallClass]*Raa_nSmallClass[iSmallClass];
//        sumForCBWC_Rbb_nClass[iClass] +=  nEntries_nSmallClass[iSmallClass]*Rbb_nSmallClass[iSmallClass];
//        sumForCBWC_Rab_nClass[iClass] +=  nEntries_nSmallClass[iSmallClass]*Rab_nSmallClass[iSmallClass];
//        sumForCBWC_NDy_nClass[iClass] +=  nEntries_nSmallClass[iSmallClass]*NDy_nSmallClass[iSmallClass];
//      } else {
//          if(iSmallClass == 0 && iClass == 0 ) {
//            // cout<<"iSmallClass = "<<iSmallClass<<"  :: its Class = "<<iClass<<" :: "<<smallClassLow[iSmallClass]<<" : "<<smallClassUp[iSmallClass]<<endl;
//            countSum += nEntries_nSmallClass[iSmallClass];
//            sumForCBWC_Raa_nClass[iClass] +=  nEntries_nSmallClass[iSmallClass]*Raa_nSmallClass[iSmallClass];
//            sumForCBWC_Rbb_nClass[iClass] +=  nEntries_nSmallClass[iSmallClass]*Rbb_nSmallClass[iSmallClass];
//            sumForCBWC_Rab_nClass[iClass] +=  nEntries_nSmallClass[iSmallClass]*Rab_nSmallClass[iSmallClass];
//            sumForCBWC_NDy_nClass[iClass] +=  nEntries_nSmallClass[iSmallClass]*NDy_nSmallClass[iSmallClass];
//          }else {
//            cout<<"iSmallClass = "<<iSmallClass <<" ERROR NO BIN "<<endl;
//          }
//      }
//    }//smallClass loop is over
//    Raa_CBWC[iClass] = sumForCBWC_Raa_nClass[iClass]/countSum;
//    Rbb_CBWC[iClass] = sumForCBWC_Rbb_nClass[iClass]/countSum;
//    Rab_CBWC[iClass] = sumForCBWC_Rab_nClass[iClass]/countSum;
//    NDy_CBWC[iClass] = sumForCBWC_NDy_nClass[iClass]/countSum;
//  }



}
