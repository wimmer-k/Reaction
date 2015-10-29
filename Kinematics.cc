#include "Kinematics.hh"
double deg2rad=PI/180.;
double rad2deg=180./PI;

using namespace std;

Kinematics::Kinematics(Nucleus* projectile, Nucleus* target, double ebeam){
  fParticle[0] = projectile;
  fParticle[1] = target;
  fParticle[2] = NULL;
  fParticle[3] = NULL;
  fM[0]=fParticle[0]->GetMass();
  fM[1]=fParticle[1]->GetMass();
  fEBeam = ebeam;
  fQValue =0;
  Initial();
  FinalCm();
}
Kinematics::Kinematics(Nucleus* projectile, Nucleus* target, Nucleus* recoil, Nucleus* ejectile, double ebeam, double ex3){
  fParticle[0] = projectile;
  fParticle[1] = target;
  fParticle[2] = recoil;
  fParticle[3] = ejectile;
  for(int i=0;i<4;i++)
    fM[i]=fParticle[i]->GetMass();
    
  fEBeam = ebeam;
  fQValue = (fM[0]+fM[1])-(fM[2]+fM[3])-ex3;
  Initial();
  FinalCm();
}

TSpline3* Kinematics::Evslab(double thmin, double thmax, double size, int part){
  //cout << "maximum scattering angle: " << GetMaxAngle(fVcm[part])*180./PI  << endl;
  //cout << "max " << thmax << " min " << thmin << " steps " << (int)((thmax-thmin)/size)+1 << endl;
  double* energy = new double[(int)((thmax-thmin)/size)+1];
  double* angle = new double[(int)((thmax-thmin)/size)+1];
  int number =0;
  for(int i=0;i<((thmax-thmin)/size);i++){
    Final((thmin+i*size)*deg2rad,2);
    angle[i]=GetThetalab(part)*rad2deg;
    energy[i]=GetTlab(part)*1000;
    if(energy[i]<1e15||energy[i]>0.0)
      number++;
    else
      break;
    //cout << setprecision(4) << GetThetacm(3)/deg2rad << "\t" << setprecision(4) << GetThetacm(2)/deg2rad << "\t" << setprecision(4) << GetThetalab(3)/deg2rad << "\t" << setprecision(4) << GetThetalab(2)/deg2rad << "\t" << setprecision(4) << GetTlab(3) << "\t" << setprecision(4) << GetTlab(2) << endl;  
    //cout << "theta " <<(thmin+i*size)<<" res "<< angle[i] << " energy " << energy[i] << endl;
  }
  TGraph* graph = new TGraph(number, angle, energy);
  TSpline3* spline = new TSpline3("ETh_lab",graph);
  delete graph;
  delete[] angle;
  delete[] energy;
  return spline;
}
TSpline3* Kinematics::Evscm(double thmin, double thmax, double size, int part){
  double* energy = new double[(int)((thmax-thmin)/size)+1];
  double* angle = new double[(int)((thmax-thmin)/size)+1];
  int number =0;
  for(int i=0;i<((thmax-thmin)/size);i++){
    Final((thmin+i*size)*deg2rad,2);
    angle[i]=GetThetacm(part)*rad2deg;
    energy[i]=GetTlab(part)*1000;
    number++;
  }
  TGraph* graph = new TGraph(number, angle, energy);
  TSpline3* spline = new TSpline3("ETh_cm",graph);
  delete graph;
  delete[] angle;
  delete[] energy;
  return spline;
}

double Kinematics::GetExcEnergy(TLorentzVector recoil){
  TLorentzVector ejectile;
  recoil.Boost(0,0,-GetBetacm()); //boost to cm system
  
  ejectile.SetVect( -recoil.Vect() ); //pr = -pe
  ejectile.SetE(GetCmEnergy()*1000. - recoil.E()); //Ee=Ecm-Er
  
  ejectile.Boost(0,0,GetBetacm()); //boost to lab system

  double eex = ejectile.M() - fParticle[3]->GetMass()*1000.;

  return eex;							    
}
double Kinematics::GetBeamEnergy(double LabAngle, double LabEnergy){
  double ProjectileMass=fM[0]*1000;
  double TargetMass=fM[1]*1000;

  double ts = pow(TargetMass,2);
  double ps = pow(ProjectileMass,2);
  double cs = pow(cos(LabAngle),2);
  double es = pow(LabEnergy,2);
  double te = TargetMass*LabEnergy;

  return (-8*ProjectileMass - 4*TargetMass + LabEnergy/cs - 2*LabEnergy*tan(LabAngle) 
	  + sqrt(16*ts*cs + LabEnergy*(LabEnergy/cs*pow(cos(LabAngle) - sin(LabAngle),4) + 8*TargetMass*(3 + sin(2*LabAngle))))/cos(LabAngle) 
	  + sqrt((24*ts*LabEnergy*cs + 
		  32*ps*TargetMass*pow(cos(LabAngle),4) + 
		  2*TargetMass*es*pow(cos(LabAngle) - sin(LabAngle),4) + 
		  16*ts*LabEnergy*pow(cos(LabAngle),3)*sin(LabAngle) - 
		  8*ps*LabEnergy*cs*(sin(2*LabAngle) - 1) + 
		  2*ps*cos(3*LabAngle)*sqrt(2*(4*ts + 12*te + es) + 
					    (8*ts - 2*es)*cos(2*LabAngle) + 
					    LabEnergy*(LabEnergy/cs + 8*TargetMass*sin(2*LabAngle) - 4*LabEnergy*tan(LabAngle))) + 
		  2*cos(LabAngle)*(3*ps + te - te*sin(2*LabAngle))*sqrt(2*(4*ts + 12*te + es) + (8*ts - 2*es)*cos(2*LabAngle) + LabEnergy*(LabEnergy/cs + 8*TargetMass*sin(2*LabAngle) - 4*LabEnergy*tan(LabAngle))))/(pow(cos(LabAngle),4)*TargetMass)))/8.;
}


void Kinematics::Initial(){
  fT[0]=fEBeam;
  fT[1]=0;
  fE[0]=E_tm(fT[0],fM[0]);
  fE[1]=E_tm(fT[1],fM[1]);
  fP[0]=sqrt(fT[0]*fT[0]+2*fT[0]*fM[0]);
  fP[1]=0;
  fV[0]=V_pe(fP[0],fE[0]);
  fV[1]=V_pe(fP[1],fE[1]);

  fEcm[0]=GetCmEnergy(fEBeam)/2+(fM[0]*fM[0]-fM[1]*fM[1])/(2*GetCmEnergy(fEBeam));
  fEcm[1]=GetCmEnergy(fEBeam)/2-(fM[0]*fM[0]-fM[1]*fM[1])/(2*GetCmEnergy(fEBeam));
  fTcm[0]=fEcm[0]-fM[0];
  fTcm[1]=fEcm[1]-fM[1];
  fPcm[0]=Pcm_em(fEcm[0],fM[0]);
  fPcm[1]=Pcm_em(fEcm[1],fM[1]);
  fVcm[0]=V_pe(fPcm[0],fEcm[0]);
  fVcm[1]=V_pe(fPcm[1],fEcm[1]);

  fBeta_cm = (fP[0]-fP[1])/(fE[0]+fE[1]);
  fBetacm[0] = betacm_tm(fTcm[0],fM[0]);
  fBetacm[1] = -betacm_tm(fTcm[1],fM[1]);
  fGamma_cm = 1/sqrt(1-fBeta_cm*fBeta_cm);
  fTCm_i = GetCmEnergy(fEBeam)-fM[0]-fM[1];
  fTCm_f = fTCm_i + fQValue;

}
void Kinematics::FinalCm(){
  if(fParticle[2]==NULL && fParticle[3]==NULL){
    cout << "warning: outgoing particles not defined! Assuming elastic scattering" << endl;
    cout << "recoil = target" << endl;
    cout << "ejectile = projectile" << endl;
    fM[2]=fParticle[1]->GetMass();
    fM[3]=fParticle[0]->GetMass();     
  }
  fTcm[2]=fTCm_f/2*(fTCm_f+2*fM[3])/GetCmEnergy(fEBeam);
  fTcm[3]=fTCm_f/2*(fTCm_f+2*fM[2])/GetCmEnergy(fEBeam);
  fEcm[2]=E_tm(fTcm[2],fM[2]);
  fEcm[3]=E_tm(fTcm[3],fM[3]);
  fPcm[2]=Pcm_em(fEcm[2],fM[2]);
  fPcm[3]=Pcm_em(fEcm[3],fM[3]);
  fVcm[2]=V_pe(fPcm[2],fEcm[2]);
  fVcm[3]=V_pe(fPcm[3],fEcm[3]);
  fBetacm[2]=-betacm_tm(fTcm[2],fM[2]);
  fBetacm[3]=betacm_tm(fTcm[3],fM[3]);

  /*
  for(int i=0;i<4;i++){
    cout << "fBetacm["<<i<<"] = " << fBetacm[i]<< "\t";
  }
  cout <<endl;
  for(int i=0;i<4;i++){
    cout << "fPcm["<<i<<"] = " << fPcm[i]<< "\t";
  }
  cout <<endl;
  cout << "fBeta_cm = " << fBeta_cm<<endl;
  */
}
void Kinematics::Final(double angle, int part){//angle of proton in lab system
  if(angle>GetMaxAngle(fVcm[part]))
    SetAngles(0, part);
  else
    SetAngles(angle, part);
  fE[2]=E_final(2);
  fE[3]=E_final(3);
  fT[2]=T_final(2);
  fT[3]=T_final(3);
  //fP[2]=fGamma_cm*(fPcm[2]+fBeta_cm*fEcm[2]);
  //cout << "Vr = " << V_pe(fP[2],fE[2]) << endl; 
  //fP[3]=fGamma_cm*(fPcm[3]+fBeta_cm*fEcm[3]);
  //cout << "Ve = " << V_pe(fP[3],fE[3]) << endl; 
  fP[2]=P_tm(fT[2],fM[2]);
  //fP[2]=fGamma_cm*fPcm[2]*(cos(fThetacm[2])+fBeta_cm/fVcm[2]);
  //cout << fTheta[2]*180./PI << "\t" << GetTlab(2)  << "\tVr = " << V_pe(fP[2],fE[2]) << "\tPr = " << fP[2] << "\t"; 
  fP[3]=P_tm(fT[3],fM[3]);
  //fP[3]=fGamma_cm*fPcm[3]*(cos(fThetacm[3])+fBeta_cm/fVcm[3]);
  //cout << fTheta[3]*180./PI << "\t" << GetTlab(3)  << "\tVe = " << V_pe(fP[3],fE[3]) << "\tPe = " << fP[3] << endl; 
  fV[2]=V_pe(fP[2],fE[2]);
  fV[3]=V_pe(fP[3],fE[3]);

  fQ = fPcm[0]*fPcm[0] + fPcm[2]*fPcm[2] - 2*fPcm[0]*fPcm[2]*cos(GetThetacm(2));
}
double Kinematics::ELab(double angle_lab, int part){
  Final(angle_lab, part);
  return GetTlab(part);
}
void Kinematics::SetAngles(double angle, int part, bool upper){
  int given;
  int other;
  if(part==2){
    given =2;
    other =3;
  }
  else if(part==3){
    given =3;
    other =2;
  }
  else{
    cout << " error in Kinematics::SetAngles("<<angle<<", "<<part<<") "<<endl;
    cout << " part must be 2 or 3 " << endl;
    exit(4);
  } 
  fTheta[given]=angle;
  fThetacm[given]=Angle_lab2cm(fVcm[given],fTheta[given]);
  if(given==3&&(fParticle[0]->GetMass()>fParticle[1]->GetMass())){
    //cout << "inverse kinematics" << endl;
    fThetacm[given]=Angle_lab2cminverse(fVcm[given],fTheta[given],upper);
  }
  fThetacm[other]=PI-fThetacm[given];
  if(fTheta[given]==0)
    fTheta[other]=PI/given;
  else
    fTheta[other]=Angle_cm2lab(fVcm[other],fThetacm[other]);

//  fTheta[3]=angle;
//  fThetacm[3]=Angle_lab2cm(fVcm[3],fTheta[3]);
//  fThetacm[2]=PI-fThetacm[3];
//  if(fTheta[3]==0)
//    fTheta[2]=PI/2;
//  else
//    fTheta[2]=Angle_cm2lab(fVcm[2],fThetacm[2]);
}


double Kinematics::GetCmEnergy(double ebeam){
  double ecm;
  ecm = sqrt(fM[0]*fM[0]+fM[1]*fM[1]+2.*fM[1]*(fM[0]+ebeam));
  return ecm;
}
double Kinematics::GetCmEnergy(){
  return GetCmEnergy(fEBeam);
}
double Kinematics::NormalkinEnergy(){
  double ENorm = (GetCmEnergy(fEBeam)*GetCmEnergy(fEBeam)-fM[0]*fM[0]-fM[1]*fM[1])/(2*fM[0])-fM[1];
  return ENorm;
}

double Kinematics::GetMaxAngle(double vcm){
  double x;
  x = fBeta_cm/vcm;
  if(x*x<1)
    return PI;
  else
    return atan2(sqrt(1/(x*x-1)),fGamma_cm);
}
double Kinematics::GetMaxAngle(int part){
  return GetMaxAngle(fVcm[part]);
}
bool Kinematics::CheckMaxAngle(double angle, int part){
  return angle <= GetMaxAngle(fVcm[part]);
}
double Kinematics::Angle_lab2cm(double vcm, double angle_lab){
  double tan_lab, gtan,x;
  tan_lab = tan(angle_lab);
  gtan = tan_lab*tan_lab*fGamma_cm*fGamma_cm;
  x = fBeta_cm/vcm;

  if(tan_lab>=0){
    //cout << "tan_lab>=0" << endl;
    return acos( (-x*gtan+sqrt( 1+gtan*(1-x*x) ))/(1+gtan) );
  }
  else{
    //cout << "tan_lab<0" << endl;
    return acos( (-x*gtan-sqrt( 1+gtan*(1-x*x) ))/(1+gtan) );
  }
}
double Kinematics::Angle_lab2cminverse(double vcm, double angle_lab, bool upper){
  double tan_lab, gtan,x;
  tan_lab = tan(angle_lab);
  gtan = tan_lab*tan_lab*fGamma_cm*fGamma_cm;
  x = fBeta_cm/vcm;
  //cout << "angle_lab " << angle_lab*180./TMath::Pi() << endl;

  if(upper){
    return acos( (-x*gtan+sqrt( 1+gtan*(1-x*x) ))/(1+gtan) );
  }
  else{
    return acos( (-x*gtan-sqrt( 1+gtan*(1-x*x) ))/(1+gtan) );
  }
}
void Kinematics::AngleErr_lab2cm(double angle, double &err){
  double angle_lab = angle;
  angle = Angle_lab2cm(fVcm[2], angle_lab);
  err =  fabs(Angle_lab2cm(fVcm[2], angle_lab+err)-Angle_lab2cm(fVcm[2], angle_lab-err))/2.;
  /*
  double tang, tang2, gtang, g2,x;
  tang = tan(angle_lab);
  gtang = tang*tang*fGamma_cm*fGamma_cm;
  tang2 = tang*tang;
  g2= fGamma_cm*fGamma_cm;
  x = fBeta_cm/fVcm[2];
  
  if(tang>=0){
    err *= fabs((-2*x*tang*g2*(1*tang2)+tang*g2*(1-x*x)*(1+tang2)/sqrt(1+gtang*(1-x*x)))/(1-gtang) + 2*(-x*gtang+sqrt(1+gtang*(1-x*x))*(1+tang2)*tang*g2)/(1-gtang)/(1-gtang) )/sqrt(1-((-x*gtang+sqrt( 1+gtang*(1-x*x) ))/(1+gtang))*((-x*gtang+sqrt( 1+gtang*(1-x*x) ))/(1+gtang)));
  }
  else{
    err *= fabs((-2*x*tang*g2*(1*tang2)-tang*g2*(1-x*x)*(1+tang2)/sqrt(1+gtang*(1-x*x)))/(1-gtang) + 2*(-x*gtang-sqrt(1+gtang*(1-x*x))*(1+tang2)*tang*g2)/(1-gtang)/(1-gtang) )/sqrt(1-((-x*gtang-sqrt( 1+gtang*(1-x*x) ))/(1+gtang))*((-x*gtang-sqrt( 1+gtang*(1-x*x) ))/(1+gtang)));
    if(err>1){
      cout << "part1 " << (-2*x*tang*g2*(1*tang2)-tang*g2*(1-x*x)*(1+tang2)/sqrt(1+gtang*(1-x*x)))/(1-gtang) << endl;
      cout << "part2 " << 2*(-x*gtang-sqrt(1+gtang*(1-x*x))*(1+tang2)*tang*g2)/(1-gtang)/(1-gtang) << endl;
    }
  }
  */
}
double Kinematics::Angle_cm2lab(double vcm, double angle_cm){
  double x;
  x = fBeta_cm/vcm;
  return atan2(sin(angle_cm),fGamma_cm*(cos(angle_cm)+x));
  /*
  cout << " old " << atan2(sin(angle_cm),fGamma_cm*(cos(angle_cm)+x)) << " x " << x << endl;
  double gam2 = fM[0]*fM[2]/fM[1]/fM[3]*fTCm_i/fTCm_f;
  gam2 = sqrt(gam2);
  double gam3 = fM[0]*fM[3]/fM[1]/fM[2]*fTCm_i/fTCm_f;
  gam3 = sqrt(gam3);
  double y2=1.+gam2*gam2+2.*gam2*cos(angle_cm);
  double y3=1.+gam3*gam3+2.*gam3*cos(PI-angle_cm);
  cout << "y2 " << y2 << " y3 " << y3 << "\n";
  if(asin(sin(PI-angle_cm)/sqrt(y2)) < PI-acos(-gam2))
    cout << " new " << PI - asin(sin(PI-angle_cm)/sqrt(y2)) << " = asin("<<sin(PI-angle_cm)<<"/"<<sqrt(y2)<<")" << endl;
  else 
    cout << " new " << asin(sin(PI-angle_cm)/sqrt(y2)) << " = asin("<<sin(PI-angle_cm)<<"/"<<sqrt(y2)<<")" << endl;
  */

}
//x = sqrt(fM[0]*fM[3]/fM[1]/fM[2]*fTCm_i/fTCm_f);
//cout << "thorsten\t" << asin(sin(angle_cm)/sqrt(1+x*x+2*x*cos(angle_cm)))*180./PI << endl;
//cout << "ich\t" << atan2(sin(angle_cm),fGamma_cm*(cos(angle_cm)+x))*180./PI << endl;  
TSpline3* Kinematics::labvscm(double thmin, double thmax, double size, int part){
  double* cm = new double[(int)((thmax-thmin)/size)+1];
  double* lab = new double[(int)((thmax-thmin)/size)+1];
  int nr =0;
  for(int i=0;i<((thmax-thmin)/size);i++){
    cm[nr] = i;
    lab[nr] = Angle_cm2lab(fVcm[part],cm[nr]*deg2rad)*rad2deg;
    if(lab[nr]>0.01 && lab[nr]<179.99)
      nr++;
  }
  TGraph* graph = new TGraph(nr, cm, lab);
  TSpline3* spline = new TSpline3("Th_cmvslab",graph);
  delete graph;
  delete[] lab;
  delete[] cm;
  return spline;
}
TSpline3* Kinematics::cmvslab(double thmin, double thmax, double size, int part){
  double* cm = new double[(int)((thmax-thmin)/size)+1];
  double* lab = new double[(int)((thmax-thmin)/size)+1];
  int nr =0;
  for(int i=0;i<((thmax-thmin)/size);i++){
    cm[nr] = i;
    lab[nr] = Angle_cm2lab(fVcm[part],cm[nr]*deg2rad)*rad2deg;
    if(lab[nr]>0.01 && lab[nr]<179.99)
      nr++;
  }
  TGraph* graph = new TGraph(nr, lab, cm);
  TSpline3* spline = new TSpline3("Th_cmvslab",graph);
  delete graph;
  delete[] lab;
  delete[] cm;
  return spline;
}

double Kinematics::Sigma_cm2lab(double angle_cm, double sigma_cm){
  double gam2 = fM[0]*fM[2]/fM[1]/fM[3]*fTCm_i/fTCm_f;
  gam2 = sqrt(gam2);
  double wurzel=1.+gam2*gam2+2.*gam2*cos(PI-angle_cm);
  wurzel = sqrt(wurzel);
  return sigma_cm*(1/fGamma_cm*wurzel*wurzel*wurzel/(1+gam2*cos(PI-angle_cm)));
}
/*
double Kinematics::Sigma_cm2lab(double angle_cm, double sigma_cm){
  //old
  double x;
  x = fBeta_cm/fVcm[2];
  double wurzel;
  wurzel = sqrt(sin(PI-angle_cm)*sin(PI-angle_cm) + fGamma_cm*fGamma_cm*(cos(PI-angle_cm)+x)*(cos(PI-angle_cm)+x));
  return sigma_cm*(1/fGamma_cm*wurzel*wurzel*wurzel/(1+x*cos(PI-angle_cm)));  
}
*/
double Kinematics::Sigma_lab2cm(double angle_cm, double sigma_lab){
  double gam2 = fM[0]*fM[2]/fM[1]/fM[3]*fTCm_i/fTCm_f;
  gam2 = sqrt(gam2);
  //test
  //double x;
  //x = fBeta_cm/fVcm[2];
  //cout << "x " << x << " gam2 "<< gam2 << " cm "<< angle_cm << " pi - cm "<< PI-angle_cm << endl;
  //test
  
  //double angle_cm = Angle_lab2cm(fVcm[2], angle_lab);
  //orig
  double wurzel=1.+gam2*gam2+2.*gam2*cos(PI-angle_cm);
  wurzel = sqrt(wurzel);
  return sigma_lab/(1/fGamma_cm*wurzel*wurzel*wurzel/(1+gam2*cos(PI-angle_cm)));  
  //test no pi-
  //double wurzel=1.+gam2*gam2+2.*gam2*cos(angle_cm);
  //wurzel = sqrt(wurzel);
  //return sigma_lab/(1/fGamma_cm*wurzel*wurzel*wurzel/(1+gam2*cos(angle_cm)));  
}
void Kinematics::SigmaErr_lab2cm(double angle, double err, double &sigma, double &errsigma){
  double g = fM[0]*fM[2]/fM[1]/fM[3]*fTCm_i/fTCm_f;
  g = sqrt(g);
  double w=1.+g*g+2.*g*cos(PI-angle);
  w = sqrt(w);
  errsigma = fGamma_cm/pow(w,1.5) * sqrt( pow(sigma*g*sin(PI-angle)*(-2+g*g-g*cos(PI-angle))/w * err,2) + pow((1+g*cos(PI-angle))*errsigma,2  ) ); 
  //sigma/=(1/fGamma_cm*w*w*w/(1+g*cos(PI-angle)));
}
void Kinematics::Transform2cm(double &angle, double &sigma){
  //double angle_lab = angle;
  angle = PI-Angle_lab2cm(fVcm[2], angle);
  sigma = Sigma_lab2cm(angle, sigma);
  return;
}
void Kinematics::Transform2cm(double &angle, double &errangle, double &sigma, double &errsigma){
  AngleErr_lab2cm(angle, errangle);
  Transform2cm(angle,sigma);
  SigmaErr_lab2cm(angle, errangle, sigma, errsigma);
  return;
}
double Kinematics::RutherfordMilliBarn(double angle_cm){
  return Rutherford(angle_cm);
}
double Kinematics::Rutherford(double angle_cm){
  double a = 0.5*1.43997649*fParticle[0]->GetZ()*fParticle[1]->GetZ()/fTCm_i;
  double b = sin(angle_cm/2.)*sin(angle_cm/2.);
  b=b*b;
  return a*a/b*0.0025;//1b=0.01fm
}
TSpline3* Kinematics::Ruthvscm(double thmin, double thmax, double size){
  double* cross = new double[(int)((thmax-thmin)/size)+1];
  double* angle = new double[(int)((thmax-thmin)/size)+1];
  int number =0;
  for(int i=0;i<((thmax-thmin)/size);i++){
    angle[i]=thmin+i*size;
    if(angle[i]>179.99||angle[i]<0.01)
      break;
    cross[i]=Rutherford(angle[i]*deg2rad);
    number++;

    //cout << angle[i] << "   " << Rutherford(GetThetacm(2))<<endl;
    //cout << setprecision(4) << GetThetacm(3)/deg2rad << "\t" << setprecision(4) << GetThetacm(2)/deg2rad << "\t" << setprecision(4) << GetThetalab(3)/deg2rad << "\t" << setprecision(4) << GetThetalab(2)/deg2rad << "\t" << setprecision(4) << Rutherford(GetThetacm(3)) << "\t" << setprecision(4) << Rutherford(GetThetacm(2)) << endl;  
    //cout << (thmin+i*size)*PI/180. << "  max angle: " << GetMaxAngle(fVcm[2]) << endl;
  }
  TGraph* graph = new TGraph(number, angle, cross);
  TSpline3* spline = new TSpline3("sigmaTh_cm",graph);
  delete graph;
  delete[] angle;
  delete[] cross;
  return spline;
}
TSpline3* Kinematics::Ruthvslab(double thmin, double thmax, double size, int part){
  double* cross = new double[(int)((thmax-thmin)/size)+1];
  double* angle = new double[(int)((thmax-thmin)/size)+1];
  int number =0;
  for(int i=0;i<((thmax-thmin)/size);i++){
    if(part==3||part==2)
      angle[i]=thmin+i*size; //angle[i] is in cm system
    else{
      cout << "error " << endl;
      exit(1);
    }
    if(angle[i]>179.99||angle[i]<0.01)
      break;
    cross[i]=Rutherford(angle[i]*deg2rad);
    number++;
    //cout << " angle cm " << angle[i]; 
    //cout << " cs cm " << cross[i];
    
    cross[i]=Sigma_cm2lab(angle[i]*deg2rad, Rutherford(angle[i]*deg2rad));
    if(part==2){
      angle[i]=180-angle[i];
    }
    angle[i]=Angle_cm2lab(fVcm[part],angle[i]*deg2rad)*rad2deg;
    //cout << fVcm[part] << "fVcm[part]" << endl;
    //cout << "\t\tangle lab " << angle[i]; 
    //cout << " cs lab " << cross[i] << endl;; 

  /* for(int i=0;i<(thmax-thmin)/size;i++){
    Final((thmin+i*size)*PI/180.);
    angle[i]=(thmin+i*size);
    if(CheckMaxAngle(angle[i]*PI/180., part)){
      cout << " angle lab " << angle[i] << " GetThetacm(2) " << GetThetacm(2)*180./PI << " GetThetacm(3) " << GetThetacm(3)*180./PI << endl;
      if(angle[i]>179.999||angle[i]<0.001)
	continue;
      if( fabs(GetMaxAngle(fVcm[part])*180./PI-angle[i]) <size)
	continue;
      
      cross[i]= Sigma_cm2lab(fVcm[part], GetThetacm(2), Rutherford(GetThetacm(2)));
      cout << " cs2 cm " <<  Rutherford(GetThetacm(2));
      cout << " cs3 cm " <<  Rutherford(GetThetacm(3));
      cout << " cs lab " << cross[i] << endl;
      if(cross[i]<1e-9)
	continue;
      if(cross[i]>1e+9)
	continue;
      number++;
    }
    else
      continue;
*/
  }
  TGraph* graph = new TGraph(number, angle, cross);
  TSpline3* spline = new TSpline3("sigmaTh_lab",graph);
  delete graph;
  delete[] angle;
  delete[] cross;
  return spline;
}

double Kinematics::Pcm_em(double e, double m){
  return sqrt(e*e-m*m);
}
double Kinematics::P_tm(double t, double m){
  return sqrt(t*t+2.*t*m);
}
double Kinematics::E_tm(double t, double m){
  return t+m;
}
double Kinematics::T_em(double e, double m){
  return e-m;
}
double Kinematics::betacm_tm(double t, double m){
  return sqrt(t*t+2*t*m)/(t+m);
}
double Kinematics::V_pe(double p, double e){
  return p/e;
}
double Kinematics::E_final(int i){
  return fGamma_cm*(fEcm[i]+fBeta_cm*fPcm[i]);
}
double Kinematics::T_final(int i){
  //cout <<(fGamma_cm-1)*fM[i]<<" + "<<fGamma_cm*fTcm[i]<<" + "<<fGamma_cm*fBeta_cm*fPcm[i]*cos(fThetacm[i])<<endl;
  //cout <<fGamma_cm<<"*"<<fBeta_cm<<"*"<<fPcm[i]<<"*"<<cos(fThetacm[i])
  return (fGamma_cm-1)*fM[i]+fGamma_cm*fTcm[i]+fGamma_cm*fBeta_cm*fPcm[i]*cos(fThetacm[i]); 
}
