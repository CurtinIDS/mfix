#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <cstring>
#include <iomanip>

#include <sys/types.h>
#include <sys/stat.h>

using namespace std;

#include "MfixData.h"

// MFIX restart helper functions


namespace
{
   int IMIN1 =  2;
   int JMIN1 =  2;
   int KMIN1 =  2;

   int IMAX = 40;
   int JMAX = 700;
   int KMAX = 60;

   int IMAX1 = IMAX + 1;
   int JMAX1 = JMAX + 1;
   int KMAX1 = KMAX + 1;

   int IMAX2 = IMAX + 2;
   int JMAX2 = JMAX + 2;
   int KMAX2 = KMAX + 2;

   int IJMAX2  = IMAX2 * JMAX2;
   int IJKMAX2 = IMAX2 * JMAX2 * KMAX2;

   int ISTART = 1;
   int JSTART = 1;
   int KSTART = 1;


template <typename T>
void CopyArrays(T * a1 , int imax2  , int jmax2  , int kmax2  ,
                T * a2 , int IMAX   , int JMAX   , int KMAX   ,
                         int ISTART , int JSTART , int KSTART)
{
         int flag_count_old = 0;
         int flag_count_new = 0;

         for (int k=0; k<kmax2; ++k)
         {
            for (int j=0; j<jmax2; ++j)
            {
               for (int i=0; i<imax2; ++i)
               {
                   if (i>=ISTART-1 && i<ISTART+IMAX+1)
                   {
                      if (j>=JSTART-1 && j<JSTART+JMAX+1)
                      {
                         if (k>=KSTART-1 && k<KSTART+KMAX+1)
                         {
                            a2[flag_count_new++] = a1[flag_count_old];
                         }
                      }
                   }
                  ++flag_count_old;
              }
            }
         }
}

    template<typename T1 , typename T2>
    T2 Convert(const T1 & t1)
    {
       T2 t2;
       stringstream ss;
       ss << t1;
       ss >> t2;
       return t2;
    }


    void SWAP_DOUBLE(double & value)
    {
    
        static char Swapped[8];

        double * Addr = &value;

        Swapped[0]=*((char*)Addr+7);
        Swapped[1]=*((char*)Addr+6);
        Swapped[2]=*((char*)Addr+5);
        Swapped[3]=*((char*)Addr+4);
        Swapped[4]=*((char*)Addr+3);
        Swapped[5]=*((char*)Addr+2);
        Swapped[6]=*((char*)Addr+1);
        Swapped[7]=*((char*)Addr  );

        value = *(reinterpret_cast<double*>(Swapped));
    
    }

    void SWAP_FLOAT(float & value)
    {
    
        static char Swapped[4];

        float * Addr = &value;

        Swapped[0]=*((char*)Addr+3);
        Swapped[1]=*((char*)Addr+2);
        Swapped[2]=*((char*)Addr+1);
        Swapped[3]=*((char*)Addr  );

        value = *(reinterpret_cast<float*>(Swapped));
    
    }

    void SWAP_INT(int & value)
    {
    
        static char Swapped[4];

        int * Addr = &value;

        Swapped[0]=*((char*)Addr+3);
        Swapped[1]=*((char*)Addr+2);
        Swapped[2]=*((char*)Addr+1);
        Swapped[3]=*((char*)Addr  );

        value = *(reinterpret_cast<int*>(Swapped));
    
    }

    void GetInt(istream& in, int &val)
    {
        in.read( (char*)&val,sizeof(int));
        SWAP_INT(val);
    }

    void PutInt(ostream& out, int val)
    {
        SWAP_INT(val);
        out.write( (char*)&val,sizeof(int));
    }

    void PutDouble(ostream& out, double val)
    {
        SWAP_DOUBLE(val);
        out.write( (char*)&val,sizeof(double));
    }

    void PutFloat(ostream& out, float val)
    {
        SWAP_FLOAT(val);
        out.write( (char*)&val,sizeof(float));
    }

    void GetDouble(istream& in, double& val)
    {
        in.read( (char*)&val,sizeof(double));
        SWAP_DOUBLE(val);
    }

    bool GetFloat(istream& in, float& val)
    {
        in.read( (char*)&val,sizeof(float));

        if (in.good()) SWAP_FLOAT(val);

        return !in.fail();
    }

    static char ext [] = "123456789AB";

    char buffer[513];

    std::string version;
    double version_number;

    double P_ref;
    double P_scale;
    int    DIM_IC;
    int    DIM_BC;
    int    DIM_C;
    int    DIM_IS;
    int    next_reca;

    double c_e;
    double c_f;
    double phi;
    double phi_w;

    double C_e, C_f, Phi, Phi_w;

    int spx_records_per_timestep;

    map<string,int> varNameToSPX;
    vector<double> C;
   
    double dt , xmin ;

    void SkipBytes(istream& in, int n)
    {
       in.read(buffer,n); // maybe seekg instead
    }





}


  

MfixData::MfixData()
{
   P_ref     = 0;
   P_scale   = 1;
   DIM_IC    = 5;
   DIM_BC    = 5;
   DIM_C     = 5;
   DIM_IS    = 5;
   c_e       = 1.0;
   c_f       = 0.0;
   phi       = 0.0;
   phi_w     = 0.0;
   nspx_use  = 9;
   NScalar   = 0;
   bKepsilon = false;  
   nRR       = 0;
   nTimes    = 0;

   res_option = 0;
}

MfixData::~MfixData() { }



void MfixData::RestartVersionNumber(char* buffer)
{
   string s;

   stringstream ss(buffer);
   ss >> s >> s >> version_number;  // RES = 01.10
   version = buffer;
   version.resize(11);
}


void MfixData::IN_BIN_512(istream& in, double* v, int n)
{
   const int nr = 512/sizeof(double);

   double array[nr];

   int num_records;

   if ( n%nr == 0)
      num_records = n/nr;
   else
      num_records = 1 + n/nr;

   int c = 0;
   for (int i=0; i<num_records; ++i)
   {
       in.read( reinterpret_cast<char*>(&array) , 512 );
       for (int j=0; j<nr; ++j)
       {
           if (c < n) 
           {
              v[c] = array[j];
              SWAP_DOUBLE(v[c]);
              ++c;
           }
       }
   }
}



void MfixData::OUT_BIN_512(ostream& out, double* v, int n)
{
   const int nr = 512/sizeof(double);

   double array[nr];

   int num_records;

   if ( n%nr == 0)
      num_records = n/nr;
   else
      num_records = 1 + n/nr;

   int c = 0;
   for (int i=0; i<num_records; ++i)
   {
       for (int j=0; j<nr; ++j) array[j] = 0;
       for (int j=0; j<nr; ++j)
       {
           if (c < n) 
           {
              array[j] = v[c];
              SWAP_DOUBLE(array[j]);
              ++c;
           }
       }

       out.write( (char*)(&array[0]), 512 );
   }
}

void MfixData::OUT_BIN_512R(ostream& out, float * v, int n)
{
   const int nr = 512/sizeof(float);

   float array[nr];

   int num_records;

   if ( n%nr == 0)
      num_records = n/nr;
   else
      num_records = 1 + n/nr;

   int c = 0;
   for (int i=0; i<num_records; ++i)
   {
       for (int j=0; j<nr; ++j) array[j] = 0;

       for (int j=0; j<nr; ++j)
       {
           if (c < n) 
           {
              array[j] = v[c];

 //             if (i==0 && j==0) cout << "array[0] = " << array[0] << "\n";
              SWAP_FLOAT(array[j]);
              ++c;
           }
       }

       out.write( (char*)(&array[0]), 512 );
   }
}

void MfixData::OUT_BIN_512I(ostream& out, int * v, int n)
{
   const int nr = 512/sizeof(int);

   int array[nr];

   int num_records;

   if ( n%nr == 0)
      num_records = n/nr;
   else
      num_records = 1 + n/nr;

   int c = 0;
   for (int i=0; i<num_records; ++i)
   {
       for (int j=0; j<nr; ++j) array[j] = 0;
       for (int j=0; j<nr; ++j)
       {
           if (c < n) 
           {
              array[j] = v[c];
              SWAP_INT(array[j]);
              ++c;
           }
       }

       out.write( (char*)(&array[0]), 512 );
   }
}

void MfixData::IN_BIN_512R(istream& in, float* v, int n)
{
   const int nr = 512/sizeof(float);
   float array[nr];

   int c = 0;
   for (int i=0; i<spx_records_per_timestep; ++i)
   {
       in.read( reinterpret_cast<char*>(&array) , 512 );
       for (int j=0; j<nr; ++j)
       {
           if (c < n) 
           {
              v[c] = array[j];
              SWAP_FLOAT(v[c]);
              ++c;
           }
       }
   }
}

void MfixData::IN_BIN_512I(istream& in, int* v, int n)
{
   const int nr = 512/sizeof(int);

   int array[nr];

   int num_records;

   if ( n%nr == 0)
      num_records = n/nr;
   else
      num_records = 1 + n/nr;

   int c = 0;
   for (int i=0; i<num_records; ++i)
   {
       in.read( reinterpret_cast<char*>(&array) , 512 );
       for (int j=0; j<nr; ++j)
       {
           if (c < n) 
           {
              v[c] = array[j];
              SWAP_INT(v[c]);
              ++c;
           }
       }
   }
}

void MfixData::ReadRes0()
{
   variable_names.clear();
   varNameToSPX.clear();
   index_map.clear();
   timeMap2.clear();
   times.clear();

   int i,l,n,lc;
    
   // init some variables that are not read in in older versions
   
   int DIMENSION_USR = 5;
   
   string fname = run_name + ".RES";
//   fstream in(fname.c_str(),ios::binary|std::ios::in|std::ios::out);
   ifstream in(fname.c_str(),ios::binary);

   if (!in)
   {
      cout << "could not open file" << endl;
      return;
   }

   ofstream out("SUBSET.RES",ios::binary);

   // version : record 1
   memset(buffer,0,513);
   in.read(buffer,512);
   RestartVersionNumber(buffer);
   
   out.write(buffer,512);
   

   // skip 2 lines : records 2 and 3

   in.read(buffer,512);
   out.write(buffer,512);

   in.read(buffer,512);
   out.write(buffer,512);


   // imin1 etc : record 4
   memset(buffer,0,513);

   FILE_POSITION pos_xlength;
   
   if (version == "RES = 01.00")
   {
      GetInt(in,imin1);  GetInt(in,jmin1);   GetInt(in,kmin1);
      GetInt(in,imax);   GetInt(in,jmax);    GetInt(in,kmax);
      GetInt(in,imax1);  GetInt(in,jmax1);   GetInt(in,kmax1);
      GetInt(in,imax2);  GetInt(in,jmax2);   GetInt(in,kmax2);
      GetInt(in,ijmax2); GetInt(in,ijkmax2); GetInt(in,MMAX);
      GetDouble(in,dt);
      GetDouble(in,xlength);  GetDouble(in,ylength);  GetDouble(in,zlength);
    
      // 15 ints ... 4 doubles = 92 bytes
      SkipBytes(in,420);
   }
   else if(version == "RES = 01.01" || version == "RES = 01.02")
   {
      GetInt(in,imin1);  GetInt(in,jmin1);   GetInt(in,kmin1);
      GetInt(in,imax);   GetInt(in,jmax);    GetInt(in,kmax);
      GetInt(in,imax1);  GetInt(in,jmax1);   GetInt(in,kmax1);
      GetInt(in,imax2);  GetInt(in,jmax2);   GetInt(in,kmax2);
      GetInt(in,ijmax2); GetInt(in,ijkmax2); GetInt(in,MMAX);
      GetInt(in,DIM_IC);    GetInt(in,DIM_BC);
      GetDouble(in,dt);
      GetDouble(in,xlength);  GetDouble(in,ylength);  GetDouble(in,zlength);
    
      // 17 ints ... 4 doubles = 100 bytes
      SkipBytes(in,412);
   }
   else if(version == "RES = 01.03")
   {
      GetInt(in,imin1);  GetInt(in,jmin1);   GetInt(in,kmin1);
      GetInt(in,imax);   GetInt(in,jmax);    GetInt(in,kmax);
      GetInt(in,imax1);  GetInt(in,jmax1);   GetInt(in,kmax1);
      GetInt(in,imax2);  GetInt(in,jmax2);   GetInt(in,kmax2);
      GetInt(in,ijmax2); GetInt(in,ijkmax2); GetInt(in,MMAX);
      GetInt(in,DIM_IC); GetInt(in,DIM_BC);
      GetDouble(in,dt);
      GetDouble(in,xmin);
      GetDouble(in,xlength);  GetDouble(in,ylength);  GetDouble(in,zlength);
    
      // 17 ints ... 5 doubles = 108 bytes
      SkipBytes(in,404);
   }
   else if(version == "RES = 01.04")
   {
      GetInt(in,imin1);  GetInt(in,jmin1);   GetInt(in,kmin1);
      GetInt(in,imax);   GetInt(in,jmax);    GetInt(in,kmax);
      GetInt(in,imax1);  GetInt(in,jmax1);   GetInt(in,kmax1);
      GetInt(in,imax2);  GetInt(in,jmax2);   GetInt(in,kmax2);
      GetInt(in,ijmax2); GetInt(in,ijkmax2); GetInt(in,MMAX);
      GetInt(in,DIM_IC); GetInt(in,DIM_BC);  GetInt(in,DIM_C);
      GetDouble(in,dt);
      GetDouble(in,xmin);
      GetDouble(in,xlength);  GetDouble(in,ylength);  GetDouble(in,zlength);
    
      // 18 ints ... 5 doubles = 112 bytes
      SkipBytes(in,400);
   }
   else if(version == "RES = 01.05")
   {
      GetInt(in,imin1);  GetInt(in,jmin1);   GetInt(in,kmin1);
      GetInt(in,imax);   GetInt(in,jmax);    GetInt(in,kmax);
      GetInt(in,imax1);  GetInt(in,jmax1);   GetInt(in,kmax1);
      GetInt(in,imax2);  GetInt(in,jmax2);   GetInt(in,kmax2);
      GetInt(in,ijmax2); GetInt(in,ijkmax2); GetInt(in,MMAX);
      GetInt(in,DIM_IC); GetInt(in,DIM_BC);  GetInt(in,DIM_C);
      GetInt(in,DIM_IS);
      GetDouble(in,dt);
      GetDouble(in,xmin);
      GetDouble(in,xlength);  GetDouble(in,ylength);  GetDouble(in,zlength);
    
      // 19 ints ... 5 doubles = 116 bytes
      SkipBytes(in,396);
   }
   else
   {
      GetInt(in,imin1);  GetInt(in,jmin1);   GetInt(in,kmin1);
      GetInt(in,imax);   GetInt(in,jmax);    GetInt(in,kmax);
      GetInt(in,imax1);  GetInt(in,jmax1);   GetInt(in,kmax1);
      GetInt(in,imax2);  GetInt(in,jmax2);   GetInt(in,kmax2);
      GetInt(in,ijmax2); GetInt(in,ijkmax2); GetInt(in,MMAX);
      GetInt(in,DIM_IC); GetInt(in,DIM_BC);  GetInt(in,DIM_C);
      GetInt(in,DIM_IS);
      GetDouble(in,dt);
      GetDouble(in,xmin);
      GetDouble(in,xlength);  GetDouble(in,ylength);  GetDouble(in,zlength);
      GetDouble(in,C_e); GetDouble(in,C_f); GetDouble(in,Phi); GetDouble(in,Phi_w);
      // 19 ints ... 9 doubles = 148 bytes
      SkipBytes(in,364);

      PutInt(out,IMIN1);  PutInt(out,JMIN1);   PutInt(out,KMIN1);   // 0 - 11
      PutInt(out,IMAX);   PutInt(out,JMAX);    PutInt(out,KMAX);    // 12 - 23
      PutInt(out,IMAX1);  PutInt(out,JMAX1);   PutInt(out,KMAX1);   // 24 - 35
      PutInt(out,IMAX2);  PutInt(out,JMAX2);   PutInt(out,KMAX2);   // 36 - 47
      PutInt(out,IJMAX2); PutInt(out,IJKMAX2); PutInt(out,MMAX);    // 48 - 59
      PutInt(out,DIM_IC); PutInt(out,DIM_BC);  PutInt(out,DIM_C);   // 60 - 71
      PutInt(out,DIM_IS);                                           // 72 - 75
      PutDouble(out,dt);                                            // 76 - 83
      PutDouble(out,xmin);   // ???                                 // 84 - 91

      pos_xlength = out.tellp();
      PutDouble(out,xlength);  PutDouble(out,ylength);  PutDouble(out,zlength);   // 92 - 115
      PutDouble(out,C_e); PutDouble(out,C_f); PutDouble(out,Phi); PutDouble(out,Phi_w);
      // 19 ints ... 9 doubles = 148 bytes   cout << " zzz : " << in.tellg() << "\n";

      out.write(buffer,364);
   }

   /*
   cout << "imin1 = " << imin1 << "\n";
   cout << "jmin1 = " << jmin1 << "\n";
   cout << "kmin1 = " << kmin1 << "\n";

   cout << "imax  = " << imax  << "\n";
   cout << "jmax  = " << jmax  << "\n";
   cout << "kmax  = " << kmax  << "\n";

   cout << "imax1  = " << imax  << "\n";
   cout << "jmax1  = " << jmax  << "\n";
   cout << "kmax1  = " << kmax  << "\n";

   */

   cout << "imax2 = " << imax2 << "\n";
   cout << "jmax2 = " << jmax2 << "\n";
   cout << "kmax2 = " << kmax2 << "\n";

 //  cout << "ijmax2  = " << ijmax2 << "\n";
 //  cout << "ijkmax2 = " << ijkmax2 << "\n";

 //  cout << "xmin = " << xmin << "\n";
   


   
   const int nr = 512/sizeof(float);

   if ( ijkmax2%nr == 0)
      spx_records_per_timestep = ijkmax2/nr;
   else
      spx_records_per_timestep = 1 + ijkmax2/nr;
      
 
   var.resize(ijkmax2);

   // C , C_name and nmax

   NMAX.resize(MMAX+1);
   for (lc=0; lc<MMAX+1; ++lc) NMAX[lc] = 1;
   
   C.resize(DIM_C);

   if (version_number > 1.04)
   {
      IN_BIN_512 (in,&C[0],DIM_C);
      OUT_BIN_512 (out,&C[0],DIM_C);

      for (lc=0; lc<DIM_C; ++lc) 
      {
         in.read(buffer,512);  // c_name[] 
         out.write(buffer,512);
      }
 
      if (version_number < 1.12)
      {
          IN_BIN_512I(in,&NMAX[0],MMAX+1);
      }
      else
      {
         IN_BIN_512I(in,&NMAX[0],MMAX+1); 
         OUT_BIN_512I(out,&NMAX[0],MMAX+1); 
      }
   }

   if (NMAX[0] > 5000)
   {
      cout << "\n";
      cout << "*******************************************\n";
      cout << "\n";
      cout << "NMAX(0) = " << NMAX[0] << "\n";
      cout << "\n";
      cout << "Is that correct ?\n";
      cout << "If not, see mfix/tools/IO_tools/doc/fix_RES_file.pdf\n";
      cout << "\n";
      cout << "kill program , or hit enter to continue > ";
      cin.ignore(1000,'\n');
      cout << "\ncontinuing ...\n";
      cout << "*******************************************\n";
      cout << "\n";
   }

   DX.resize(imax2);
   DY.resize(jmax2);
   DZ.resize(kmax2);

   double xlength_new = 0.0;
   double ylength_new = 0.0;
   double zlength_new = 0.0;

   IN_BIN_512(in,&DX[0],imax2);
   IN_BIN_512(in,&DY[0],jmax2);
   IN_BIN_512(in,&DZ[0],kmax2);

   for (int i=0; i<IMAX; ++i) xlength_new += DX[ISTART+i];   // ??
   for (int i=0; i<JMAX; ++i)
   {
      ylength_new += DY[JSTART+i];
   }
   for (int i=0; i<KMAX; ++i) zlength_new += DZ[KSTART+i];


//   cout << "xlength_new = " << xlength_new << "\n";
//   cout << "ylength_new = " << ylength_new << "\n";
//   cout << "zlength_new = " << zlength_new << "\n";


   OUT_BIN_512(out,&DX[ISTART-1],IMAX2);	// ??
   OUT_BIN_512(out,&DY[JSTART-1],JMAX2);
   OUT_BIN_512(out,&DZ[KSTART-1],KMAX2);
   
 
   // run_name etc.
   
   memset(units,0,17);
   memset(coordinates,0,17);
   
   in.read(buffer,120);      // run_name , description
   out.write(buffer,120);

   in.read(units,16);        // units
   out.write(units,16);

   in.read(buffer,16);       // run_type
   out.write(buffer,16);

   in.read(coordinates,16);  // coordinates 
   out.write(coordinates,16);

   SkipBytes(in,512-168);
   out.write(buffer,512-168);

   char tmp[17];
   
   memset(tmp,0,17);

   int ic = 0;
   for (i=0; i<17; ++i)
   {
        if (units[i] != ' ') tmp[ic++] = units[i];
   }

   memset(tmp,0,17);

   ic = 0;
   for (i=0; i<17; ++i)
   {
        if (coordinates[i] != ' ') tmp[ic++] = coordinates[i];
   }
   strcpy(coordinates,tmp);

   vector<int>    tmpI;
   vector<double> tmpD;

   if (version_number >= 1.04)
   {
      tmpD.resize(NMAX[0]);
      IN_BIN_512(in,&tmpD[0],NMAX[0]);             // MW_g
      OUT_BIN_512(out,&tmpD[0],NMAX[0]);
      for (i=0; i<MMAX; ++i) 
      {
          in.read(buffer,512);  // MW_s
          out.write(buffer,512);
      }
   }
   in.read(buffer,512);  // D_p etc.
   out.write(buffer,512);

   // read in the "DIM_IC" variables (and ignore ... not used by ani_mfix)
   tmpI.resize(DIM_IC);
   tmpD.resize(DIM_IC);




   IN_BIN_512(in,&tmpD[0],DIM_IC);  // ic_x_w
   OUT_BIN_512(out,&tmpD[0],DIM_IC);  // ic_x_w
   IN_BIN_512(in,&tmpD[0],DIM_IC);  // ic_x_e
   OUT_BIN_512(out,&tmpD[0],DIM_IC);  // ic_x_e
   IN_BIN_512(in,&tmpD[0],DIM_IC);  // ic_y_s
   OUT_BIN_512(out,&tmpD[0],DIM_IC);  // ic_y_s
   IN_BIN_512(in,&tmpD[0],DIM_IC);  // ic_y_n
   OUT_BIN_512(out,&tmpD[0],DIM_IC);  // ic_y_n
   IN_BIN_512(in,&tmpD[0],DIM_IC);  // ic_z_b
   OUT_BIN_512(out,&tmpD[0],DIM_IC);  // ic_z_b
   IN_BIN_512(in,&tmpD[0],DIM_IC);  // ic_z_t ...
   OUT_BIN_512(out,&tmpD[0],DIM_IC);  // ic_z_t ...
   IN_BIN_512I(in,&tmpI[0],DIM_IC);  // ic_i_w
   OUT_BIN_512I(out,&tmpI[0],DIM_IC);  // ic_i_w
   IN_BIN_512I(in,&tmpI[0],DIM_IC);  // ic_i_e
   OUT_BIN_512I(out,&tmpI[0],DIM_IC);  // ic_i_e
   IN_BIN_512I(in,&tmpI[0],DIM_IC);  // ic_j_s
   OUT_BIN_512I(out,&tmpI[0],DIM_IC);  // ic_j_s
   IN_BIN_512I(in,&tmpI[0],DIM_IC);  // ic_j_n
   OUT_BIN_512I(out,&tmpI[0],DIM_IC);  // ic_j_n
   IN_BIN_512I(in,&tmpI[0],DIM_IC);  // ic_k_b
   OUT_BIN_512I(out,&tmpI[0],DIM_IC);  // ic_k_b
   IN_BIN_512I(in,&tmpI[0],DIM_IC);  // ic_k_t ...
   OUT_BIN_512I(out,&tmpI[0],DIM_IC);  // ic_k_t ...
   IN_BIN_512(in,&tmpD[0],DIM_IC);  // ic_ep_g
   OUT_BIN_512(out,&tmpD[0],DIM_IC);  // ic_ep_g
   IN_BIN_512(in,&tmpD[0],DIM_IC);  // ic_p_g
   OUT_BIN_512(out,&tmpD[0],DIM_IC);  // ic_p_g
   IN_BIN_512(in,&tmpD[0],DIM_IC);  // ic_t_g
   OUT_BIN_512(out,&tmpD[0],DIM_IC);  // ic_t_g

   if (version_number < 1.15)
   {
      IN_BIN_512(in,&tmpD[0],DIM_IC);  // ic_t_s(1,1)
      IN_BIN_512(in,&tmpD[0],DIM_IC);  // ic_t_s(1,2) or ic_tmp 
   }

   if (version_number >= 1.04)
   {
      for (int i=0; i<NMAX[0]; ++i) 
      {
        IN_BIN_512(in,&tmpD[0],DIM_IC); // ic_x_g
        OUT_BIN_512(out,&tmpD[0],DIM_IC); // ic_x_g
      }
   }

   IN_BIN_512(in,&tmpD[0],DIM_IC); // ic_u_g
   OUT_BIN_512(out,&tmpD[0],DIM_IC);  // ic_u_g
   IN_BIN_512(in,&tmpD[0],DIM_IC); // ic_v_g
   OUT_BIN_512(out,&tmpD[0],DIM_IC);  // ic_v_g
   IN_BIN_512(in,&tmpD[0],DIM_IC); // ic_w_g
   OUT_BIN_512(out,&tmpD[0],DIM_IC);  // ic_w_g

   for (lc=0; lc<MMAX; ++lc)
   {
      IN_BIN_512(in,&tmpD[0],DIM_IC); // ic_rop_s
      OUT_BIN_512(out,&tmpD[0],DIM_IC);  // ic_w_g
      IN_BIN_512(in,&tmpD[0],DIM_IC); // ic_u_s
      OUT_BIN_512(out,&tmpD[0],DIM_IC);  // ic_w_g
      IN_BIN_512(in,&tmpD[0],DIM_IC); // ic_v_s
      OUT_BIN_512(out,&tmpD[0],DIM_IC);  // ic_w_g
      IN_BIN_512(in,&tmpD[0],DIM_IC); // ic_w_s
      OUT_BIN_512(out,&tmpD[0],DIM_IC);  // ic_w_g

      if (version_number >= 1.15)
      {
         IN_BIN_512(in,&tmpD[0],DIM_IC); // ic_t_s
         OUT_BIN_512(out,&tmpD[0],DIM_IC);  // ic_w_g
      }
      
      if (version_number >= 1.04)
      {
         for (n=0; n<NMAX[lc+1]; ++n) // ...
         {
            IN_BIN_512(in,&tmpD[0],DIM_IC); // ic_x_s
            OUT_BIN_512(out,&tmpD[0],DIM_IC);  // ic_w_g
         }
      }
   }


   // read in the "DIM_BC" variables (and ignore ... not used by ani_mfix)
   tmpI.resize(DIM_BC);
   tmpD.resize(DIM_BC);

   IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_x_w
   OUT_BIN_512(out,&tmpD[0],DIM_BC);  // ic_w_g
   IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_x_e
   OUT_BIN_512(out,&tmpD[0],DIM_BC);  // ic_w_g
   IN_BIN_512(in,&tmpD[0],DIM_BC); // bc y s
   OUT_BIN_512(out,&tmpD[0],DIM_BC);  // ic_w_g
   IN_BIN_512(in,&tmpD[0],DIM_BC); // bc y n
   OUT_BIN_512(out,&tmpD[0],DIM_BC);  // ic_w_g
   IN_BIN_512(in,&tmpD[0],DIM_BC); // bc z b
   OUT_BIN_512(out,&tmpD[0],DIM_BC);  // ic_w_g
   IN_BIN_512(in,&tmpD[0],DIM_BC);  // bc z t
   OUT_BIN_512(out,&tmpD[0],DIM_BC);  // ic_w_g
   IN_BIN_512I(in,&tmpI[0],DIM_BC);  // bc i w
   OUT_BIN_512I(out,&tmpI[0],DIM_BC);  // ic_w_g
   IN_BIN_512I(in,&tmpI[0],DIM_BC); // bc i e
   OUT_BIN_512I(out,&tmpI[0],DIM_BC);  // ic_w_g
   IN_BIN_512I(in,&tmpI[0],DIM_BC); // bc j s
   OUT_BIN_512I(out,&tmpI[0],DIM_BC);  // ic_w_g
   IN_BIN_512I(in,&tmpI[0],DIM_BC); // bc j n
   OUT_BIN_512I(out,&tmpI[0],DIM_BC);  // ic_w_g
   IN_BIN_512I(in,&tmpI[0],DIM_BC); // bc k b
   OUT_BIN_512I(out,&tmpI[0],DIM_BC);  // ic_w_g
   IN_BIN_512I(in,&tmpI[0],DIM_BC); // bc k t
   OUT_BIN_512I(out,&tmpI[0],DIM_BC);  // ic_w_g
   IN_BIN_512(in,&tmpD[0],DIM_BC); // bc ep g
   OUT_BIN_512(out,&tmpD[0],DIM_BC);  // ic_w_g
   IN_BIN_512(in,&tmpD[0],DIM_BC); // bc p g
   OUT_BIN_512(out,&tmpD[0],DIM_BC);  // ic_w_g
   IN_BIN_512(in,&tmpD[0],DIM_BC); // bc t g
   OUT_BIN_512(out,&tmpD[0],DIM_BC);  // ic_w_g

   if (version_number < 1.15)
   {
      IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_t_s(1,1)
      IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_t_s(1,1) or bc_tmp
   }

   if (version_number >= 1.04)
   {
      for (int i=0; i<NMAX[0]; ++i) 
      {
         IN_BIN_512(in ,&tmpD[0],DIM_BC); // bc_x_g
        OUT_BIN_512(out,&tmpD[0],DIM_BC);  // ic_w_g
      }
   }

   IN_BIN_512(in,&tmpD[0],DIM_BC); // bc u g
   OUT_BIN_512(out,&tmpD[0],DIM_BC);  // ic_w_g
   IN_BIN_512(in,&tmpD[0],DIM_BC); // bc v g
   OUT_BIN_512(out,&tmpD[0],DIM_BC);  // ic_w_g
   IN_BIN_512(in,&tmpD[0],DIM_BC); // bc w g
   OUT_BIN_512(out,&tmpD[0],DIM_BC);  // ic_w_g
   IN_BIN_512(in,&tmpD[0],DIM_BC); // bc ro g
   OUT_BIN_512(out,&tmpD[0],DIM_BC);  // ic_w_g
   IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_rop_g
   OUT_BIN_512(out,&tmpD[0],DIM_BC);  // ic_w_g
   IN_BIN_512(in,&tmpD[0],DIM_BC); // bc volflow g
   OUT_BIN_512(out,&tmpD[0],DIM_BC);  // ic_w_g
   IN_BIN_512(in,&tmpD[0],DIM_BC); // bc massflow g
   OUT_BIN_512(out,&tmpD[0],DIM_BC);  // ic_w_g

   for (lc=0; lc<MMAX; ++lc)
   {
      IN_BIN_512(in,&tmpD[0],DIM_BC); // bc rop s
      OUT_BIN_512(out,&tmpD[0],DIM_BC);  // ic_w_g
      IN_BIN_512(in,&tmpD[0],DIM_BC); // bc u s
      OUT_BIN_512(out,&tmpD[0],DIM_BC);  // ic_w_g
      IN_BIN_512(in,&tmpD[0],DIM_BC); // bc v s
      OUT_BIN_512(out,&tmpD[0],DIM_BC);  // ic_w_g
      
      if (version_number >= 1.04)
      {
         IN_BIN_512(in,&tmpD[0],DIM_BC); // bc w s
         OUT_BIN_512(out,&tmpD[0],DIM_BC);  // ic_w_g

         if (version_number >= 1.15)
         {
            IN_BIN_512(in,&tmpD[0],DIM_BC); // bc t s
            OUT_BIN_512(out,&tmpD[0],DIM_BC);  // ic_w_g
         }
         for (n=0; n<NMAX[lc+1]; ++n) // ...
         {      
            IN_BIN_512(in,&tmpD[0],DIM_BC); // bc x s
            OUT_BIN_512(out,&tmpD[0],DIM_BC);  // ic_w_g
         }
      }
      IN_BIN_512(in,&tmpD[0],DIM_BC); // bc volflow s
      OUT_BIN_512(out,&tmpD[0],DIM_BC);  // ic_w_g
      IN_BIN_512(in,&tmpD[0],DIM_BC); // bc massflow s
      OUT_BIN_512(out,&tmpD[0],DIM_BC);  // ic_w_g
   }


   if (version == "RES = 01.00")
      l = 10;
   else
      l = DIM_BC;

   for (lc=0; lc<l; ++lc) 
   {
      in.read(buffer,512); // BC TYPE
      out.write(buffer,512);
   }
   
   FILE_POSITION flag_pos = in.tellg();

   FLAG.resize(ijkmax2);
   IN_BIN_512I(in,&FLAG[0],ijkmax2);


   vector<int> FLAG_new(IMAX2*JMAX2*KMAX2,0);

   int flag_count_old = 0;
   int flag_count_new = 0;

   for (int k=0; k<kmax2; ++k)
   {
       for (int j=0; j<jmax2; ++j)
       {
           for (int i=0; i<imax2; ++i)
           {
               if (i>=ISTART-1 && i<ISTART+IMAX+1)
               {
                  if (j>=JSTART-1 && j<JSTART+JMAX+1)
                  {
                     if (k>=KSTART-1 && k<KSTART+KMAX+1)
                     {
                        FLAG_new[flag_count_new++] = FLAG[flag_count_old];
                     }
                  }
               }
               ++flag_count_old;
           }
       }
   }

   OUT_BIN_512I(out,&FLAG_new[0],IJKMAX2);


   // DIM_IS varibles (not needed by ani_mfix)
   tmpI.resize(DIM_IS);
   tmpD.resize(DIM_IS);

   if (version_number >= 1.04)
   {
      IN_BIN_512(in,&tmpD[0],DIM_IS); // is x w
      OUT_BIN_512(out,&tmpD[0],DIM_IS);  // ic_w_g
      IN_BIN_512(in,&tmpD[0],DIM_IS); // is x e
      OUT_BIN_512(out,&tmpD[0],DIM_IS);  // ic_w_g
      IN_BIN_512(in,&tmpD[0],DIM_IS); // is y s
      OUT_BIN_512(out,&tmpD[0],DIM_IS);  // ic_w_g
      IN_BIN_512(in,&tmpD[0],DIM_IS); // is y n
      OUT_BIN_512(out,&tmpD[0],DIM_IS);  // ic_w_g
      IN_BIN_512(in,&tmpD[0],DIM_IS); // is z b
      OUT_BIN_512(out,&tmpD[0],DIM_IS);  // ic_w_g
      IN_BIN_512(in,&tmpD[0],DIM_IS); // is z t
      OUT_BIN_512(out,&tmpD[0],DIM_IS);  // ic_w_g
      IN_BIN_512I(in,&tmpI[0],DIM_IS); // is i w
      OUT_BIN_512I(out,&tmpI[0],DIM_IS);  // ic_w_g
      IN_BIN_512I(in,&tmpI[0],DIM_IS); // is i e
      OUT_BIN_512I(out,&tmpI[0],DIM_IS);  // ic_w_g
      IN_BIN_512I(in,&tmpI[0],DIM_IS); // is j s
      OUT_BIN_512I(out,&tmpI[0],DIM_IS);  // ic_w_g
      IN_BIN_512I(in,&tmpI[0],DIM_IS); // is j n
      OUT_BIN_512I(out,&tmpI[0],DIM_IS);  // ic_w_g
      IN_BIN_512I(in,&tmpI[0],DIM_IS); // is k b
      OUT_BIN_512I(out,&tmpI[0],DIM_IS);  // ic_w_g
      IN_BIN_512I(in,&tmpI[0],DIM_IS); // is k t
      OUT_BIN_512I(out,&tmpI[0],DIM_IS);  // ic_w_g
      IN_BIN_512(in,&tmpD[0],DIM_IS);  // is_pc(1,1)
      OUT_BIN_512(out,&tmpD[0],DIM_IS);  // ic_w_g
      IN_BIN_512(in,&tmpD[0],DIM_IS);  // is_pc(1,2)
      OUT_BIN_512(out,&tmpD[0],DIM_IS);  // ic_w_g
     
      if (version_number >= 1.07)
      {
         for (l=0; l<MMAX; ++l) 
         {
             IN_BIN_512(in ,&tmpD[0],DIM_IS); // is_vel_s
            OUT_BIN_512(out,&tmpD[0],DIM_IS);  // ic_w_g
         }
      }

      for (lc=0; lc<DIM_IS; ++lc) 
      {
        in.read(buffer,512); // is_type
        out.write(buffer,512);
      }
   }

   if (version_number >= 1.08) 
   {
     in.read(buffer,512);
     out.write(buffer,512);
   }

   if (version_number >= 1.09) 
   {
      in.read(buffer,512);
      out.write(buffer,512);
      
      if (version_number >= 1.5)
      {
         GetInt(in,nspx_use);
         PutInt(out,nspx_use);
         SkipBytes(in,508);
         out.write(buffer,508);
      }

      for (lc=0; lc< nspx_use; ++lc) 
      {
         in.read(buffer,512); // spx_dt
         out.write(buffer,512);
      }
      
      for (lc=0; lc<MMAX+1; ++lc) 
      {
         in.read(buffer,512);    // species_eq
         out.write(buffer,512);
      }

      tmpD.resize(DIMENSION_USR);
      
      IN_BIN_512(in,&tmpD[0],DIMENSION_USR); // usr_dt
      OUT_BIN_512(out,&tmpD[0],DIMENSION_USR); // usr_dt
      IN_BIN_512(in,&tmpD[0],DIMENSION_USR); // usr x w
      OUT_BIN_512(out,&tmpD[0],DIMENSION_USR); // usr_dt
      IN_BIN_512(in,&tmpD[0],DIMENSION_USR); // usr x e
      OUT_BIN_512(out,&tmpD[0],DIMENSION_USR); // usr_dt
      IN_BIN_512(in,&tmpD[0],DIMENSION_USR); // usr y s
      OUT_BIN_512(out,&tmpD[0],DIMENSION_USR); // usr_dt
      IN_BIN_512(in,&tmpD[0],DIMENSION_USR); // usr y n
      OUT_BIN_512(out,&tmpD[0],DIMENSION_USR); // usr_dt
      IN_BIN_512(in,&tmpD[0],DIMENSION_USR); // usr z b
      OUT_BIN_512(out,&tmpD[0],DIMENSION_USR); // usr_dt
      IN_BIN_512(in,&tmpD[0],DIMENSION_USR); // usr z t
      OUT_BIN_512(out,&tmpD[0],DIMENSION_USR); // usr_dt

      for (lc=0; lc<DIMENSION_USR; ++lc) 
      {
         in.read(buffer,512);    // usr_ext etc.
         out.write(buffer,512);
      }


      tmpD.resize(DIM_IC);      
      IN_BIN_512(in,&tmpD[0],DIM_IC); // ic_p_star
      OUT_BIN_512(out,&tmpD[0],DIM_IC); // ic_p_star
      IN_BIN_512(in,&tmpD[0],DIM_IC); // ic_l_scale
      OUT_BIN_512(out,&tmpD[0],DIM_IC); // ic_p_star
      for (lc=0; lc<DIM_IC; ++lc) 
      {
         in.read(buffer,512);    // ic_type
         out.write(buffer,512);
      }

      tmpD.resize(DIM_BC);      
      IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_dt_0
      OUT_BIN_512(out,&tmpD[0],DIM_BC); // bc_dt_0
      IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_jet_g0
      OUT_BIN_512(out,&tmpD[0],DIM_BC); // bc_dt_0
      IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_dt_h
      OUT_BIN_512(out,&tmpD[0],DIM_BC); // bc_dt_0
      IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_jet_gh
      OUT_BIN_512(out,&tmpD[0],DIM_BC); // bc_dt_0
      IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_dt_l
      OUT_BIN_512(out,&tmpD[0],DIM_BC); // bc_dt_0
      IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_jet_gl
      OUT_BIN_512(out,&tmpD[0],DIM_BC); // bc_dt_0
   }

      
   if (version_number >= 1.1)  
   {   
      in.read(buffer,512);  // mu_gmax
      out.write(buffer,512);
   }
   if (version_number >= 1.11) 
   {
      in.read(buffer,512);  // x_ex , model_b
      out.write(buffer,512);
   }
      
   if (version_number >= 1.12)
   {
      in.read(buffer,512);   // p_ref , etc.
      out.write(buffer,512);
      in.read(buffer,512);   // leq_it , leq_method
      out.write(buffer,512);
    
      IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_hw_g
      OUT_BIN_512(out,&tmpD[0],DIM_BC); // bc_hw_g
      
      IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_uw_g
      OUT_BIN_512(out,&tmpD[0],DIM_BC); // bc_hw_g
      IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_vw_g
      OUT_BIN_512(out,&tmpD[0],DIM_BC); // bc_hw_g
      IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_ww_g
      OUT_BIN_512(out,&tmpD[0],DIM_BC); // bc_hw_g
    
      for (lc=0; lc<MMAX; ++lc)
      {
         IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_hw_s
         OUT_BIN_512(out,&tmpD[0],DIM_BC); // bc_hw_g
         IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_uw_s
         OUT_BIN_512(out,&tmpD[0],DIM_BC); // bc_hw_g
         IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_vw_s
         OUT_BIN_512(out,&tmpD[0],DIM_BC); // bc_hw_g
         IN_BIN_512(in,&tmpD[0],DIM_BC); // bc_ww_s
         OUT_BIN_512(out,&tmpD[0],DIM_BC); // bc_hw_g
      }
   }
      
   if (version_number >= 1.13) 
   {
      in.read(buffer,512);    // momentum_x_eq , etc.
      out.write(buffer,512);
   }
   if (version_number >= 1.14) 
   {
      in.read(buffer,512);    // detect_small
      out.write(buffer,512);
   }
      
   if (version_number >= 1.15)
   {
      in.read(buffer,512);    // k_g0 , etc.
      out.write(buffer,512);
     
      tmpD.resize(DIM_IC);   
     
      IN_BIN_512(in,&tmpD[0],DIM_IC); // ic_gama_rg
      OUT_BIN_512(out,&tmpD[0],DIM_IC); // bc_hw_g
      IN_BIN_512(in,&tmpD[0],DIM_IC); // ic_t_rg
      OUT_BIN_512(out,&tmpD[0],DIM_IC); // bc_hw_g
        
      for (lc=0; lc<MMAX; ++lc)
      {
         IN_BIN_512(in,&tmpD[0],DIM_IC); // ic_gama_rs
         OUT_BIN_512(out,&tmpD[0],DIM_IC); // bc_hw_g
         IN_BIN_512(in,&tmpD[0],DIM_IC); // ic_t_rs
         OUT_BIN_512(out,&tmpD[0],DIM_IC); // bc_hw_g
      }
  }
     
  if (version_number >= 1.2) 
  {
     in.read(buffer,512); // norm_g , norm_s
     out.write(buffer,512);
  }
 
  if (version_number >= 1.3)
  {
     GetInt(in,NScalar);
     PutInt(out,NScalar);
     SkipBytes(in,sizeof(double)); // tol_resid_scalar
     out.write(buffer,sizeof(double));

     int DIM_tmp;
     GetInt(in,DIM_tmp);
     PutInt(out,DIM_tmp);
     SkipBytes(in,512-sizeof(double)-2*sizeof(int));
     out.write(buffer,512-sizeof(double)-2*sizeof(int));
    
     tmpI.resize(DIM_tmp);
     IN_BIN_512I(in,&tmpI[0],DIM_tmp);  // Phase4Scalar;
     OUT_BIN_512I(out,&tmpI[0],DIM_tmp);
  }
        
  if (version_number >= 1.5)
  {
     GetInt(in,nRR);
     PutInt(out,nRR);
     SkipBytes(in,508);
     out.write(buffer,508);
  }
     
     
  if (version_number >= 1.5999)
  {
     int tmp;
     GetInt(in,tmp);
     PutInt(out,tmp);
     SkipBytes(in,508);
     out.write(buffer,508);
    
     if (tmp != 0) bKepsilon = true;
  }

  FILE_POSITION end_pos = out.tellp();

  FILE_POSITION nrec3_res = end_pos/(FILE_POSITION)512 + (FILE_POSITION)1;

  int rec3_res = nrec3_res;

  double d_time , d_dt;
  int    nstep;
  GetDouble(in,d_time);
  GetDouble(in,d_dt);
  GetInt(in,nstep);
  in.read(buffer,512-sizeof(double)-sizeof(double)-sizeof(int));

  PutDouble(out,d_time);
  PutDouble(out,d_dt);
  PutInt(out,nstep);
  out.write(buffer,512-sizeof(double)-sizeof(double)-sizeof(int));

  vector<double>  v(IMAX2*JMAX2*KMAX2,0.0);
  vector<double> vi(imax2*jmax2*kmax2,0.0);

  n = IMAX2*JMAX2*KMAX2;
  int ni = imax2 * jmax2 * kmax2;

  IN_BIN_512(in,&vi[0],ni);
  CopyArrays<double>(&vi[0]  , imax2  , jmax2  , kmax2  ,
                     &v[0]   , IMAX   , JMAX   , KMAX   ,
                     ISTART  , JSTART , KSTART);
  OUT_BIN_512(out,&v[0],n); // EP_g

  IN_BIN_512(in,&vi[0],ni);
  CopyArrays<double>(&vi[0]  , imax2  , jmax2  , kmax2  ,
                     &v[0]   , IMAX   , JMAX   , KMAX   ,
                     ISTART  , JSTART , KSTART);
  OUT_BIN_512(out,&v[0],n); // P_g

  IN_BIN_512(in,&vi[0],ni);
  CopyArrays<double>(&vi[0]  , imax2  , jmax2  , kmax2  ,
                     &v[0]   , IMAX   , JMAX   , KMAX   ,
                     ISTART  , JSTART , KSTART);
  OUT_BIN_512(out,&v[0],n); // P_star

  IN_BIN_512(in,&vi[0],ni);
  CopyArrays<double>(&vi[0]  , imax2  , jmax2  , kmax2  ,
                     &v[0]   , IMAX   , JMAX   , KMAX   ,
                     ISTART  , JSTART , KSTART);
  OUT_BIN_512(out,&v[0],n); // RO_g

  IN_BIN_512(in,&vi[0],ni);
  CopyArrays<double>(&vi[0]  , imax2  , jmax2  , kmax2  ,
                     &v[0]   , IMAX   , JMAX   , KMAX   ,
                     ISTART  , JSTART , KSTART);
  OUT_BIN_512(out,&v[0],n); // ROP_g

  IN_BIN_512(in,&vi[0],ni);
  CopyArrays<double>(&vi[0]  , imax2  , jmax2  , kmax2  ,
                     &v[0]   , IMAX   , JMAX   , KMAX   ,
                     ISTART  , JSTART , KSTART);
  OUT_BIN_512(out,&v[0],n); // T_g

  for (int i=0; i<NMAX[0]; ++i)
  {

  IN_BIN_512(in,&vi[0],ni);
  CopyArrays<double>(&vi[0]  , imax2  , jmax2  , kmax2  ,
                     &v[0]   , IMAX   , JMAX   , KMAX   ,
                     ISTART  , JSTART , KSTART);
     OUT_BIN_512(out,&v[0],n); // X_g
  }


  IN_BIN_512(in,&vi[0],ni);
  CopyArrays<double>(&vi[0]  , imax2  , jmax2  , kmax2  ,
                     &v[0]   , IMAX   , JMAX   , KMAX   ,
                     ISTART  , JSTART , KSTART);
  OUT_BIN_512(out,&v[0],n); // U_g

  IN_BIN_512(in,&vi[0],ni);
  CopyArrays<double>(&vi[0]  , imax2  , jmax2  , kmax2  ,
                     &v[0]   , IMAX   , JMAX   , KMAX   ,
                     ISTART  , JSTART , KSTART);
  OUT_BIN_512(out,&v[0],n); // V_g

  IN_BIN_512(in,&vi[0],ni);
  CopyArrays<double>(&vi[0]  , imax2  , jmax2  , kmax2  ,
                     &v[0]   , IMAX   , JMAX   , KMAX   ,
                     ISTART  , JSTART , KSTART);
  OUT_BIN_512(out,&v[0],n); // W_g

  for (int i=0; i<MMAX; ++i)
  {

  IN_BIN_512(in,&vi[0],ni);
  CopyArrays<double>(&vi[0]  , imax2  , jmax2  , kmax2  ,
                     &v[0]   , IMAX   , JMAX   , KMAX   ,
                     ISTART  , JSTART , KSTART);
     OUT_BIN_512(out,&v[0],n); // ROP_s

  IN_BIN_512(in,&vi[0],ni);
  CopyArrays<double>(&vi[0]  , imax2  , jmax2  , kmax2  ,
                     &v[0]   , IMAX   , JMAX   , KMAX   ,
                     ISTART  , JSTART , KSTART);
     OUT_BIN_512(out,&v[0],n); // T_s

  IN_BIN_512(in,&vi[0],ni);
  CopyArrays<double>(&vi[0]  , imax2  , jmax2  , kmax2  ,
                     &v[0]   , IMAX   , JMAX   , KMAX   ,
                     ISTART  , JSTART , KSTART);
     OUT_BIN_512(out,&v[0],n); // U_s

  IN_BIN_512(in,&vi[0],ni);
  CopyArrays<double>(&vi[0]  , imax2  , jmax2  , kmax2  ,
                     &v[0]   , IMAX   , JMAX   , KMAX   ,
                     ISTART  , JSTART , KSTART);
     OUT_BIN_512(out,&v[0],n); // V_s

  IN_BIN_512(in,&vi[0],ni);
  CopyArrays<double>(&vi[0]  , imax2  , jmax2  , kmax2  ,
                     &v[0]   , IMAX   , JMAX   , KMAX   ,
                     ISTART  , JSTART , KSTART);
     OUT_BIN_512(out,&v[0],n); // W_s

  IN_BIN_512(in,&vi[0],ni);
  CopyArrays<double>(&vi[0]  , imax2  , jmax2  , kmax2  ,
                     &v[0]   , IMAX   , JMAX   , KMAX   ,
                     ISTART  , JSTART , KSTART);
     OUT_BIN_512(out,&v[0],n); // THETA_m

     for (int j=1; j<=NMAX[i]; ++j)
     {

  IN_BIN_512(in,&vi[0],ni);
  CopyArrays<double>(&vi[0]  , imax2  , jmax2  , kmax2  ,
                     &v[0]   , IMAX   , JMAX   , KMAX   ,
                     ISTART  , JSTART , KSTART);
        OUT_BIN_512(out,&v[0],n); // X_s
     }
  }


  for (int i=0; i<NScalar; ++i)
  {

  IN_BIN_512(in,&vi[0],ni);
  CopyArrays<double>(&vi[0]  , imax2  , jmax2  , kmax2  ,
                     &v[0]   , IMAX   , JMAX   , KMAX   ,
                     ISTART  , JSTART , KSTART);
     OUT_BIN_512(out,&v[0],n); // W_s
  }


  IN_BIN_512(in,&vi[0],ni);
  CopyArrays<double>(&vi[0]  , imax2  , jmax2  , kmax2  ,
                     &v[0]   , IMAX   , JMAX   , KMAX   ,
                     ISTART  , JSTART , KSTART);
  OUT_BIN_512(out,&v[0],n); // gamma_rg

  IN_BIN_512(in,&vi[0],ni);
  CopyArrays<double>(&vi[0]  , imax2  , jmax2  , kmax2  ,
                     &v[0]   , IMAX   , JMAX   , KMAX   ,
                     ISTART  , JSTART , KSTART);
  OUT_BIN_512(out,&v[0],n); // T_rg
  
  for (int i=0; i<MMAX; ++i)
  {

  IN_BIN_512(in,&vi[0],ni);
  CopyArrays<double>(&vi[0]  , imax2  , jmax2  , kmax2  ,
                     &v[0]   , IMAX   , JMAX   , KMAX   ,
                     ISTART  , JSTART , KSTART);
     OUT_BIN_512(out,&v[0],n); // gamma_rs

  IN_BIN_512(in,&vi[0],ni);
  CopyArrays<double>(&vi[0]  , imax2  , jmax2  , kmax2  ,
                     &v[0]   , IMAX   , JMAX   , KMAX   ,
                     ISTART  , JSTART , KSTART);
     OUT_BIN_512(out,&v[0],n); // T_rs
  }
  
  for (int i=0; i<nRR; ++i)
  {

  IN_BIN_512(in,&vi[0],ni);
  CopyArrays<double>(&vi[0]  , imax2  , jmax2  , kmax2  ,
                     &v[0]   , IMAX   , JMAX   , KMAX   ,
                     ISTART  , JSTART , KSTART);
     OUT_BIN_512(out,&v[0],n); // ReactionRates
  }

  if (bKepsilon)
  {

  IN_BIN_512(in,&vi[0],ni);
  CopyArrays<double>(&vi[0]  , imax2  , jmax2  , kmax2  ,
                     &v[0]   , IMAX   , JMAX   , KMAX   ,
                     ISTART  , JSTART , KSTART);
     OUT_BIN_512(out,&v[0],n); // Eturb_g

  IN_BIN_512(in,&vi[0],ni);
  CopyArrays<double>(&vi[0]  , imax2  , jmax2  , kmax2  ,
                     &v[0]   , IMAX   , JMAX   , KMAX   ,
                     ISTART  , JSTART , KSTART);
     OUT_BIN_512(out,&v[0],n); // Kturb_g
  }


  out.seekp(512*2,ios::beg);
  PutInt(out,rec3_res);




  out.seekp(pos_xlength,ios::beg);
  PutDouble(out,xlength_new);
  PutDouble(out,ylength_new);
  PutDouble(out,zlength_new);
}


void MfixData::ReadTimeValues(ifstream & in , FILE_POSITION offset , int spxNum)
{
    in.clear();
    in.seekg( 3*512, ios::beg ); // first time

    float time;

    map<float,FILE_POSITION> tmap2;

    while (in.read( (char*)&time,sizeof(float) ) )
    {
        SWAP_FLOAT(time);

        FILE_POSITION pos = in.tellg();
        pos -= sizeof(float);

        tmap2.insert( make_pair(time,pos) );
        times.insert(time);

        in.seekg(offset,ios::cur);
    }

    timeMap2[spxNum] = tmap2;
}


void MfixData::GetTimes()
{
    timeMap2.resize(nspx_use);

    for (int i=0; i<nspx_use; ++i)
    {
        int nvars = 0;

        string fname = run_name + ".SP" + ext[i];

        ifstream in(fname.c_str(),ios::binary);

        if (in) // file exists
        {

            switch (i+1)
            {

            case 1: nvars = 1; break;

            case 2: nvars = 2; break; 

            case 3: nvars = 3; break;

            case 4: nvars = 3*MMAX; break;

            case 5: nvars = MMAX; break;

            case 6:
                {
                    if (version_number <= 1.15)
                        nvars = 3;
                    else
                        nvars = MMAX + 1;

                    break;
                }

            case 7:
                {
                    nvars = NMAX[0];

                    for (int m=1; m<=MMAX; ++m) // bug
                    {
                        nvars += NMAX[m];
                    }

                    break;
                }


            case 8: nvars = MMAX; break;


            case 9: nvars = NScalar; break;

            case 10: nvars = nRR; break; 

            case 11:
                {
                    if (bKepsilon) nvars = 2;

                    break;
                }


            default:
                {
                    cout << "unknown SPx file : " << i << "\n";
                    break;
                }
            }

            if (nvars > 0) 
            {
                FILE_POSITION offset = (FILE_POSITION)512 - (FILE_POSITION)sizeof(float) + 
                                       (FILE_POSITION)512 * 
                                     ( (FILE_POSITION)nvars * (FILE_POSITION)spx_records_per_timestep );
                ReadTimeValues( in , offset , i );
            }
        }
    }

    times.Sort(); // added SimpleSet

    nTimes = times.size();

}


void MfixData::CreateVariableNames()
{

    for (int i=0; i<nspx_use; ++i)
    {
        string fname = run_name + ".SP" + ext[i];

        ifstream in(fname.c_str(),ios::binary);

        if (in) // file exists
        {
            switch (i+1)
            {

            case 1:
                {
                    variable_names.push_back("EP_g");
                    varNameToSPX.insert(make_pair(string("EP_g"),1));
                    index_map.push_back(variable_names.size()-1);
                    break;
                }

            case 2:
                {
                    variable_names.push_back("P_g");
                    index_map.push_back(variable_names.size()-1);
                    variable_names.push_back("P_star");
                    index_map.push_back(variable_names.size()-1);
                    varNameToSPX.insert(make_pair(string("P_g"),2));
                    varNameToSPX.insert(make_pair(string("P_star"),2));
                    break;
                }

            case 3:
                {
                    variable_names.push_back("U_g");
                    index_map.push_back(variable_names.size()-1);
                    variable_names.push_back("V_g");
                    variable_names.push_back("W_g");
                    varNameToSPX.insert(make_pair(string("U_g"),3));
                    varNameToSPX.insert(make_pair(string("V_g"),3));
                    varNameToSPX.insert(make_pair(string("W_g"),3));
                    break;
                }


            case 4:
                {
                    string us = "U_s_";
                    string vs = "V_s_";
                    string ws = "W_s_";

                    for (int i=0; i<MMAX; ++i)
                    {
                        variable_names.push_back(us+Convert<int,string>(i+1));
                        index_map.push_back(variable_names.size()-1);
                        variable_names.push_back(vs+Convert<int,string>(i+1));
                        variable_names.push_back(ws+Convert<int,string>(i+1));
                        varNameToSPX.insert(make_pair(us+Convert<int,string>(i+1),4));
                        varNameToSPX.insert(make_pair(vs+Convert<int,string>(i+1),4));
                        varNameToSPX.insert(make_pair(ws+Convert<int,string>(i+1),4));

                    }

                    break;
                }

            case 5:
                {
                    string rops = "ROP_s_";

                    for (int i=0; i<MMAX; ++i)
                    {
                        variable_names.push_back(rops+Convert<int,string>(i+1));
                        index_map.push_back(variable_names.size()-1);
                        varNameToSPX.insert(make_pair(rops+Convert<int,string>(i+1) , 5));
                    }

                    break;
                }

            case 6:
                {
                    variable_names.push_back("T_g");
                    index_map.push_back(variable_names.size()-1);
                    varNameToSPX.insert(make_pair(string("T_g"),6));

                    if (version_number <= 1.15)
                    {
                        variable_names.push_back("T_s_1");
                        index_map.push_back(variable_names.size()-1);
                        varNameToSPX.insert(make_pair(string("T_s_1"), 6));

                        if (MMAX > 1) 
                            variable_names.push_back("T_s_2");
                        else
                            variable_names.push_back("T_s_2_not_used");

                        index_map.push_back(variable_names.size()-1);
                        varNameToSPX.insert(make_pair(string("T_s_2"),6));
                    }
                    else
                    {
                        for (int i=0; i<MMAX; ++i)
                        {
                            string ts = "T_s_";

                          //  for (int i=0; i<MMAX; ++i)  // 6-14-2010
                          //  {
                                variable_names.push_back(ts+Convert<int,string>(i+1));
                                index_map.push_back(variable_names.size()-1);
                                varNameToSPX.insert(make_pair(ts+Convert<int,string>(i+1) , 6));
                          //  }
                        }
                    }

                    break;
                }

            case 7:
                {
                    string var = "X_g_";

                    for (int i=0; i<NMAX[0]; ++i)
                    {
                        variable_names.push_back(var+Convert<int,string>(i+1));
                        index_map.push_back(variable_names.size()-1);
                        varNameToSPX.insert(make_pair(var+Convert<int,string>(i+1) , 7));
                    }

                    var = "X_s_";

                    for (int m=0; m<MMAX; ++m)
                    {
                        for (int i=0; i<NMAX[m+1]; ++i)
                        {
                            string name = var + Convert<int,string>(m+1) + "_"
                                              + Convert<int,string>(i+1);

                            variable_names.push_back(name);
                            index_map.push_back(variable_names.size()-1);
                            varNameToSPX.insert(make_pair(name, 7));
                        }
                    }


                    break;
                }


            case 8:
                {
                    string var = "Theta_m_";

                    for (int i=0; i<MMAX; ++i)
                    {
                        variable_names.push_back(var+Convert<int,string>(i+1));
                        index_map.push_back(variable_names.size()-1);
                        varNameToSPX.insert(make_pair(var+Convert<int,string>(i+1),8));
                    }

                    break;
                }


            case 9:
                {
                    string var = "Scalar_";

                    for (int i=0; i<NScalar; ++i)
                    {
                        variable_names.push_back(var+Convert<int,string>(i+1));
                        index_map.push_back(variable_names.size()-1);
                        varNameToSPX.insert(make_pair(var+Convert<int,string>(i+1), 9));
                    }

                    break;
                }


            case 10:
                {
                    string var = "RRates_";

                    for (int i=0; i<nRR; ++i)
                    {
                        variable_names.push_back(var+Convert<int,string>(i+1));
                        index_map.push_back(variable_names.size()-1);
                        varNameToSPX.insert(make_pair(var+Convert<int,string>(i+1), 10));
                    }

                    break;
                }

            case 11:
                {
                    if (bKepsilon)
                    {
                        variable_names.push_back("k_turb_g");
                        index_map.push_back(variable_names.size()-1);
                        variable_names.push_back("e_turb_g");
                        index_map.push_back(variable_names.size()-1);
                        varNameToSPX.insert(make_pair(string("k_turb_g"),11));
                        varNameToSPX.insert(make_pair(string("e_turb_g"),11));
                    }

                    break;
                }


            default:
                {
                    cout << "unknown SPx file : " << i << "\n";
                    break;
                }


            }


        }
    }

}


void MfixData::GetVariableAtTimestep(int vari , int tstep)
{
    // This routine opens and closes the file for each request.
    // Maybe keep all SPX files open, and just perform relative
    // moves to get to the correct location in the file

    // get filename that vaiable # vari is located in

    string & vname = variable_names[vari];

    int spx = varNameToSPX[vname];

    string fname = run_name + ".SP" + ext[spx-1];

    // spx_timeMap2 ... map of times and file position for the
    //                 specified SPX file

    map<float,FILE_POSITION> & spx_timeMap2 = timeMap2[spx-1];

    // find where the requested time  , times[tstep] , is located in the map

    map<float,FILE_POSITION>::iterator mit = spx_timeMap2.lower_bound(times[tstep]);

    if (mit == spx_timeMap2.end()) // time is greater than last time in file
    {
        mit = --spx_timeMap2.end(); // point to last element
    }
    else
    {
        // mit either points to correct time if the time exists
        // in the file , or the next larger time. For now, I want
        // to return the first time before the requested time
        // if it is not in the file 

        // could do time interpolation instead
        // could use the closest time instead

        if (mit->first > times[tstep] && mit->second != 0) --mit;
    }

    int file_position = mit->second;

    int nskip = 0;  // number of variables to skip to get to desired variable
                    // example : request : U_g ... nskip = 0
                    // example : request : V_g ... nskip = 1
                    // example : request : W_g ... nskip = 2

    switch (spx)
    {

    case 1: nskip = 0; break;

    case 2:
        {
            if (vname == "P_g")
                nskip = 0;
            else
                nskip = 1;

            break;
        }

    case 3:
        {
            if (vname == "U_g")
                nskip = 0;
            else if (vname == "V_g")
                nskip = 1;
            else
                nskip = 2;

            break;
        }

    case 4:
        {
            // form of variable : "A_s_b"
            // where "A" = U / V / W
            // and "b" = solid phase index ( 1 - based )

            // examples : U_s_1 , W_s_2

            int m = Convert<string,int>(vname.substr(4));

            nskip = 3*(m-1);
            
            if (vname[0] == 'V') ++nskip;
            if (vname[0] == 'W') nskip += 2;

            break;
        }

    case 5:
        {
            // form of variable : "ROP_s_b"
            // "b" = solid phase index ( 1 - based )

            int m = Convert<string,int>(vname.substr(6));
            nskip = m - 1;

            break;
        }

    case 6:
        {
            // form of variable : "T_g" or "T_s_b"
            // "b" = solid phase index ( 1 - based )

            if (vname == "T_g")
            {
                nskip = 0;
            }
            else
            {
                int m = Convert<string,int>(vname.substr(4));
                nskip = m; // 1 + (m-1)
            }

            break;
        }

    case 7:
        {
            // form of variable : "X_g_a" or "X_s_b_c"
            // "a" = gas species index ( 1 - based )
            // "b" = solid phase index ( 1 - based )
            // "c" = solid species index ( 1 -based )

            if (vname.substr(0,3) == "X_g")
            {
                int ng = Convert<string,int>(vname.substr(4));
                nskip = ng - 1;
            }
            else
            {

                string::size_type pos1 = vname.find_first_of('_');
                string::size_type pos2 = vname.find_last_of('_');

                string str_m = vname.substr( pos1+3 , pos2-pos1-3 );
                int m = Convert<string,int>(str_m);

                int ns = Convert<string,int>(vname.substr(pos2+1));

                nskip = NMAX[0];

                for (int mm=1; mm<m; ++mm)
                {
                    nskip += NMAX[mm];
                }

                nskip += ns - 1;
            }

            break;
        }

    case 8:
        {
            // form of variable : "Theta_m_b"
            // "b" = solid phase index ( 1 - based )

            int m = Convert<string,int>(vname.substr(8));
            nskip = m - 1;

            break;
        }

    case 9:
        {
            // form of variable : "Scalar_s"
            // "s" = scalar index ( 1 - based )

            int m = Convert<string,int>(vname.substr(7));
            nskip = m - 1;

            break;
        }

    case 10:
        {
            // form of variable : "RRates_r"
            // "r" = reaction rate index ( 1 - based )

            int m = Convert<string,int>(vname.substr(7));

             nskip = m - 1;

            break;
        }

    case 11:
        {
            if (vname == "k_turb_g")
                nskip = 0;
            else
                nskip = 1;

            break;
        }
    } // switch

    // get the requested variable/time data

    FILE_POSITION nBytesSkip = mit->second;

    nBytesSkip += 512; // record with time in it

    nBytesSkip += nskip * (FILE_POSITION)spx_records_per_timestep * (FILE_POSITION)512;

    ifstream in(fname.c_str(),ios::binary);

    in.seekg(nBytesSkip,ios::beg);

    IN_BIN_512R (in,&var[0],ijkmax2);
}

void MfixData::Split_SP1(int SPX_file)
{
//    for (int i=0; i<nspx_use; ++i)
    for (int i=SPX_file; i<SPX_file+1; ++i)   // just do one file
    {
        int nvars = 0;

        string fname = run_name + ".SP" + ext[i];

	string fname_a = run_name + ".SP" + ext[i] + "a";
	string fname_b = run_name + ".SP" + ext[i] + "b";

        ifstream in(fname.c_str(),ios::binary);

        if (in) // file exists
        {

            switch (i+1)
            {

            case 1: nvars = 1; break;

            case 2: nvars = 2; break; 

            case 3: nvars = 3; break;

            case 4: nvars = 3*MMAX; break;

            case 5: nvars = MMAX; break;

            case 6:
                {
                    if (version_number <= 1.15)
                        nvars = 3;
                    else
                        nvars = MMAX + 1;

                    break;
                }

            case 7:
                {
                    nvars = NMAX[0];

                    for (int m=1; m<=MMAX; ++m) // bug
                    {
                        nvars += NMAX[m];
                    }

                    break;
                }


            case 8: nvars = MMAX; break;


            case 9: nvars = NScalar; break;

            case 10: nvars = nRR; break; 

            case 11:
                {
                    if (bKepsilon) nvars = 2;

                    break;
                }


            default:
                {
                    cout << "unknown SPx file : " << i << "\n";
                    break;
                }
            }

            if (nvars > 0) 
            {
//                FILE_POSITION offset = (FILE_POSITION)512 - (FILE_POSITION)sizeof(float) + 
//                                       (FILE_POSITION)512 *
//                                     ( (FILE_POSITION)nvars * (FILE_POSITION)spx_records_per_timestep );
                ReadWriteValues( in , i , nvars , fname_a , fname_b);
            }
        }
    }
}


void MfixData::ReadWriteValues(ifstream & in , int spxNum , int nvars ,
             const string & fname_a , const string & fname_b)
{
    ofstream out1(fname_a.c_str(),ios::binary);
    ofstream out2(fname_b.c_str(),ios::binary);

    unsigned int c = 3;

    in.clear();
    in.seekg( 0, ios::beg ); 

    char buf[512];

    // read & write the first 3 records (header)
    for (int i=0; i<3; ++i)
    {
       in.read( buf , 512 );
       out1.write( buf , 512 );
       out2.write( buf , 512 );
    }

    int rec1 = 4;
    int rec2 = 4;

    bool bSwitched = false;

    float time , time_split , time_start;

    cout << "enter the time that should be used to split the file > ";
    cin  >> time_split;

    cout << "enter the time to start writing to the file > ";
    cin  >> time_start;

    while (in.read( (char*)&time,sizeof(float) ) )
    {

        SWAP_FLOAT(time);

        float time_save = time;

        if (time <= time_split)
        {
           if (time >= time_start) ++rec1;
        }
        else
        {
           if (time >= time_start) ++rec2;
        }

        ++c;

        cout << "-processing time = " << time << "   (record # " << c << ")\n";


        in.read( buf , 512-sizeof(float) );
        SWAP_FLOAT(time);

        if (time_save <= time_split)
        {
           if (time_save >= time_start)
           {
              out1.write( (char*)(&time),sizeof(float) );
              out1.write( buf , 512-sizeof(float) );
           }
        }
        else
        {
           if (!bSwitched)
           {
              cout << "switching to second file\n";
              bSwitched = true;
              out1.close();

              cout << fname_a << " ... updating record 3 ... n = " << rec1 << "\n";
              SWAP_INT(rec1);
              fstream out(fname_a.c_str(),ios::binary|ios::in|ios::out);
              out.seekp(512*2,ios::beg);
              out.write( (char*)(&rec1) , 4);
              out.close();
              SWAP_INT(rec1);

           }

           if (time_save >= time_start)
           {
              out2.write( (char*)(&time),sizeof(float) );
              out2.write( buf , 512-sizeof(float) );
           }
        }
 
        if (time_save >= time_start)
        {
		cout << " in write mode\n";
        	for (int i=0; i<nvars*spx_records_per_timestep; ++i)
        	{
           		++c;
           		in.read( buf , 512 );
           		if (time_save <= time_split)

           		{
             			++rec1;
             			out1.write( buf , 512 );
           		}
           		else
           		{
             			++rec2;
             			out2.write( buf , 512 );
           		}
        	}
        }
        else
        {
             FILE_POSITION pos = (FILE_POSITION)512 * (FILE_POSITION)nvars * 
                                 (FILE_POSITION)spx_records_per_timestep;

             in.seekg( pos , ios::cur );
             c += nvars*spx_records_per_timestep;
        }
    }

    if (!bSwitched)
    {
        out1.close();

        cout << fname_a << " ... updating record 3 ... n = " << rec1 << "\n";
        SWAP_INT(rec1);
        fstream out(fname_a.c_str(),ios::binary|ios::in|ios::out);
        out.seekp(512*2,ios::beg);
        out.write( (char*)(&rec1) , 4);
        out.close();
        SWAP_INT(rec1);
    }

    out2.close();

    cout << fname_b << " ... updating record 3 ... n = " << rec2 << "\n";
    fstream out(fname_b.c_str(),ios::binary|ios::in|ios::out);
    out.seekp(512*2,ios::beg);
    SWAP_INT(rec2);
    out.write( (char*)(&rec2) , 4);
    out.close();

   
   ifstream in1(fname_a.c_str(),ios::binary);
   in1.seekg(2*512,ios::beg);
   in1.read( (char*)(&rec1),4);
   SWAP_INT(rec1);
   cout << "first file ... check ... rec4 = " << rec1 << "\n";

   ifstream in2(fname_b.c_str(),ios::binary);
   in2.seekg(2*512,ios::beg);
   in2.read( (char*)(&rec2),4);
   SWAP_INT(rec2);
   cout << "second file ... check ... rec4 = " << rec2 << "\n";
}

int MfixData::NVARS_SPX_file(int SPX_file)
{
   int nvars = -1;

   switch (SPX_file+1)
  {

            case 1: nvars = 1; break;

            case 2: nvars = 2; break; 

            case 3: nvars = 3; break;

            case 4: nvars = 3*MMAX; break;

            case 5: nvars = MMAX; break;

            case 6:
                {
                    if (version_number <= 1.15)
                        nvars = 3;
                    else
                        nvars = MMAX + 1;

                    break;
                }

            case 7:
                {
                    nvars = NMAX[0];

                    for (int m=1; m<=MMAX; ++m) // bug
                    {
                        nvars += NMAX[m];
                    }

                    break;
                }


            case 8: nvars = MMAX; break;


            case 9: nvars = NScalar; break;

            case 10: nvars = nRR; break; 

            case 11:
                {
                    if (bKepsilon) nvars = 2;

                    break;
                }


            default:
                {
                    cout << "unknown SPx file : " << SPX_file+1 << "\n";
                    break;
                }
   }
   return nvars;
}

float MfixData::GetLastTime(int SPX_file)
{
	string fname = run_name + ".SP" + ext[SPX_file];
        ifstream in(fname.c_str(),ios::binary);

        FILE_POSITION pos = (FILE_POSITION)512 + (FILE_POSITION)512 * 
                            (FILE_POSITION)NVARS_SPX_file(SPX_file) *
                            (FILE_POSITION)spx_records_per_timestep;

        in.seekg(0,ios::end);

        in.seekg(-pos,ios::cur);

        float t;
        GetFloat(in,t);

	return t;
}

float MfixData::GetLastTime(fstream & in , int nvars)
{
        FILE_POSITION pos = (FILE_POSITION)512   + (FILE_POSITION)512 * 
                            (FILE_POSITION)nvars * (FILE_POSITION)spx_records_per_timestep;

	in.clear();

        in.seekg(0,ios::end);

        in.seekg(-pos,ios::cur);

        float t;
        GetFloat(in,t);

	return t;
}

int MfixData::GetRec3_value(fstream & in)
{
	in.clear();

	in.seekg(512*2,ios::beg);

        int rec3;

        GetInt(in,rec3);

	return rec3;
}



void MfixData::CombineSPX(int SPX_file , MfixData & data2)
{
	if ( NVARS_SPX_file(SPX_file)  != data2.NVARS_SPX_file(SPX_file) )
	{
		cout << "number of variables in the SPx file are not equal\n\n";
		return;
	}

	string fname = run_name + ".SP" + ext[SPX_file];
        fstream in(fname.c_str(),ios::binary | ios::in | ios::out);

        if (!in)
	{
		cout << "could not open file : " << fname << "\n";
		return;
	}

	int nvars =  NVARS_SPX_file(SPX_file);

	int   rec3      = GetRec3_value(in);
	float last_time = GetLastTime(in,nvars);

   //     cout << "\n";
  //      cout << "rec3      = " << rec3      << "\n";
        cout << "last time (first file) = " << last_time << "\n";

	// open second file and get to first time record
	string fname2 = data2.run_name + ".SP" + ext[SPX_file];
        ifstream in2(fname2.c_str(),ios::binary | ios::in | ios::out);

	cout << fname  << "\n";
	cout << fname2 << "\n\n\n";


	// skip first time in second file ... it is a repeat
	// of the last time in the first file

        FILE_POSITION pos  = 3*512;  				// header
                  pos += 512;    				// first time
                  pos += 512*nvars*spx_records_per_timestep;	// data for first time

        in2.seekg(pos,ios::beg);

	// process the second file ... appending to the first

        in.seekp(0,ios::end);

	char  buf[512];
	float time2 , time_file;

        while (in2.read(buf,512))
        {
		time2 = *(reinterpret_cast<float*>(buf));
                SWAP_FLOAT(time2);
                time_file = time2 + last_time;
		cout << "processing time = " << setw(15) << time_file 
                                             << " = " 
                                             << setw(15) << last_time
                                             << " + "
                                             << time2
                                             << "\n";

 		SWAP_FLOAT(time_file);

		++rec3;
		in.write(reinterpret_cast<char*>(&time_file),sizeof(float));
		in.write(buf+4,512-sizeof(float));

		for (int i=0; i<nvars*spx_records_per_timestep; ++i)
		{
			in2.read(buf,512);
			in.write(buf,512);
                        ++rec3;
		}

        }

	// write the "next record" value in record 3

        cout << "updating record 3 value ...\n";

	in.seekp(2*512,ios::beg);
	SWAP_INT(rec3);
	in.write( reinterpret_cast<char*>(&rec3),sizeof(int) );
}


void MfixData::Subset(int spx_p)
{
   char buffer[512];
   char ext[] = "123456789ab";

   int spx = spx_p - 1;

   int nvars = NVARS_SPX_file(spx);

   cout << "run_name = " << run_name << "\n";

   string fname = string(run_name) + string(".SP") + ext[spx];

   cout << "fname = " << fname << "\n";

   ifstream in(fname.c_str(),ios::binary);
   if (!in) 
   {
     cout << "Could not open file : " << fname << "\n\n";
     return;
   }

   string fname2 = string("SUBSET.SP") + ext[spx];

   cout << "f2 = " << fname2 << "\n";

   ofstream out(fname2.c_str(),ios::binary);
   if (!out) 
   {
      cout << "out : error\n";
   }

   in.read(buffer,512);
   out.write(buffer,512);

   in.read(buffer,512);
   out.write(buffer,512);

   in.read(buffer,512);
   out.write(buffer,512);

   int nrec3 = 4;

   const int nr = 512/sizeof(float);

   int num_records;

   int n = IMAX2*JMAX2*KMAX2;

   if ( n%nr == 0)
      num_records = n/nr;
   else
      num_records = 1 + n/nr;



   float t;
   int   nrec;

   while (GetFloat(in,t))
   {
      cout << "t = " << t << "\n";
      GetInt(in,nrec);
      in.read(buffer,512-sizeof(int)-sizeof(float));
 
      PutFloat(out,t);
      PutInt(out,nrec);
      out.write(buffer,512-sizeof(int)-sizeof(float));
      ++nrec3;

      for (int i=0; i<nvars; ++i)
      {
         IN_BIN_512R(in,&var[0],ijkmax2);

         vector<float> var_new(IMAX2*JMAX2*KMAX2,0.0f);

         CopyArrays<float>(&var[0]     , imax2  , jmax2  , kmax2  ,
                           &var_new[0] , IMAX   , JMAX   , KMAX   ,
                                         ISTART , JSTART , KSTART);


         OUT_BIN_512R(out,&var_new[0],IJKMAX2);
         nrec3 += num_records;
      }
   }

   cout << "nrec3   = " << nrec3 << "\n";
   cout << "tellp() = " << out.tellp() << "\n";

   int nr_subset = 512 / sizeof(float);

   int spx_records_per_timestep_subset;
   if ( ijkmax2%nr == 0)
      spx_records_per_timestep_subset = IJKMAX2/nr_subset;
   else
      spx_records_per_timestep_subset = 1 + IJKMAX2/nr_subset;

   spx_records_per_timestep_subset *= nvars;

   ++spx_records_per_timestep_subset;   // for time record
   out.seekp(2*512,ios::beg);
   PutInt(out,nrec3);
   PutInt(out,spx_records_per_timestep_subset);
}
void MfixData::SetRange(int is , int ni , int js , int nj , int ks , int nk)
{
   ISTART = is;
   IMAX   = ni;
   JSTART = js;
   JMAX   = nj;
   KSTART = ks;
   KMAX   = nk;

   IMIN1 =  2;
   JMIN1 =  2;
   KMIN1 =  2;

   IMAX1 = IMAX + 1;
   JMAX1 = JMAX + 1;
   KMAX1 = KMAX + 1;

   IMAX2 = IMAX + 2;
   JMAX2 = JMAX + 2;
   KMAX2 = KMAX + 2;

   IJMAX2  = IMAX2 * JMAX2;
   IJKMAX2 = IMAX2 * JMAX2 * KMAX2;

}




