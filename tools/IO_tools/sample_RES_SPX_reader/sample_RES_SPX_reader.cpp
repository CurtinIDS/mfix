#include "MfixData.h"

#include <fstream>  // for debugging pursposes
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <string>

using namespace std;

#define for if (0); else for


struct IJK
{
   int i,j,k;
   IJK(int ii , int jj , int kk) : i(ii) , j(jj) , k(kk) {}

   bool operator < (const IJK & rhs) const
   {
      if (i < rhs.i) return true;
      if (i > rhs.i) return false;

      if (j < rhs.j) return true;
      if (j > rhs.j) return false;

      return k < rhs.k;
   }
};

void PrintVariable_bottom_to_top(MfixData & data , int var_index , int time_index , ostream & out ,
                                 bool bPrintIndices = false , bool bPrintFlags = false ,
                                 bool bOneBased = true)
{
   data.GetVariableAtTimestep(var_index,time_index);

   out << "\n\n*****************************************************\n\n";
   out << data.variable_names[var_index] << " at time = " << data.times[time_index] << "\n\n";
   out << "bottom to top\n\n";

   int ijk = 0;

   for (int k=0; k<data.kmax2; ++k)
   {
      if (bOneBased)
         out << "k= " << k+1 << "\n\n";
      else
         out << "k= " << k << "\n\n";

      if (bPrintIndices)
      {
         out << "           ";
         for (int i=0; i<data.imax2; ++i)
         {
             if (bOneBased)
                out << "  i = " << setw(5) << i+1 << "   ";
             else
                out << "  i = " << setw(5) << i << "   ";

             if (bPrintFlags) out << "      ";
         }
         out << "\n";
         out << "           ";
         for (int i=0; i<data.imax2; ++i)
         {
             out << " -----"  << "-----" << "---";
             if (bPrintFlags) out << "      ";
         }
         out << "\n";
      }

      for (int j=0; j<data.jmax2; ++j)
      {
         if (bPrintIndices)                out << "[j=" << setw(5) << j << "]  ";

         {
            if (bOneBased)
               out << "[j=" << setw(5) << j+1 << "]  ";
            else
               out << "[j=" << setw(5) << j << "]  ";
         }
         for (int i=0; i<data.imax2; ++i)
         {
             out << setw(14) << data.var[ijk++];
             if (bPrintFlags) out << "(" << setw(4) << data.FLAG[ijk] << ")";
         }
         out << "\n";
      }
   }
}

void PrintVariable_top_to_bottom(MfixData & data , int var_index , int time_index , ostream & out ,
                                 bool bPrintIndices = false , bool bPrintFlags = false ,
                                 bool bOneBased = true)
{

   static map<IJK,int> ijk_map;

   if (ijk_map.size() == 0)
   {
      int ijk = 0;
      for (int k=0; k<data.kmax2; ++k)
      {
         for (int j=0; j<data.jmax2; ++j)
         {
            for (int i=0; i<data.imax2; ++i)
            {
                ijk_map[IJK(i,j,k)] =ijk;
                ++ijk;
            }
         }
      }
   }

   data.GetVariableAtTimestep(var_index,time_index);


   out << "\n\n*****************************************************\n\n";
   out << data.variable_names[var_index] << " at time = " << data.times[time_index] << "\n";
   out << "top to bottom\n\n";

   for (int k=0; k<data.kmax2; ++k)
   {
      if (bOneBased)
         out << "k= " << k+1 << "\n\n";
      else
         out << "k= " << k << "\n\n";

      if (bPrintIndices)
      {
         out << "           ";
         for (int i=0; i<data.imax2; ++i)
         {
             if (bOneBased)
                out << "  i = " << setw(5) << i+1 << "   ";
             else
                out << "  i = " << setw(5) << i << "   ";

             if (bPrintFlags) out << "      ";
         }
         out << "\n";
         out << "           ";
         for (int i=0; i<data.imax2; ++i)
         {
             out << " -----"  << "-----" << "---";
             if (bPrintFlags) out << "      ";
         }
         out << "\n";
      }


      for (int j=data.jmax2-1; j>=0; --j)
      {
         if (bPrintIndices) 
         {
            if (bOneBased)
               out << "[j=" << setw(5) << j+1 << "]  ";
            else
               out << "[j=" << setw(5) << j << "]  ";
         }
         for (int i=0; i<data.imax2; ++i)
         {
             int ijk = ijk_map[IJK(i,j,k)];
             out << setw(14) <<  data.var[ ijk ];
             if (bPrintFlags) out << "(" << setw(4) << data.FLAG[ijk] << ")";
         }
         out << "\n";
      }
   }



}



int main(int argc , char **argv)
{
   int i;

   if (argc != 2)
   {
      cout << "usage:  mfix_reader.exe RUN_NAME\n";
      return 1;
   }


   MfixData data;

   string run_name;

   run_name = argv[1];

   data.SetName(run_name.c_str());
     
   data.ReadRes0();
   
   data.CreateVariableNames();

   data.GetTimes();
   
   string fname = run_name + "_info.txt";
   ofstream out(fname.c_str());  // for debugging purposes


   // rest is for debugging / information
   
   out << " size : " ;
   out << data.imax2 << " * " << data.jmax2 << " * "
       << data.kmax2 << " = " << data.ijkmax2 << "\n";

   // output variable names

   out << "There are " << data.Get_NVARS() << " variables defined\n\n\t";
   
   for (i=0; i<data.variable_names.size(); ++i)
   {
       out << data.variable_names[i] << "\n\t";
   }

   // output time step data

   out << "\n\n\nThere are " << data.Get_NTIMES() << " times\n\n\t";


   for (i=0; i<data.times.size(); ++i)
   {
       out << "(" << setw(6) << i << ")    " << data.times[i] << "\n\t";
   }


   // output timestep / file position information ...
   

   out << "\n\n";

   char ext [] = "123456789ab";

   for (i=0; i<data.nspx_use; ++i)
   {
       out << data.run_name << ".SP" << ext[i] << "\n";

       map<float,FILE_POSITION> & tmap = data.timeMap2[i];
       map<float,FILE_POSITION>::iterator it;
       for (it=tmap.begin(); it!=tmap.end(); ++it)
       {
           out << "\ttime = " << it->first << " at position = " << it->second << "\n";
       }

   }

   // output ... paraview index variable map to MfixData internal variable index

   out << "\n\n" << setw(10) << "Paraview" 
                 << setw(10) << "MfixData"  << "\t"
                 << setw(10) << "variable*"
                 << "\n\n";

   for (int j=0; j<data.index_map.size(); ++j)
   {
       out << setw(10) << j 
           << setw(10) << data.index_map[j] << "\t"
           << setw(10) << data.variable_names[ data.index_map[j] ]
           << "\n";
   }

   out << "XLENGTH = " << data.xlength << "\n";
   out << "YLENGTH = " << data.ylength << "\n";
   out << "ZLENGTH = " << data.zlength << "\n";

   out << "\n\n****************************************************\n\n";

   int var_index  = 0;
   int time_index = 0;

   bool bPrintIndices = true;
   bool bPrintFlags   = true;
   bool bOneBased     = true;


//   PrintVariable_bottom_to_top(data,var_index,time_index,out);

//   PrintVariable_top_to_bottom(data,var_index,time_index,out,bPrintIndices,bPrintFlags);

   bPrintFlags   = false;
   bOneBased     = false;

//   PrintVariable_top_to_bottom(data,var_index,time_index,out,bPrintIndices,bPrintFlags);

//   var_index = 1;
//   PrintVariable_top_to_bottom(data,var_index,time_index,out);


   return 0;
}
