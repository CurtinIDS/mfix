#include "MfixData.h"

#include <iostream>
#include <string>
#include <algorithm>
#include <cctype>

using namespace std;

struct ToUpper
{
    char operator () (char c)
    {
       return toupper(c);
    }
};

int main()
{
   string run_name1 , run_name2;


   cout << "\n\nnote : this code does endian byte-swaps\n\n";

   cout << "enter first  run name (without the .RES) (case sensititive) > ";
   cin  >> run_name1;

   cout << "enter second run name (without the .RES) (case sensititive) > ";
   cin  >> run_name2;

  // run_name1 = "P1/proptest";
 //  run_name2 = "P2/proptest";
 //  run_name1 = "R1/bub01";
 //  run_name2 = "R2/bub01";

//   transform(run_name1.begin() , run_name1.end() , run_name1.begin() , ToUpper());
//   transform(run_name2.begin() , run_name2.end() , run_name2.begin() , ToUpper());

   int res_option = 0;

   MfixData data1 , data2;

   data1.SetName(run_name1.c_str());
   data1.SetResOption(res_option);
   data1.ReadRes0();

   cout << "\n";

   data2.SetName(run_name2.c_str());
   data2.SetResOption(res_option);
   data2.ReadRes0();

   if ( (data1.imax2 != data2.imax2) || 
        (data1.jmax2 != data2.jmax2) ||
        (data1.kmax2 != data2.kmax2)
      )
   {
       cout << "dimensions not the same\n";
       return 0;
   }

   int SPX_file;

   cout << "enter the SPX file to combine (1-9) > ";
   cin  >> SPX_file;

   int add_time_code;

   cout << "add last time of first file to each time in second file ? (0 = no) > ";
   cin >> add_time_code;

   data1.CombineSPX( SPX_file-1 , data2 , add_time_code);

   return 0;
}
