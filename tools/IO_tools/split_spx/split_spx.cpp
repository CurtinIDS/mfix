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
   string run_name;


   cout << "\n\nnote : this code does endian byte-swaps\n\n";

   cout << "enter run name (without the .RES) > ";
   cin  >> run_name;

   transform(run_name.begin() , run_name.end() , run_name.begin() , ToUpper());

   int res_option = 0;

   /*
   while (res_option<1 || res_option>3)
   {
   	cout << "\n\nHow do you want to modify the RES file ?\n\n";
        cout << "\t1)  modify NMAX array\n";
        cout << "\t2)  set all FLAGs equal to FLUID (one)\n";
        cout << "\t3)  modify the dates in the second record\n";
        cout << "\n\nEnter choice > ";
   	cin >> res_option;
   }
   */
   res_option = 4;

   MfixData data;

   data.SetName(run_name.c_str());

   data.SetResOption(res_option);
     
   data.ReadRes0();

   int SPX_file;

   cout << "enter the SPX file to split (1-9) > ";
   cin  >> SPX_file;

   data.Split_SP1(SPX_file-1);
   
   return 0;
}
