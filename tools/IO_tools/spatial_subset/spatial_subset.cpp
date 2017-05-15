#include "MfixData.h"

#include <iostream>
#include <string>
#include <sstream>
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

   cout << "Enter run name > ";
   cin >> run_name1;
   cin.ignore(1000,'\n');

 //  run_name1 = "prop/PROPTEST";



   int istart , imax_new , i2;
   int jstart , jmax_new , j2;
   int kstart , kmax_new , k2;

   string line;

   cout << "Enter    I1   I2 > ";
   getline(cin,line);
   replace (line.begin(),line.end(),',',' ');
   stringstream ss1(line);
   ss1 >> istart >> i2;

   cout << "Enter    J1   J2 > ";
   getline(cin,line);
   replace (line.begin(),line.end(),',',' ');
   stringstream ss2(line);
   ss2 >> jstart >> j2;

   cout << "Enter    K1   K2 > ";
   getline(cin,line);
   replace (line.begin(),line.end(),',',' ');
   stringstream ss3(line);
   ss3 >> kstart >> k2;

   imax_new = i2 - istart - 1;
   jmax_new = j2 - jstart - 1;
   kmax_new = k2 - kstart - 1;

   MfixData data1;

   data1.SetRange(istart,imax_new,jstart,jmax_new,kstart,kmax_new);

   int res_option = 0;

   data1.SetName(run_name1.c_str());
   data1.SetResOption(res_option);

   data1.ReadRes0();

   data1.CreateVariableNames();
   data1.GetTimes();

   for (int i=1; i<=9; ++i)
   {
      data1.Subset(i);
   }

   return 0;
}
