#include "MfixData.h"

#include <iostream>
#include <string>
#include <algorithm>
#include <cctype>

using namespace std;

struct ToUpper
{
    char operator () (char c) const
    {
       return toupper(c);
    }
};

int main()
{
   string run_name;
   float  t0 , t1;

   cout << "Enter run name (without the .RES) > ";
   cin  >> run_name;
   
   cout << "Enter first time to keep in the data set > ";
   cin >> t0;
   
   cout << "Enter last  time to keep in the data set > ";
   cin >> t1;
   
   transform(run_name.begin() , run_name.end() , run_name.begin() , ToUpper());
 
   MfixData data;

   data.SetName(run_name.c_str());
     
   data.ReadRes0();

   data.Time_subset_SPx(t0,t1);
   
   return 0;
}
