#ifndef MFIXDATA_H_ABC_XYZ_098_123
#define MFIXDATA_H_ABC_XYZ_098_123

#include <iosfwd>
#include <string>
#include <vector>
#include <map>

using namespace std;

typedef long long FILE_POSITION;


struct Pair_string_int
{
    std::string first;
    int         second;

    Pair_string_int() {}
    Pair_string_int(const std::string & t1 , const int & t2) : first(t1) , second(t2) {}

    bool operator < (Pair_string_int & rhs)
    {
        return first < rhs.first;
    }

    bool operator <= (Pair_string_int & rhs)
    {
        return first <= rhs.first;
    }
};

struct Pair_float_int
{
    float  first;
    int    second;

    Pair_float_int() {}
    Pair_float_int(const float & t1 , const int & t2) : first(t1) , second(t2) {}

    bool operator < (Pair_float_int & rhs)
    {
        return first < rhs.first;
    }

    bool operator <= (Pair_float_int & rhs)
    {
        return first <= rhs.first;
    }
};










































class SimpleSet
{
public:

    void insert(const float & value)
    {
        // note: assumes the user calls Sort() after inserting data

        // check to see if value is already in the set 
        for (int i=0; i<size(); ++i)
        {
        if (value == v[i]) return;
        }
        v.push_back(value);
    }

    void Sort()
    {
        // do a bubble sort for now ... inefficient, but there are not that
        // many elements in the time array... can use qsort() safely, 
        // since I am only using SimpleSet for POD. (can probably use
        // qsort for non-POD objects , I have never seen it fail)

        int n = size();

        for (int i=0; i< n-1; i++)
        {
            for (int j=i+1; j< n; j++)
            {
                if ( v[ j ] < v[ i ] )
                {
                    float temp = v [ i ];
                    v[ i ] = v[ j ];
                    v[ j ] = temp;
                }
            }
        }
    }

    

    int size() { return v.size(); }

    void clear() { v.clear(); }


   float & operator[] (int i) { return v[i]; }

   const float & operator[] (int i) const { return v[i]; }

private: // making it public for debugging purposes only

    std::vector<float> v;
};












class MfixData
{
public:

   MfixData();
   ~MfixData();
   void ReadRes0();

   void SetName(const char * name) { run_name = name; }
   void CreateVariableNames();

   void RestartVersionNumber(char* buffer);
   void IN_BIN_512(std::istream& in, double* v, int n);
   void IN_BIN_512I(std::istream& in, int* v, int n);
   void IN_BIN_512R(std::istream& in, float* v, int n);


   void OUT_BIN_512(std::ostream& out, double* v, int n);
   void OUT_BIN_512R(std::ostream& out, float* v, int n);
   void OUT_BIN_512I(std::ostream& out, int* v, int n);


   void GetTimes();
   void ReadTimeValues(std::ifstream & in , FILE_POSITION offset , int spxNum);
   void ReadWriteValues(std::ifstream & in , int spxNum , int nvars ,
                const std::string & s1 , const std::string & s2);

   int Get_NVARS() { return variable_names.size(); }
   int Get_NTIMES() { return nTimes; }

   void GetVariableAtTimestep(int var , int tstep);

   void SetResOption(int option) { res_option = option; }

   void Split_SP1(int SPX_file);

   void CombineSPX(int SPX_file , MfixData & data2);

   int NVARS_SPX_file(int SPX_file);

   float GetLastTime(int SPX_file);
   float GetLastTime(fstream & in , int nvars);

   int GetRec3_value(fstream & in);

   void Subset(int spx);

   void SetRange(int is , int ni , int js , int nj , int ks , int nk);

   // variables

   std::string run_name;

   std::vector<std::string> variable_names;

   int imin1 , jmin1 , kmin1 , imax  , jmax  , kmax;
   int imax1 , jmax1 , kmax1 , imax2 , jmax2 , kmax2;
   int ijmax2 , ijkmax2 , MMAX , nspx_use;
   int res_option;

   double xlength , ylength , zlength;
   
   int  NScalar , nRR;

   int nTimes;

   bool bKepsilon;
   
   char coordinates[17];
   char units[17];

   std::vector<int>    FLAG;

   std::vector<double> DX,DY,DZ;

   std::vector<float> var;
   
   std::vector<int> NMAX;

   // The variables that follow could be removed
   // from the header file. They are here to make
   // debugging easier

   std::vector< std::map<float,FILE_POSITION> >  timeMap2;

   SimpleSet times;

   std::vector<int> index_map;
};

#endif
