#include <fstream>
#include <iostream>
#include <string>
#include <limits>
#include <sys/stat.h>

using namespace std;


namespace
{
    void SWAP_ENDIAN_FLOAT(float & value)
    {
        static char Swapped[4];

        float * Addr = &value;

        Swapped[0]=*((char*)Addr+3);
        Swapped[1]=*((char*)Addr+2);
        Swapped[2]=*((char*)Addr+1);
        Swapped[3]=*((char*)Addr  );

        value = *(reinterpret_cast<float*>(Swapped));
    }
}

typedef long long FILE_POSITION;


int main()
{	
	cout << "\nmaximum file position allowed = " << numeric_limits<FILE_POSITION>::max() << "\n\n";

	FILE_POSITION  position;
	string         fname;
	float          time;

	cout << "\nenter file name to modify (case sensitive) > ";
	cin  >> fname;

	struct stat64 stat_buffer;
	stat64(fname.c_str(),&stat_buffer);
	cout << "file size = " << stat_buffer.st_size << "\n";

	fstream file(fname.c_str(),ios::binary|ios::in|ios::out);

	if (!file)
	{
		cout << "file could not be openned\n";
		return 0;
	}

	cout << "\nenter position in file to change > ";
	cin  >> position;

	file.seekg(position,ios::beg);
	file.read( reinterpret_cast<char*>(&time), sizeof(float) );

	SWAP_ENDIAN_FLOAT(time);

	cout << "\ncurrent time = " << time << "\n\n";

	cout << "enter new time > ";
	cin  >> time;

	SWAP_ENDIAN_FLOAT(time);

	file.seekp(position,ios::beg);

	file.write( reinterpret_cast<char*>(&time), sizeof(float) );

	return 0;
}
