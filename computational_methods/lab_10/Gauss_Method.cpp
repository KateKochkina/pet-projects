#include "Gauss_Method.h"

void clear_file_gauss(std::string filename) {
	std::ofstream fout;
	fout.open(filename, std::ios_base::trunc);
	fout.close();
	fout.clear();
}
