#include"SLAE_solve.h"

void filing(std::string filename, double norm) {
	std::ofstream fout;
	fout.open(filename, std::ios::app); //std::ios_base::out | std::ios_base::trunc);
	fout << norm << '\t';
	fout << std::endl;
	fout.close();
	fout.clear();
}

void clear_file(std::string filename) {
	std::ofstream fout;
	fout.open(filename, std::ios_base::trunc);
	fout.close();
	fout.clear();
}