// Author = Thomas Davis & Dr. A Mouttura, email = txd283@bham.ac.uk / University of Birmingham
// calculates the clustering from the output from convertToCluster.py. It does not distinguish 
// between Va or Cu. Edit the output to give desired results. Currently, it only displays the
// total amount clusters (>3) and the amount of Cu left (isolated).
// to run, compile it with 'g++ clusterCalculator.cpp'

#include <cstdlib>
#include<cmath>
#include<vector>
#include<iostream>
#include<fstream>

#define MAX_REC_LEN 256

using namespace std;

vector< vector<long double> > atoms;
vector<bool> visited;
vector<int> cluster;

void look_for_neighbours(int atom_id)
{
  //cout << "  Looking for neighbours of atom " << atom_id << endl;
  long double dx, dy, dz, dist_sq;
  visited[atom_id] = 1;
  for (int i = 0; i < atoms.size(); i++) {
    if (visited[i] == 0) {
      dx = atoms[atom_id][0] - atoms[i][0];
      if (abs(dx) > 0.5) dx = dx - copysign(1,dx);
      dy = atoms[atom_id][1] - atoms[i][1];
      if (abs(dy) > 0.5) dy = dy - copysign(1,dy);
      dz = atoms[atom_id][2] - atoms[i][2];
      if (abs(dz) > 0.5) dz = dz - copysign(1,dz);
      dist_sq = dx * dx + dy * dy + dz * dz;
      if ( dist_sq <= 0.0004 ) {
        cluster[cluster.size()-1]++;
        look_for_neighbours(i);
      }
    }
  }
}

int main(int argc, char *argv[])
{
  char input_file_name[MAX_REC_LEN];
  ifstream input;
  vector<long double> pos;
  pos.push_back(0.0);
  pos.push_back(0.0);
  pos.push_back(0.0);
  if (argc == 2) {
    sprintf(input_file_name, "%s", argv[1]);
  } else {
    cout << "ERROR: you must specify a file" << endl;
    exit(1);
  }
  input.open(input_file_name);
  if (!input.is_open()) {
    cout << "ERROR: could not open file" << endl;
    exit(1);
  }
  while(!input.eof()) {
    input >> pos[0] >> pos[1] >> pos[2];
    atoms.push_back(pos);
    visited.push_back(0);
  }
  //cout << "  Found " << atoms.size() << " atoms in this file" << endl;
  for (int i=0; i<atoms.size(); i++) {
    if (visited[i] == 0) {
      cluster.push_back(1);
      look_for_neighbours(i);
    }
  }
  
  // output to command line

  int count_1=0, count_2=0,count_3=0, count_4=0, count_5=0, count_6=0, count_7=0, count_8=0, count_9=0, count_10=0;
  int count_11=0,count_12=0,count_13=0,count_14=0,count_15=0,count_16=0, count_17=0, count_18=0, count_19=0,count_20=0;
  for (int j=0; j<cluster.size(); j++){
	  
	  if(cluster[j] == 1) ++count_1;
	  if(cluster[j] == 2) ++count_2;
	  if(cluster[j] == 3) ++count_3;
	  if(cluster[j] == 4) ++count_4;
	  if(cluster[j] == 5) ++count_5;
	  if(cluster[j] == 6) ++count_6;
	  if(cluster[j] == 7) ++count_7;
	  if(cluster[j] == 8) ++count_8;
	  if(cluster[j] == 9) ++count_9;
	  if(cluster[j] == 10) ++count_10;
	  if(cluster[j] == 11) ++count_11;
	  if(cluster[j] == 12) ++count_12;
	  if(cluster[j] == 13) ++count_13;
	  if(cluster[j] == 14) ++count_14;
	  if(cluster[j] == 15) ++count_15;
	  if(cluster[j] == 16) ++count_16;
	  if(cluster[j] == 17) ++count_17;
	  if(cluster[j] == 18) ++count_18;
	  if(cluster[j] == 19) ++count_19;
	  if(cluster[j] == 20) ++count_20;
  }
  
  int total_above = 0; 
	total_above = count_3 +  count_4 +  count_5 +  count_6 +  count_7 +  count_8 +  count_9 +  count_10 + count_11 + count_12 + count_13 + count_14 + count_15 + count_16 +  count_17 +  count_18 +  count_19 + count_20;
	int reduction = 0;
	  reduction = atoms.size() - count_1 - count_2;
	  
/*	  cout << "C_Size   |  Freq  " << endl << endl;
	  cout << count_1 << endl;
	  cout << count_2 << endl;
	  cout << count_3 << endl;
	  cout << count_4 << endl;
	  cout <<  count_5 << endl;
	  cout <<  count_6 << endl;
	  cout <<  count_7 << endl;
	  cout <<  count_8 << endl;
	  cout <<  count_9 << endl;
	  cout <<  count_10 << endl;
	  cout <<  count_11 << endl;
	  cout <<  count_12 << endl;
	  cout <<  count_13 << endl;
	  cout <<  count_14 << endl;
	  cout <<  count_15 << endl;
	  cout <<  count_16 << endl;
	  cout <<  count_17 << endl;
	  cout <<  count_18 << endl;
	  cout <<  count_19 << endl;
cout <<  count_20 << endl << endl; */
	  cout << "------------------------------------------------------------------------" << endl;
	  cout << "Total amount of Clusters (>3) = " << total_above << "	" << count_1 + count_2 << endl;
	  cout << "Cu left = " << count_1 + count_2  << endl;
	  cout << "------------------------------------------------------------------------" << endl;
/*	  
  cout << "C_Size   |  Freq  " << endl << endl;
  cout << "1          	" << count_1 << endl;
  cout << "2          	" << count_2 << endl;
  cout << "3          	" << count_3 << endl;
  cout << "4          	" << count_4 << endl;
  cout << "5          	" << count_5 << endl;
  cout << "6          	" << count_6 << endl;
  cout << "7          	" << count_7 << endl;
  cout << "8          	" << count_8 << endl;
  cout << "9          	" << count_9 << endl;
  cout << "10          	" << count_10 << endl;
  cout << "11          	" << count_11 << endl;
  cout << "12          	" << count_12 << endl;
  cout << "13          	" << count_13 << endl;
  cout << "14          	" << count_14 << endl;
  cout << "15          	" << count_15 << endl;
  cout << "16          	" << count_16 << endl;
  cout << "17          	" << count_17 << endl;
  cout << "18          	" << count_18 << endl;
  cout << "19          	" << count_19 << endl;
  cout << "20          	" << count_20 << endl << endl;
  cout << "Total amount of Clusters (>3) = " << total_above << endl;
  cout << "Reduction is isolation Cu = " << reduction << endl;
  */
  //for (int j=0; j<cluster.size(); j++){
	  
	  
	  
    //cout << "            " << j << "       " << cluster[j] << endl;
  //}
  return 0;
}




