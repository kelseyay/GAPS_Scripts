//Attempting to mess with matrices in root... And sorting them...
#include <string>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <list>
#include <vector>
using namespace std;

//Instead use 15 as the xmax. Define xscore[15] and keep track of how many of each score you're getting.

void  MatrixTest(){

        //std::vector<int> myVector = {1, 2, 5};
        //auto it = myVector.begin() + 2; // Iterator pointing to the position before 5 Okay I see what I did wrong...

        int newel = 3;
        //Oh I've got it, let's save the lrms as a string, the strings are saved in the string vector and we can put that into the .txt file!


        std::vector<std::string> fruits = {"apple", "banana", "orange"};

        // Using auto for cleaner syntax (C++11 and later)
        for (auto it = fruits.begin(); it != fruits.end(); ++it) {
            std::cout << *it << " ";
        }


        std::vector<double> mag = {4,3.1,2};

        // Using auto for cleaner syntax (C++11 and later)

        int tkr = 0;

        for(double val : mag){
            cout << val << endl;
            if(newel < val) tkr++;
        }
        cout << "tkr = " << tkr << endl;
        auto itloc = mag.begin() + tkr;

        mag.insert(itloc,newel);

        for(double val : mag){
            cout << val << endl;
        }



/*
    std::vector<std::string> cars = {"Volvo", "BMW", "Ford", "Mazda"};
    cars.insert(2,"Heyoo");

    for (int i = 0; i < cars.size(); i++) {
      cout << cars[i] << "\n";
    }

    std::list<int> mylist = {10,8,3,2,1};
    double newel = 6;
    for(int i = 0; i < mylist.size();i++){
        //cout << i << endl;
        if(newel < mylist[i]){
            mylist.insert(i, newel);
        }
        }

        cout << mylist[1]  << endl;
*/


}
