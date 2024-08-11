/*

Copyright or © or Copr. Thomas Prest.

Thomas.Prest@ENS.fr

This software is a computer program which purpose is to provide to the 
research community a proof-of-concept implementation of the identity-based
encryption scheme over NTRU lattices, described in the paper
"Efficient Identity-Based Encryption over NTRU Lattices", of
Léo Ducas, Vadim Lyubashevsky and Thomas Prest, available at
homepages.cwi.nl/~ducas/ , www.di.ens.fr/~lyubash/
and www.di.ens.fr/~prest/ .

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.

*/



#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <gmp.h>
#include <chrono>
#include <fstream>

#include <iostream>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <string>
#include <json/json.h>

#include "params.h"
#include "io.h"
#include "FFT.h"
#include "Sampling.h"
#include "Random.h"
#include "Algebra.h"
#include "Scheme.h"

using namespace std;
using namespace NTL;

//==============================================================================
//==============================================================================
//                                  MAIN
//==============================================================================
//==============================================================================
// Function to process the JSON data
void processJsonData(const std::string& data_str) {
    Json::CharReaderBuilder readerBuilder;
    Json::CharReader* reader = readerBuilder.newCharReader();
    Json::Value root;
    std::string errs;

    bool parsingSuccessful = reader->parse(data_str.c_str(), data_str.c_str() + data_str.size(), &root, &errs);
    delete reader;

    if (!parsingSuccessful) {
        std::cerr << "Failed to parse the JSON data: " << errs << std::endl;
        return;
    }

    for (const auto& voter : root) {
        int voter_id = voter["voter_id"].asInt();
        int secret_share = voter["secret_share"].asInt();
        vector<pair<int, int>> shares;

        for (const auto& share : voter["shares"]) {
            int first = share[0].asInt();
            int second = share[1].asInt();
            shares.emplace_back(first, second);
        }

        // Process the data (example: print it)
        cout << "Voter ID: " << voter_id << ", Secret Share: " << secret_share << endl;
        cout << "Shares: ";
        for (const auto& share : shares) {
            cout << "(" << share.first << ", " << share.second << ") ";
        }
        cout << endl;
    }
}

// Function to save MPKD and MSKD to files
void SaveKeys(const MPK_Data &MPKD, const MSK_Data &MSKD) {
    //  unsigned int i;
    
    std::ofstream mpk_file("mpk_data.dat", std::ios::binary);
    mpk_file << MPKD.MPK;
    mpk_file.close();
    cout << "mpk_data.dat saved" << endl;

    std::ofstream msk_file("msk_data.dat", std::ios::binary);
    for (int i = 0; i < 4; i++) {
        msk_file << MSKD.MSK[i];
    }
    msk_file.close();
    cout << "msk_data.dat saved" << endl;

    // // Save SKid_FFT to a .dat file
    // std::ofstream skFile("SKid_FFT_data.dat", std::ios::binary);
    // if (skFile.is_open())
    // {
    //     skFile.write(reinterpret_cast<const char*>(SKid_FFT), sizeof(SKid_FFT));
    //     skFile.close();
    // }

    // // Save id vector to a .dat file
    // std::ofstream idFile("id_data.dat", std::ios::binary);
    // if (idFile.is_open())
    // {
    //     for (i = 0; i < id.length(); i++)
    //     {
    //         long int temp = conv<long int>(id[i]); // Convert ZZ to long int
    //         idFile.write(reinterpret_cast<const char*>(&temp), sizeof(temp));
    //     }
    //     idFile.close();
    // }
}

// Function to load MPKD and MSKD from files
// void LoadKeys(MPK_Data &MPKD, MSK_Data &MSKD) {
//     std::ifstream mpk_file("mpk_data.dat", std::ios::binary);
//     mpk_file >> MPKD.MPK;
//     // Load other MPK_Data members if needed
//     mpk_file.close();

//     std::ifstream msk_file("msk_data.dat", std::ios::binary);
//     for (int i = 0; i < 4; i++) {
//         msk_file >> MSKD.MSK[i];
//     }
//     // Load other MSK_Data members if needed
//     msk_file.close();
// }

void saveMSKDataToFile(const MSK_Data* data, const string& filename) {
    ofstream outfile(filename, ios::binary);
    outfile.write(reinterpret_cast<const char*>(data), sizeof(MSK_Data));
    outfile.close();
}

void loadMSKDataFromFile(MSK_Data* data, const string& filename) {
    ifstream infile(filename, ios::binary);
    infile.read(reinterpret_cast<char*>(data), sizeof(MSK_Data));
    infile.close();
}

int main(int argc, char* argv[])
{
    string arg1;
    int arg2;
    // if (argc > 3) {
    //     std::cerr << "Usage: " << argv[0] << " <arg1> <arg2>" << " argc=" << argc << endl;
    //     return 1;
    // } else if (argc == 3)
    // {
    //     arg1 = argv[1];
    //     arg2 = std::atoi(argv[2]);

    //     cout << "argc: " << argc << endl;
    //     cout << "Hello: " << arg1 << endl;
    //     cout << "Argument 2: " << arg2 << endl;
    // }
    
    // if (argc != 2) {
    //     std::cerr << "Usage: " << argv[0] << " <json_data>" << std::endl;
    //     return 1;
    // }

    // std::string data_str = argv[1];
    // processJsonData(data_str);

    // cout << "\n=======================================================================\n";
    // cout << "This program is a proof-of concept for efficient IBE over lattices.\n";
    // cout << "It generates a NTRU lattice of dimension 2N and associated modulus q,\n";
    // cout << "and perform benches and tests, for user key extraction and encryption/decryption.";
    // cout << "\n=======================================================================\n\n";

    ZZX MSK[4];
    ZZ_pX phiq, MPK;
    unsigned int i;
    float diff;
    MSK_Data * MSKD = new MSK_Data;
    MPK_Data * MPKD = new MPK_Data;
    clock_t t1, t2;
    const ZZX phi = Cyclo();

    // srand(rdtsc()); // initialisation of rand`
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    srand(seed);

    // if no sys input, do encryption
    if (argc == 1) {
        cout << "N = " << N0 << endl;
        cout << "q = " << q0 << endl;

        ZZ_p::init(q1);
        zz_p::init(q0);

        phiq = conv<ZZ_pX>(phi);
        ZZ_pXModulus PHI(phiq);


        cout << "\n===================================================================\n KEY GENERATION";
        cout << "\n===================================================================\n";
        t1 = clock();
        for(i=0; i<1; i++)
        {
            Keygen(MPK, MSK);
        }

        CompleteMSK(MSKD, MSK);
        CompleteMPK(MPKD, MPK);

        t2 = clock();
        diff = ((float)t2 - (float)t1)/1000000.0F;
        cout << "It took " << diff << " seconds to generate the Master Secret Key" << endl;

        std::ofstream mpk_file("../trustee_cred/mpk_data.dat", std::ios::binary);
        // mpk_file << MPKD->MPK;
        mpk_file << MPKD;
        mpk_file.close();
        // cout << "mpk_data.dat saved" << endl;

        // std::ofstream msk_file("../trustee_cred/msk_data.dat", std::ios::binary);
        // // for (int i = 0; i < 4; i++) {
        // //     msk_file << MSKD->MSK[i];
        // // }
        // msk_file << MSKD;
        // msk_file.close();

        saveMSKDataToFile(MSKD, "../trustee_cred/msk_data.dat");
        cout << "msk_data.dat saved" << endl;

        //==============================================================================
        //Key extraction bench and encryption/decryption bench
        //==============================================================================
        // const unsigned int nb_extrb = 100;
        // const unsigned int nb_crypb = 1000;

        // cout << "\n===================================================================\n RUNNING EXTRACTION BENCH FOR ";
        // cout << nb_extrb << " DIFFERENT IDENTITIES\n===================================================================\n";
        // Extract_Bench(nb_extrb, MSKD);

        // cout << "\n===================================================================\n RUNNING ENCRYPTION BENCH FOR ";
        // cout << nb_crypb << " DIFFERENT MESSAGES\n===================================================================\n";
        // Encrypt_Bench(nb_crypb, MPKD, MSKD);



        ///==============================================================================
        //Key extraction test and encryption/decryption test
        //==============================================================================
        const unsigned int nb_extrt = 100;
        const unsigned int nb_voter = 2;
        const unsigned int nb_shares = 3; // nb of trustees

        // cout << "\n===================================================================\n CHECKING EXTRACTION VALIDITY FOR ";
        // cout << nb_extrt << " DIFFERENT IDENTITIES\n===================================================================\n";
        // Extract_Test(nb_extrt, MSKD);

        cout << "\n===================================================================\n CHECKING ENCRYPTION VALIDITY FOR ";
        cout << nb_extrt << " DIFFERENT MESSAGES\n===================================================================\n";
        Encrypt_Test(nb_voter, nb_shares, MPKD, MSKD);

        free(MSKD);
        free(MPKD);
        return 0;
    } else if (argc == 2 && string(argv[1]) == "dec") {
        // decrypt commitment cipher
        unsigned ci, cn0, cj;
        long int id0[N0], Ciphertext[2][N0];
        long int message[N0], decrypted[N0];

        loadMSKDataFromFile(MSKD, "../trustee_cred/msk_data.dat");
        // cout << "MSKD:" << endl;
        // cout << MSKD->sigma << endl; 
        
        ifstream infile("../trustee_cred/cipher_comm.dat", ios::binary);  // Open file in binary mode
        for (ci = 0; ci < 2; ci++) {
            infile.read(reinterpret_cast<char*>(Ciphertext[ci]), N0 * sizeof(long int));
        }
        infile.close();

        // Check if cipher_comm.dat is retrived correctly.
        // cout << "Cipher_comm:" << endl;
        // for (ci = 0; ci < 2; ci++) {
        //     for(cj=0; cj< N0; cj++) {
        //         cout << Ciphertext[ci][cj] << " ";
        //     } cout << endl;
        // }

    }
    else if (argc == 3 && string(argv[1]) == "gen") {
        cout << "Hello Trustee/Authority" << endl; 
        // Key generation and saving
        // ZZX MSK[4];
        // ZZ_pX MPK;
        vec_ZZ id;
        ZZX SK_id[2];
        CC_t SKid_FFT[N0];

        Keygen(MPK, MSK);
        CompleteMSK(MSKD, MSK);
        CompleteMPK(MPKD, MPK);
        
        TrusteeGen(1, MPKD, MSKD);

        free(MSKD);
        free(MPKD);
        return 0;
    } else if (std::string(argv[1]) == "enc") {
        cout << "Hello voter" << endl;
        MPK_Data *MPKD = nullptr;
        MSK_Data *MSKD = nullptr;

        VoterEncrypt(argv[2]);

        free(MSKD);
        free(MPKD);
        return 0;
    }

}
