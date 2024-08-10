#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <gmp.h>
#include <fstream>
#include <iostream>
#include <vector>

#include "Sampling.h"
#include "params.h"
#include "FFT.h"
#include "Random.h"
#include "Algebra.h"

using namespace std;
using namespace NTL;

const ZZX phi = Cyclo();


//==============================================================================
//Generates from parameters N and q :
// - a public key : polynomial h
// - a private key : polynomials f,g,F,G
//==============================================================================
void Keygen(ZZ_pX& PublicKey, ZZX* PrivateKey)
{
    ZZ SqNorm;
    ZZX f,g,F,G;

    SqNorm = conv<ZZ>(1.36*q0/2);

    GenerateBasis(f, g, F, G, SqNorm);
    PrivateKey[0] = f;
    PrivateKey[1] = g;
    PrivateKey[2] = F;
    PrivateKey[3] = G;

    for(unsigned int i=0; i<4; i++)
    {
            PrivateKey[i].SetLength(N0);
    }

    PublicKey = Quotient(f, g);
}

//==============================================================================
//Computes the private basis B from private key PrivateKey and parameter N
//==============================================================================
void CompletePrivateKey(mat_ZZ& B, const ZZX * const PrivateKey)
{
    ZZX f,g,F,G;
    f = PrivateKey[0];
    g = PrivateKey[1];
    F = PrivateKey[2];
    G = PrivateKey[3];

    f = -f;
    F = -F;

    B = BasisFromPolynomials(g, f, G, F);
}

void GPV(RR_t * v, const RR_t * const c, const RR_t s, const MSK_Data * const MSKD)
{

    int i;
    unsigned j;
    RR_t ci[2*N0], zi, cip, sip, aux;

    for(j=0; j<2*N0;j++)
    {
        ci[j] = c[j];
    }

    // for(j=0; j<2*N0; j++)
    // {

    // }    
    for(i=2*N0-1; i>=0; i--)
    {
        aux = (MSKD->GS_Norms)[i];
        cip = DotProduct(ci, MSKD->Bstar[i])/(aux*aux);
        sip = s/aux;
        zi = Sample4(cip, sip*PiPrime);

        for(j=0; j<2*N0; j++)
        {
            ci[j] -= zi*(MSKD->B)[i][j];
        }
    }
    // cout << "DEBUG" << endl;
    for(j=0; j<2*N0; j++)
    {
        v[j] = c[j] - ci[j];
    }
}

//==============================================================================
//==============================================================================
//                            MAIN PROGRAMS
//==============================================================================
//==============================================================================

void CompleteMSK(MSK_Data * MSKD, ZZX * MSK)
{
    unsigned int i, j;
    mat_ZZ B0;

    for(i=0; i<4; i++)
    {
        MSKD->PrK[i] = MSK[i];
        ZZXToFFT(MSKD->PrK_fft[i], MSK[i]);
    }

    CompletePrivateKey(B0, MSK);

    for(i=0; i<2*N0; i++)
    {
        for(j=0; j<2*N0; j++)
        {
            MSKD->B[i][j] = ( (RR_t) conv<double>(B0[i][j]) );
        }
    }

    for(i=0; i<1; i++)
    {
        FastMGS(MSKD->Bstar, MSKD->B);
    }

    for(i=0; i<2*N0; i++)
    {
        MSKD->GS_Norms[i] = sqrt( DotProduct(MSKD->Bstar[i], MSKD->Bstar[i]) );
    }

    MSKD->sigma = 2*MSKD->GS_Norms[0];

}

void CompleteMPK(MPK_Data * MPKD, ZZ_pX MPK)
{
    MPKD->h = MPK;
    ZZXToFFT(MPKD->h_FFT, conv<ZZX>(MPK));
}

void IBE_Extract(ZZX SK_id[2], vec_ZZ id, const MSK_Data * const MSKD)
{
    unsigned int i;
    RR_t c[2*N0], sk[2*N0], sigma;
    ZZX f,g,aux;

    f = MSKD -> PrK[0];
    g = MSKD -> PrK[1];
    sigma = MSKD->sigma;
    SK_id[0].SetLength(N0);
    SK_id[1].SetLength(N0);

    for(i=0;i<N0;i++)
    {
        c[i] = ((RR_t) conv<double>(id[i])) ;
        c[i+N0] = 0;
    }
    GPV(sk, c, sigma, MSKD);
    for(i=0; i<N0; i++)
    {
        sk[i] = c[i] - sk[i];
        sk[i+N0] = - sk[i+N0];
    }

    for(i=0; i<N0; i++)
    {
        SK_id[0][i] = sk[i];
        SK_id[1][i] = sk[i+N0];
    }
}

unsigned long IBE_Verify_Key(const ZZX SK_id[2], const vec_ZZ id, const MSK_Data * const MSKD)
{
    unsigned int i;
    ZZX f,g,t,aux;

    f = MSKD -> PrK[0];
    g = MSKD -> PrK[1];
    
    t = conv<ZZX>(id);
    aux = ((SK_id[0] - t)*f + g*SK_id[1])%phi;

    for(i=0; i<N0; i++)
    {
        aux[i] %= q1;
    }

    if( IsZero(aux) != 0)
    {
        cout << "The signature (s1,s2) doesn't verify the required equality [ (s1 - t)*f + g*s2 = 0 ] !\nActually, (s1 - t)*f + g*s2 = " << aux << endl << endl;
    }
    return IsZero(aux);
}

void IBE_Encrypt(long C[2][N0], const long m[N0], const long id0[N0], const MPK_Data * const MPKD)
{
    unsigned long i;
    long r[N0], e1[N0], e2[N0];
    CC_t r_FFT[N0], t_FFT[N0], aux1_FFT[N0], aux2_FFT[N0];

    for(i=0; i<N0; i++)
    {
        e1[i] = (rand()%3) - 1;
        e2[i] = (rand()%3) - 1;
        r[i] = (rand()%3) - 1;
    }

    MyIntFFT(r_FFT, r);
    MyIntFFT(t_FFT, id0);

    for(i=0; i<N0; i++)
    {
        aux1_FFT[i] = r_FFT[i]*((MPKD->h_FFT)[i]);
        aux2_FFT[i] = r_FFT[i]*t_FFT[i];
    }

    MyIntReverseFFT(C[0], aux1_FFT);
    MyIntReverseFFT(C[1], aux2_FFT);

    for(i=0; i<N0; i++)
    {
        C[0][i] = (C[0][i] + e1[i]               + q0/2)%q0 - (q0/2);
        C[1][i] = (C[1][i] + e2[i] + (q0/2)*m[i] + q0/2)%q0 - (q0/2);
    } 

}

void IBE_Decrypt(long message[N0], const long C[2][N0], const CC_t * const SKid_FFT)
{
    unsigned int i;
    CC_t c0_FFT[N0], aux_FFT[N0];

    MyIntFFT(c0_FFT, C[0]);

    for(i=0; i<N0; i++)
    {
        aux_FFT[i] = c0_FFT[i]*SKid_FFT[i];
    }

    MyIntReverseFFT(message, aux_FFT);

    for(i=0; i<N0; i++)
    {
        message[i] = C[1][i] - message[i];
        message[i] = ((unsigned long)(message[i] ))%q0;
        message[i] = (message[i] + (q0>>2) )/(q0>>1);
        message[i] %= 2;
    }

}

//==============================================================================
//==============================================================================
//                             BENCHES AND TESTS
//                   FOR EXTRACTION AND ENCRYPTION/DECRYPTION
//==============================================================================
//==============================================================================

void Extract_Bench(const unsigned int nb_extr, MSK_Data * MSKD)
{
    clock_t t1, t2;
    float diff;
    unsigned int i;
    vec_ZZ id;
    ZZX SK_id[2];

    t1 = clock();

    cout << "0%" << flush;
    for(i=0; i<nb_extr; i++)
    {
        id = RandomVector();

        IBE_Extract(SK_id, id, MSKD);
        if((i+1)%(nb_extr/10)==0)
        {
            cout << "..." << (i+1)/(nb_extr/10) << "0%" << flush;
        }
    }

    t2 = clock();
    diff = ((float)t2 - (float)t1)/1000000.0F;
    cout << "\n\nIt took " << diff << " seconds to extract " << nb_extr << " keys." << endl;
    cout << "That's " << (diff/nb_extr)*1000 << " milliseconds per key." << endl << endl;
}

void Encrypt_Bench(const unsigned int nb_cryp, MPK_Data * MPKD, MSK_Data * MSKD)
{
    clock_t t1, t2;
    double diff;
    unsigned int i,j;
    vec_ZZ id;
    ZZX SK_id[2], w;
    CC_t SKid_FFT[N0];
    long int message[N0], decrypted[N0];
    long int identity[N0], Ciphertext[2][N0];


    id = RandomVector();
    IBE_Extract(SK_id, id, MSKD);
    IBE_Verify_Key(SK_id, id, MSKD);
    ZZXToFFT(SKid_FFT, SK_id[1]);
    for(i=0; i<N0; i++)
    {
        identity[i] = conv<long int>(id[i]);
    }

    t1 = clock();

    cout << "0%" << flush;
    for(i=0; i<nb_cryp; i++)
    {

        for(j=0; j<N0; j++)
        {
            message[j] = (rand()%2);
        }

        IBE_Encrypt(Ciphertext, message, identity, MPKD);
        IBE_Decrypt(decrypted, Ciphertext, SKid_FFT);

        if((i+1)%(nb_cryp/10)==0)
        {
            cout << "..." << (i+1)/(nb_cryp/10) << "0%" << flush;
        }
    }

    t2 = clock();
    diff = ((double)t2 - (double)t1)/1000000.0l;
    cout << "\n\nIt took " << diff << " seconds to do " << nb_cryp << " encryptions and decryptions." << endl;
    cout << "That's " << (diff/nb_cryp)*1000 << " milliseconds per encryption+decryption." << endl;
    cout << "That's " << (diff/nb_cryp)*1000*1024/N0 << " milliseconds per encryption+decryption per Kilobit." << endl << endl;
}

void Extract_Test(const unsigned int nb_extr, MSK_Data * MSKD)
{
    unsigned int i, rep;
    vec_ZZ id;
    ZZX SK_id[2];

    rep = 0;

    cout << "0%" << flush;
    for(i=0; i<nb_extr; i++)
    {
        id = RandomVector();

        IBE_Extract(SK_id, id, MSKD);
        rep += IBE_Verify_Key(SK_id, id, MSKD);
        if((i+1)%(nb_extr/10)==0)
        {
            cout << "..." << (i+1)/(nb_extr/10) << "0%" << flush;
        }
    }

    cout << endl;
    if(rep == 0)
    {    cout << endl << nb_extr << " extractions successfully performed!" << endl << endl;    }
    else
    {    cout << endl << rep << " out of " << nb_extr << " extractions failed miserabily!" << endl << endl;    }
}

void Encrypt_Test(const unsigned int nb_voter, const unsigned int nb_shares, MPK_Data * MPKD, MSK_Data * MSKD)
{
    unsigned int i, j, k, y, rep;
    vec_ZZ id;
    ZZX SK_id[2], m, w;
    CC_t SKid_FFT[N0];
    long int id0[N0], Ciphertext[2][N0];
    long int message[N0], decrypted[N0];
    long int int_value;


    cout << "0%" << flush << endl;

    std::vector<int> values;  // Vector to store the integers
    std::ifstream infile("../commitment.txt");  // Open the file
    double temp_value;  // Temporary variable to hold the floating-point value
    int line_count = 0;  // Add this counter

    if (infile.is_open()) {
        double value;
        while (infile >> temp_value) {  // Read each integer
            int_value = static_cast<int>(temp_value);  // Convert to int
            values.push_back(int_value);
            line_count++;  // Increment the counter
        }
        infile.close();  // Close the file
    } else {
        std::cerr << "Unable to open file" << std::endl;
    }

    std::cout << "Total lines read: " << line_count << std::endl;  // Print the count
    // Output the values to verify
    // for (const int & val : values) {
    //     std::cout << val << std::endl;
    // }

    // Convert each integer to a binary array of size N0
    std::vector<std::bitset<N0>> binary_arrays;  // Vector to store the binary arrays

    for (const int& val : values) {
        std::bitset<N0> binary_array(val);  // Convert integer to binary and pad with zeros
        binary_arrays.push_back(binary_array);
    }

    // Output the binary arrays
    // for (const auto& binary_array : binary_arrays) {
    //     std::cout << binary_array << std::endl;  // Print each binary array
    // }

    for(i=0; i<nb_voter; i++) {
        id = RandomVector();
        IBE_Extract(SK_id, id, MSKD);
        IBE_Verify_Key(SK_id, id, MSKD);
        ZZXToFFT(SKid_FFT, SK_id[1]);

        rep = 0;

        for(y=0; y<N0; y++) {
            id0[y] = conv<long int>(id[y]);
        }

        cout << "Voter: " << i << endl;

        for(k=0; k<nb_shares; k++) {
            cout << "share " << k << ": ";

            // for(j=0; j<N0; j++) {
            //     message[j] = (rand()%2);
            //     cout << message[j];
            // } cout << endl;

            for(j=0; j<N0; j++) {
                message[j] = binary_arrays[0][j];
                cout << message[j];
            } cout << endl;

            // Convert message with binary shares
            // int index = N0 - 1; // Start from the end of the array
            // while (share > 0 && index >= 0) {
            //     message_2[index] = share % 2;
            //     share /= 2;
            //     --index;
            // }

            IBE_Encrypt(Ciphertext, message, id0, MPKD);
            IBE_Decrypt(decrypted, Ciphertext, SKid_FFT);

            // for(int k=0; k<2; k++) {
            //     for(int g=0; g<N0; g++) {
            //         cout << Ciphertext[k][g];
            //     } cout << endl;
            // }
            
            for(j=0; j<N0; j++)
            {
                if(message[j] != decrypted[j])
                {
                    cout << "ERROR : Dec(Enc(m)) != m " << endl;
                    rep++;
                    break;
                }
            }
            // cout << "message_2: " << endl;
            // for(j=0; j<N0; j++)
            // {
            //     cout << message_2[j] << " ";
            // }
            // cout << endl;

            // cout << "message: " << endl;
            // for(j=0; j<N0; j++)
            // {
            //     cout << message[j] << " ";
            // }
            // cout << endl;
            
            // cout << "decrypted: " << endl;
            // for(j=0; j<N0; j++)
            // {
            //     cout << decrypted[j] << " ";
            // }
            // cout << endl;

            if((i+1)%(nb_voter/10)==0)
            {
                cout << "..." << (i+1)/(nb_voter/10) << "0%" << flush;
            }
            cout << endl;
            if(rep == 0)
            {    cout << endl << nb_shares << " encryptions+decryptions successfully performed!" << endl << endl;    }
            else
            {    cout << endl << rep << " out of " << nb_voter << " encryptions+decryptions failed miserabily!" << endl << endl;    }
        }
    }    
}

void TrusteeGen(const unsigned int nb_cryp, MPK_Data * MPKD, MSK_Data * MSKD) {
    unsigned int i;
    vec_ZZ id;
    ZZX SK_id[2], m;
    CC_t SKid_FFT[N0];
    long int id0[N0], Ciphertext[2][N0];
    long int message[N0], decrypted[N0];

    std::ofstream mpk_file("mpk_data.dat", std::ios::binary);
    // mpk_file << MPKD->MPK;
    mpk_file << MPKD;
    mpk_file.close();
    // cout << "mpk_data.dat saved" << endl;

    std::ofstream msk_file("msk_data.dat", std::ios::binary);
    // for (int i = 0; i < 4; i++) {
    //     msk_file << MSKD->MSK[i];
    // }
    msk_file << MSKD;
    msk_file.close();
    // cout << "msk_data.dat saved" << endl;

    // id = RandomVector();
    // IBE_Extract(SK_id, id, MSKD);
    // IBE_Verify_Key(SK_id, id, MSKD);
    // ZZXToFFT(SKid_FFT, SK_id[1]);

    // // Save id vector to a .dat file
    // std::ofstream idFile("id_data.dat", std::ios::binary);
    // if (idFile.is_open()) {
    //     for (i = 0; i < id.length(); i++) {
    //         long int temp = conv<long int>(id[i]); // Convert ZZ to long int
    //         idFile.write(reinterpret_cast<const char*>(&temp), sizeof(temp));
    //     }
    //     idFile.close();
    //     // cout << "id_data.dat saved" << endl;
    //     for (i = 0; i < id.length(); i++) {
    //         cout << id[i];
    //     } cout << endl;
    // }

    // // Save SKid_FFT to a .dat file
    // std::ofstream skFile("SKid_FFT_data.dat", std::ios::binary);
    // if (skFile.is_open())
    // {
    //     skFile.write(reinterpret_cast<const char*>(SKid_FFT), sizeof(SKid_FFT));
    //     skFile.close();
    //     // cout << "SKid_FFT_data.dat saved" << endl;
    // }
}

void VoterEncrypt(string msg) {
    ZZX MSK[4];
    ZZ_pX phiq, MPK;
    unsigned int i, j;
    float diff;
    MSK_Data * MSKD = new MSK_Data;
    MPK_Data * MPKD = new MPK_Data;
    
    vec_ZZ id;
    CC_t SKid_FFT[N0];
    long int id0[N0], Ciphertext[2][N0];
    long int message[N0], decrypted[N0];

    // cout << msg << endl;

    // Load the keys before performing encryption
    std::ifstream mpk_file("mpk_data.dat", std::ios::binary);
    if (!mpk_file.is_open()) {
        std::cerr << "Failed to open mpk_data.dat" << std::endl;
        return;  // Exit the function if the file can't be opened
    }
    mpk_file >> MPKD->MPK;
    // mpk_file >> MPKD;
    // Load other MPK_Data members if needed
    mpk_file.close();
    // cout << "mpk_data.dat loaded" << endl;

    std::ifstream msk_file("msk_data.dat", std::ios::binary);
    for (int i = 0; i < 4; i++) {
        msk_file >> MSKD->MSK[i];
    }
    // msk_file >> MSKD;
    // Load other MSK_Data members if needed
    msk_file.close();
    // cout << "msk_data.dat loaded" << endl;

    // Load id vector from the .dat file
    std::ifstream idFile("id_data.dat", std::ios::binary);
    if (idFile.is_open()) {
        long int temp;
        while (idFile.read(reinterpret_cast<char*>(&temp), sizeof(temp))) {
            id.append(to_ZZ(temp)); // Convert long int back to ZZ and append to id vector
        }
        idFile.close();
        // cout << "id_data.dat loaded" << endl;
    }
    else {
        cout << "Failed to open id_data.dat" << std::endl;
    }

    // Load SKid_FFT array from the .dat file
    std::ifstream skFile("SKid_FFT_data.dat", std::ios::binary);
    if (skFile.is_open()) {
        skFile.read(reinterpret_cast<char*>(SKid_FFT), sizeof(SKid_FFT));
        skFile.close();
        // cout << "SKid_FFT_data.dat loaded" << endl;
    }
    else {
        std::cerr << "Failed to open SKid_FFT_data.dat" << std::endl;
    }

    // // (Optional) Print or use loaded data to verify
    // std::cout << "id vector:" << std::endl;
    // for (i = 0; i < id.length(); i++)
    // {
    //     std::cout << id[i] << " ";
    // }
    // std::cout << std::endl;

    // std::cout << "SKid_FFT array:" << std::endl;
    // for (i = 0; i < N0; i++)
    // {
    //     std::cout << SKid_FFT[i] << " ";
    // }
    // std::cout << std::endl;

    for(i=0; i<N0; i++)
    {
        id0[i] = conv<long int>(id[i]);
    }
    // cout << "0%" << flush;
    // cout << endl;

    // cout << msg << " ";
    // Step 1: Copy the binary string to the message array
    int msg_length = msg.length();
    for (j = 0; j < msg_length; j++) {
        message[j] = msg[j] - '0'; // Convert char '0'/'1' to int 0/1 and store in the array
    }
    // Step 2: Pad the remaining part of the message array with 0 (or another padding value)
    for (int j = msg_length; j < N0; j++) {
        message[j] = 0; // Padding with 0
    }

    for(j=0; j<N0; j++) {
        message[j] = (rand()%2);
    }

    // Optional: Print the message array to verify
    cout << "msg: ";
    for (int j = 0; j < N0; j++) {
        cout << message[j];
    } cout << endl;

    IBE_Encrypt(Ciphertext, message, id0, MPKD);
    IBE_Decrypt(decrypted, Ciphertext, SKid_FFT);

    cout << "dec: ";
    for(j=0; j<N0; j++) {
        cout << decrypted[j];
    } cout << endl;
    
    for(j=0; j<N0; j++)
    {
        if(message[j] != decrypted[j])
        {
            cout << "ERROR : Dec(Enc(m)) != m " << endl;
            // rep++;
            break;
        }
    }
}
