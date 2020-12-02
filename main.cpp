#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

int main() {
    using namespace std;

    int ndx = 100; // pocet uzlu na x
    double Termi = 1.; // doba reseni
    double Cour = 0.1; // courantovo c. ulohy (1 = mez stability pro linearni ulohu, zde je stabilita cca na 0,1)
    double h_prehrady = 2.; // pocatecni vyska hladiny prehrady
    double h_reky = 1.; // pocatecni vyska hladiny reky za prehradou
    double l_prehrady = 10.; // delka prehrady
    double l_reky = 10.; // delka reky za prehradou
    double u_prehrady = 0.; // horizontalni rychlost prehrady
    double u_reky = 1.; // horizontalni rychlost reky

    const double g = 9.81; // tihove zrychleni

//______________________________SOLVER______________________________________
// urceni casoveho kroku s ohledem na stabilitu
    double dx = (l_prehrady + l_reky)/ndx;
    double dt = Cour * dx;
// pocet casovych kroku
   int ndt = round(Termi / dt);
// pole reseni konzerv. promennych
    long double fi[ndt][ndx]={0};
    long double m[ndt][ndx]={0};
// pocatecni podminky
    double ndx_prehrady = round((l_prehrady)/(l_prehrady + l_reky) * ndx);
    double ndx_reky = ndx - ndx_prehrady;
// zapis pocatecnich podminek
    for (int i=0; i < ndx; i++) {
        if (i<ndx_prehrady) {
            fi[0][i] = g * h_prehrady;
            m[0][i] = g * h_prehrady * u_prehrady;
        }
        else {
            fi[0][i] = g * h_reky;
            m[0][i] = g * h_reky * u_reky;
        }
    }

// RESENI
   for (int n=0; n < ndt-1; n++) {     //cyklus pres cas
       fi[n + 1][0] = fi[n][0];
       fi[n + 1][ndx-1] = fi[n][ndx-1];
       m[n + 1][0] = m[n][0];
       m[n + 1][ndx-1] = m[n][ndx-1];

       for (int k = 1; k < ndx - 1; k++) {
        fi[n + 1][k] = (1./2)*(fi[n][k + 1] + fi[n][k - 1]) - (Cour/2.)*(m[n][k + 1] - m[n][k - 1]);
        long double Fp = (m[n][k + 1]*m[n][k + 1]/fi[n][k + 1]) + (fi[n][k + 1]*fi[n][k + 1]/2);
        long double Fl = (m[n][k - 1]*m[n][k - 1]/fi[n][k - 1]) + (fi[n][k - 1]*fi[n][k - 1]/2);
        m[n + 1][k] = (1./2)*(m[n][k + 1] + m[n][k - 1]) - (Cour/2)*(Fp - Fl);
       }
   }

    long double H[ndt][ndx] = {0};
    long double U[ndt][ndx] = {0};
    for(int i = 0; i < 50; i++) {
        for (int j = 0; j < 100; j++) {
            H[i][j] = fi[i][j]/g;
            U[i][j] = m[i][j]/H[i][j]/g;
        }
   }

    ofstream soubor1("vyska_hladiny.txt");
    for (int i = 0; i < 50; i++)
    {
        for (int j = 0; j < 100; j++)
            {
                soubor1 << std::right << std::setw(8) << H[i][j] ;
            }
        soubor1 << std::endl;
    }
    soubor1.close();

    ofstream soubor2("horizontalni_rychlost.txt");
    for (int i = 0; i < 50; i++)
    {
        for (int j = 0; j < 100; j++)
        {
            soubor2 << std::right << std::setw(13) << U[i][j] ;
        }
        soubor2 << std::endl;
    }
    soubor2.close();









// kontrolni vypis matic
    for (int i = 0; i < 50; i++)
    {
        for (int j = 0; j < 100; j++)
        {
            {
                std::cout << std::right << std::setw(8) << fi[i][j];
            }
        }

        std::cout << std::endl;
    }
    std::cout << std::endl;
    for (int i = 0; i < 50; i++)
    {
        for (int j = 0; j < 100; j++)
        {
            std::cout << std::right << std::setw(13) << m[i][j];
        }

        std::cout << std::endl;
    }
    std::cout << std::endl;
    for (int i = 0; i < 50; i++)
    {
        for (int j = 0; j < 100; j++)
        {
            {
                std::cout << std::right << std::setw(8) << H[i][j];
            }
        }

        std::cout << std::endl;
    }
    std::cout << std::endl;
    for (int i = 0; i < 50; i++)
    {
        for (int j = 0; j < 100; j++)
        {
            {
                std::cout << std::right << std::setw(8) << U[i][j];
            }
        }

        std::cout << std::endl;
    }
    return 0;
}
