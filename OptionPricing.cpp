#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <random>
#include <array>
#include <map>
using namespace std;


// Generic functions
double Normal_cdf(double x) // N(x)
{
    return erfc(-x / sqrt(2)) / 2;
}

double Normal_pdf(double x) // N'(x)
{
    return exp(-pow(x,2)*0.5)/sqrt(2*M_PI);
}

float Round_x(float x, int places) // Round to x decimal places 
{
    return round(x*pow(10,places))/pow(10,places);
}

// Financial instruments
class FixedIncome {
    private:
        float r;

    public:
    // WIP - Triangle arbitrage for forwards and spot currency exchange rates
    // WIP - FRA
    // WIP - IRS

};

class Option {
    private:
        float St, t, T, K, sigma, r;

    public:
        void Set_inputs(float SharePrice, float timet, float timeT, float Strike, float Volatility, float Rate) {
            St = SharePrice;
            t = timet;
            T = timeT;
            K = Strike;
            sigma = Volatility;
            r = Rate;
        }

        void Describe() {
            cout << "Basic description: \nShare price: " << St << " Strike Price: " << K << " From t= " << t << " to T= " << T;
        }

        std::map<std::string, float> BSM_value_call() {
            float call, d1, d2;
            d1 = (1/sigma * sqrt(T - t))*(log(St/K)+(r+pow(sigma,2)/2)*(T-t));
            d2 = d1 - sigma * sqrt(T - t);

            // Implement normal distribution for d1 and d2
            call = Normal_cdf(d1)*St-Normal_cdf(d2)*K*(::exp(-r*(T-t)));

            // Greeks
            float delta, gamma, vega, theta, rho;
            delta = Normal_cdf(d1);
            gamma = Normal_pdf(d1)/(St*sigma*sqrt(T-t));
            vega = St*Normal_pdf(d1)*sqrt(T-t);
            theta = -St*Normal_pdf(d1)*sigma/(2*sqrt(T-t))-r*K*exp(-r*(T-t)*Normal_cdf(d2));
            rho = K*(T-t)*exp(-r*(T-t))*Normal_cdf(d2);

            // Return list
            std::map<std::string, float> return_dict;
                return_dict["Call value"]=call;
                return_dict["Delta"]=delta;
                return_dict["Gamma"]=gamma;
                return_dict["Vega"]=vega;
                return_dict["Theta"]=theta;
                return_dict["Rho"]=rho;

            return return_dict;
        }

        std::map<std::string, float> BSM_value_put() {
            float put, d1, d2;
            d1 = (1/sigma * sqrt(T - t))*(log(St/K)+(r+pow(sigma,2)/2)*(T-t));
            d2 = d1 - sigma * sqrt(T - t);

            // Implement normal distribution for d1 and d2
            put = Normal_cdf(-d2)*K*(::exp(-r*(T-t)))-Normal_cdf(-d1)*St;

                        // Greeks
            float delta, gamma, vega, theta, rho;
            delta = Normal_cdf(d1)-1;
            gamma = Normal_pdf(d1)/(St*sigma*sqrt(T-t));
            vega = St*Normal_pdf(d1)*sqrt(T-t);
            theta = -St*Normal_pdf(d1)*sigma/(2*sqrt(T-t))+r*K*exp(-r*(T-t)*Normal_cdf(-d2));
            rho = -K*(T-t)*exp(-r*(T-t))*Normal_cdf(-d2);

            // Return list
            std::map<std::string, float> return_dict;
                return_dict["Put value"]=put;
                return_dict["Delta"]=delta;
                return_dict["Gamma"]=gamma;
                return_dict["Vega"]=vega;
                return_dict["Theta"]=theta;
                return_dict["Rho"]=rho;

            return return_dict;
        }
};


int main()
{
    // General variables
    string CURRENCY;
    float Rf;

    Rf = 0.05;
    CURRENCY = "$";

    // Creating a standard european call
    Option opt_eur;
    opt_eur.Set_inputs(45,0,0.25,43,0.2,Rf);
    opt_eur.Describe();

    std::map<string,float> call_bsm;
    call_bsm = opt_eur.BSM_value_call();

    cout << "\n\nCall (BSM value and Greeks): \n";
    map <string,float>::iterator i_c;
    for (i_c=call_bsm.begin();i_c!=call_bsm.end();i_c++) {
        cout << "\n" << (*i_c).first << " " << (*i_c).second;
    }

    std::map<string,float> put_bsm;
    put_bsm = opt_eur.BSM_value_put();

    cout << "\n\nPut (BSM value and Greeks): \n";
    map <string,float>::iterator i_p;
    for (i_p=put_bsm.begin();i_p!=put_bsm.end();i_p++) {
        cout << "\n" << (*i_p).first << " " << (*i_p).second;
    }

    return 0;
}