/******************************************************************************
 $Id: vegetation.cc,v 1.9 2004/09/23 12:41:36 schatz Exp $
 ******************************************************************************/

#include <math.h>

#include "globals.h"
#include "func.h"
#include "initsurf.h"
#include "vegetation.h"

//*****************************************************************************
//  class vegetation

vegetation::vegetation(const dunepar& p) : dunedata(p)
{
    m_veget.Create( duneglobals::nx(), duneglobals::ny(), duneglobals::dx() );

    /*Vegetation parameters*/
    m_Lveg = p.getdefault("veget.xmin", 0.0);
    m_xmin0 = m_Lveg/duneglobals::dx();
    m_xmin = m_xmin0;
    m_zmin = p.getdefault("veget.zmin", 0.0);

    // New edits ODV This is not time anymore! but a length in m (kindo of ratio of vegetation volume to cover area: effective vegetation height)
    m_Tveg = p.getdefault("veget.Tveg", 1.0);       // (m)
    cout << "grass constructor: Tveg = " << m_Tveg << endl;

    // initial value
    m_veget_init0 = p.getdefault("veget.0", 0.0); // fix to 10%

    // erosion/acc
    m_sens = p.getdefault("veget.sensitivity", 1.0);       // species sensitivity for erosion/acc. rate
    m_Hveg = p.getdefault("veget.Hveg", 1.0/*0.2*/);       // maximum plant height (~ 1m)
    m_sens /= m_Hveg;

    // time conversion
    m_wind_factor =  duneglobals::timefrac();

    // different species
    m_spec1 = p.getdefault("veget.spec.1", 1);
    m_spec2 = p.getdefault("veget.spec.2", 0);

    // lateral growth
    m_lateral = p.getdefault("veget.lateralgrowth", false);
    m_Vlateral_factor = p.getdefault("veget.Vlateral.factor", 1.0);

    // extra
    m_rho_max = p.getdefault("veget.rho.max", 1.0);
    m_rho_min = p.getdefault("veget.rho.min", 0.0);

    m_survive = p.getdefault("veget.survive", false);

}

/*! Initialize vegetation */
void vegetation::init(const dunepar& par)
{
    arrayinit *init_veget;

    if( par.exists("veget.Init-Surf") )
        init_veget= arrayinit::create(par, "veget.");
    else
        init_veget= new CInitSurfPlain(0.0);

    init_veget->init_2d_vec( m_veget );

    delete init_veget;

//	// INIT VEGET
//	for( int y= 0; y< duneglobals::ny(); ++y ){
//         for( int x= 60; x< duneglobals::nx(); ++x ){
//         m_veget(x,y)[0]= 0.1;
//        }
//         }


}

/*! return cover fraction*/
void vegetation::getcover(TFktScal& rho_veget)
{
    // NORMALIZATION
    for( int x= 0; x< duneglobals::nx(); ++x )
        for( int y= 0; y< duneglobals::ny(); ++y ){
            // added density:
            rho_veget(x,y) = m_veget(x,y)[0] + 0*m_veget(x,y)[1];
            if (rho_veget(x,y) > m_rho_max) {
                rho_veget(x,y) = m_rho_max;
            }
        }

}

/*!  Computes the evolution of the vegetation.  */
int vegetation::evol(TFktScal& rho_veget, const double time, const double timestep, const double shoreline, const TFktScal& h, const TFktScal& dh_dt)
{
    int veget_X0 = 0, veget_X1 = 0;

    // SPEC 1
    veget_X0 = evolspec(time, timestep, shoreline, h, dh_dt, 0);
    // SPEC 2
    //veget_X1 = evolspec(time, timestep, shoreline, h, dh_dt, 1);

    getcover(rho_veget);

    // set overwash to zero
    //overwash.SetAll(0.0);
    return veget_X0;
}

/*!  Computes the evolution of each species.  */
int vegetation::evolspec(const double time, const double timestep, const double shoreline, const TFktScal& h, const TFktScal& dh_dt, int species)
{
    // VEGET GRAD (for lateral propagation)
    TFktVec grad_veget;
    grad_veget.Create(duneglobals::nx(), duneglobals::ny(), duneglobals::dx(),vec2(0.0,0.0));

    TFktScal veget_aux;
    veget_aux.Create(duneglobals::nx(), duneglobals::ny(), duneglobals::dx());
    for( int y= 0; y< duneglobals::ny(); ++y ){
        for( int x= 0; x< duneglobals::nx(); ++x ){
            veget_aux(x,y) = m_veget(x,y)[species];
        }
    }
    grad_veget.GradMin(veget_aux);

    // CALCULATION OF VEG LIMIT RELAXATION

    // repose angle for normalization
    double m_angle_ref= 15.; // 15.
    double m_tan_ref= tan(m_angle_ref * M_PI / 180.0);

    // calculate gradient
    TFktVec grad_h;
    grad_h.Create(duneglobals::nx(), duneglobals::ny(), duneglobals::dx());
    grad_h.GradMid(h);

    // General variables
    double growthrate, proprate, reprod_rate, mortal_rate;
    int vegetpoints = 0;

    double m_shore_HMWL = duneglobals::HMWL();


    // Species growth parameters:
    // GENERIC (1)
    double H_V_inverse = 1./m_Tveg;
    double Veg_BETA = m_Vlateral_factor;
    //double dhdt_c = 1e-5; //originally this was 1e-7...
    double dhdt_c = 1;


    for(int y = 0; y< duneglobals::ny(); ++y ){
        for(int x = 0; x< duneglobals::nx()-1; ++x ){

            // AUXILIAR
            double dhdt = dh_dt(x,y); // erosion rate
            double erosion = (dhdt < 0 ? 1 : 0);

            // absolute value gradient
            double abs_grad_veget = sqrt(grad_veget(x,y)[0]*grad_veget(x,y)[0]+grad_veget(x,y)[1]*grad_veget(x,y)[1]);

            // define an alternative cover density encoding competition
            double rho_competition = m_veget(x,y)[0] + 0*m_veget(x,y)[1];
            rho_competition = (rho_competition > m_rho_max ? m_rho_max : rho_competition);
            // species limiting factors
            double shorefactor = (x < shoreline + m_xmin ? 0 : 1);  // 0: maximum effect; 1: no effect

            // needed for lateral propagation
            double abs_grad_h = sqrt(grad_h(x,y)[0]*grad_h(x,y)[0]+grad_h(x,y)[1]*grad_h(x,y)[1]);
            double dhdxfactor = (1-abs_grad_h/m_tan_ref/*/0.4*/);

            //  VERTICAL GROWTH

            double G0 = (species == 0 ? H_V_inverse : 0) * (dhdt/dhdt_c) * (1. - erosion) * shorefactor;

            reprod_rate =  G0 * (m_veget(x,y)[species]) * (1 - rho_competition);

            // LATERAL PROPAGATION
            double CC = Veg_BETA * (dhdt/dhdt_c) *  (1. - erosion) * (dhdxfactor > 0 ? 1 : 0) * (m_veget(x,y)[species] < 1 ? 1 : 0);

            proprate =  CC * abs_grad_veget;

            //DEATH

            mortal_rate =  (1/m_sens) * (m_veget(x,y)[species]) * fabs(dhdt) * erosion;

            //growthrate = proprate + reprod_rate - mortal_rate;
            growthrate = proprate + reprod_rate;

            // evolution of cover fraction
            m_veget(x,y)[species]+= timestep * growthrate; // * m_wind_factor;

            // limiting conditions
            if(m_veget(x,y)[species] > 1 ){
                m_veget(x,y)[species] = 1;
            }

            if(m_veget(x,y)[species] < 0 || shorefactor == 0){
                m_veget(x,y)[species] = 0;
            }

            // initial condition
            if (m_veget(x,y)[species] == 0 && shorefactor * dhdxfactor > 0 /*&& rand() % N == 1*/){ // 50000 ref. case for X=4
                m_veget(x,y)[species] = m_veget_init0;
            }


        }
    }

    return m_xmin;
}

/* Shift back*/
void vegetation::shiftback(const int plusminus)
{
    m_veget.ShiftOne(plusminus);
}
/*!  Saves the arrays m_u and m_rho.  */
void vegetation::save_arrays(){
    save_2d_vecarray( "veget", m_veget );
}
