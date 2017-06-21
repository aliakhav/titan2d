/*******************************************************************
 * Copyright (C) 2003 University at Buffalo
 *
 * This software can be redistributed free of charge.  See COPYING
 * file in the top distribution directory for more details.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Author:
 * Description:
 *
 *******************************************************************
 */

#ifndef TITAN2D_SIMULATION_H
#define TITAN2D_SIMULATION_H

#include <string>
#include <vector>
#include <algorithm>

#include "hd5calls.h"

#include "../gisapi/GisApi.h"

#include "integrators.h"

#include "properties.h"

class ElementsHashTable;
class NodeHashTable;


/**
 * cxxTitanSimulation
 */
class cxxTitanSimulation
{
public:
    cxxTitanSimulation();
    ~cxxTitanSimulation();

    void set_short_speed(bool short_speed);
    void set_geoflow_tiny(double _geoflow_tiny);
    //!>Process input and initiate dependencies, replacing Read_data
    void process_input(bool start_from_restart);

    void run(bool start_from_restart);
    void input_summary();

    void save_restart_file();

    void load_restart(const std::string restartFilename, const std::string new_gis_main);

    void save_vizoutput_file(const int mode);

    int myid;
    int numprocs;

    TiScale scale_;


    //! adapt
    int adapt;

    //!use a GIS material map
    bool use_gis_matmap;
    /**
     * vizoutput is used to determine which viz output to use
     * nonzero 1st bit of viz_flag means output tecplotxxxx.tec
     * nonzero 2nd bit of viz_flag means output mshplotxxxx.tec (debug purposes)
     * nonzero 3rd bit of viz_flag means output Paraview/XDMF format
     * nonzero 4th bit of viz_flag means output grass_sites files
     */
    int vizoutput;

    std::string vizoutput_prefix;
    std::string restart_prefix;
    bool restart_keep_all;
    bool restart_keep_redundant_data;
    bool restart_enabled;

    bool overwrite_output;

    //!Integrator
    Integrator *integrator;


    //!>Piles
    PileProps* pileprops;
    //PileProps pileprops_single_phase;
    //PilePropsTwoPhases pileprops_two_phases;

    //!>Flux sources
    FluxProps fluxprops;


    //!>Discharge planes
    DischargePlanes discharge_planes;
    //std::vector<cxxTitanDischargePlane> discharge_planes;

    //!>MatProps
    MatProps* matprops;
    //MatProps matprops_single_phase;
    //MatPropsTwoPhases matprops_two_phases;

    StatProps *statprops;
    TimeProps timeprops;
    MapNames mapnames;
    OutLine outline;

    NodeHashTable* NodeTable;
    ElementsHashTable* ElemTable;

    PileProps* get_pileprops(){return pileprops;}
    void set_pileprops(PileProps* m_pileprops);


    MatProps* get_matprops(){return matprops;}
    void set_matprops(MatProps* m_matprops);

    Integrator* get_integrator(){return integrator;}
    void set_integrator(Integrator* m_integrator);



    FluxProps* get_fluxprops(){return &fluxprops;}
    DischargePlanes* get_discharge_planes(){return &discharge_planes;}
    StatProps* get_statprops(){return statprops;}
    TimeProps* get_timeprops(){return &timeprops;}
    MapNames* get_mapnames(){return &mapnames;}
    OutLine* get_outline(){return &outline;}

    NodeHashTable* get_HT_Node(){return NodeTable;}
    ElementsHashTable* get_HT_Elem(){return ElemTable;}

    void set_element_type(const ElementType m_elementType);
    const ElementType& get_element_type() const{return elementType;}

protected:
    /** this function intializes the piles, by commenting/uncommenting define statements you can switch from
     * parabaloid to elliptical cylinder shaped piles, or even a hard coded pileshapes written to match particular
     * experiments.  Adaptive remeshing and pile reinitialization helps detect small piles and refine around pile
     * edges to obtain a more accurate initial solution and speed up the first few timesteps before adaptive
     * refinement and unrefinement would otherwise occur.
     */
    void init_piles();
    /**
     * write some of the content to h5group
     */
    void h5write(H5::CommonFG *parent);
    void h5read(const H5::CommonFG *parent);

    //if Quad9==true write Quad9 otherwise Quad4
    void xmdfWrite(const char *xmdf_filename,const char *hf5_filename,bool Quad9);

    ElementType elementType;
    std::string prev_hf5_filename;
    std::string prev_xmdf_snapshot;

    std::vector<double> N_elements;
    std::vector<double> Time_el;
    std::vector<double> Area;
    std::vector<double> Finertx;
    std::vector<double> Finerty;
    std::vector<double> Fconvx;
    std::vector<double> Fconvy;
    std::vector<double> Fgravx;
    std::vector<double> Fgravy;
    std::vector<double> Fbx;
    std::vector<double> Fby;
    std::vector<double> Fix;
    std::vector<double> Fiy;

    std::vector<double> Pinert;
    std::vector<double> Pconv;
    std::vector<double> Pgrav;
    std::vector<double> Pbed;
    std::vector<double> Pint;

    std::vector<double> FVinert;
    std::vector<double> FVconv;
    std::vector<double> FVgrav;
    std::vector<double> FVbed;
    std::vector<double> FVint;

    std::vector<double> XX;
    std::vector<double> YY;
    std::vector<double> X_el;
    std::vector<double> Y_el;

    double T_Finertx;
    double T_Finerty;
    double T_Fconvx;
    double T_Fconvy;
    double T_Fgravx;
    double T_Fgravy;
    double T_Fbx;
    double T_Fby;
    double T_Fix;
    double T_Fiy;

    double T_Pinert;
    double T_Pconv;
    double T_Pgrav;
    double T_Pbed;
    double T_Pint;

};

#endif
