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

#ifdef HAVE_CONFIG_H
# include <titan_config.h>
#endif

#include <sstream>
#include <iostream>
#include <fstream>

#include "../header/hpfem.h"
#include "../header/hadapt.h"
#include "../header/stats.hpp"

#if HAVE_HDF5_H
#include "../header/hd5calls.h"
#endif

//#define DEBUG
//#define LOAD_BAL_DEBUG  //turns on a whole mess of mpi barriers so it makes run time more sensitive to load imbalances i.e. more sensitive to the load balance weights, it just makes it easier to adjust the constants.
//#define PERFTEST
#define TARGETPROC  -1
#define TARGETPROCA -1

int REFINE_LEVEL = 3;


#include "../header/elenode.hpp"
#include "../header/titan2d_utils.h"

#include "../header/titan_simulation.h"

TitanTimings titanTimings;
TitanTimings titanTimingsAlongSimulation;
TitanProfiling titanProfiling;
int threads_number;


#include <stdio.h>
#include <string.h>
#include "../header/ticore/omp_mpi.hpp"
#include "../header/titan_simulation.h"
#include <math.h>
#include "../header/constant.h"
#include "../header/properties.h"

int NUM_STATE_VARS;
bool SHORTSPEED;
double GEOFLOW_TINY;

cxxTitanSimulation::cxxTitanSimulation() :
        integrator(nullptr),
        matprops(nullptr),
        pileprops(nullptr),
        NodeTable(nullptr),
        ElemTable(nullptr)
{
    elementType=ElementType::UnknownElementType;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    IF_MPI(MPI_Barrier (MPI_COMM_WORLD));



    SHORTSPEED=false;
    GEOFLOW_TINY=0.0001;

    set_element_type(ElementType::SinglePhase);

    NodeTable=new NodeHashTable();
    ElemTable=new ElementsHashTable(NodeTable);

    pileprops=nullptr;

    set_element_type(ElementType::SinglePhase);

    statprops=new StatProps(ElemTable,NodeTable);

    vizoutput_prefix="";
    restart_prefix="";

    prev_hf5_filename="";
    prev_xmdf_snapshot="";

    restart_keep_all=false;
    restart_keep_redundant_data=false;
    restart_enabled=true;

    outline.scale=&scale_;

    IF_MPI(MPI_Barrier (MPI_COMM_WORLD));
}
cxxTitanSimulation::~cxxTitanSimulation()
{
    FREE_VAR_IF_NOT_NULLPTR(integrator);
    FREE_VAR_IF_NOT_NULLPTR(ElemTable);
    FREE_VAR_IF_NOT_NULLPTR(NodeTable);
    FREE_VAR_IF_NOT_NULLPTR(matprops);
    FREE_VAR_IF_NOT_NULLPTR(pileprops);
    FREE_VAR_IF_NOT_NULLPTR(statprops);
}
void cxxTitanSimulation::set_element_type(const ElementType m_elementType)
{
    elementType=m_elementType;

    if(elementType==ElementType::SinglePhase)
    {
        NUM_STATE_VARS = 3;
    }
    else if(elementType==ElementType::TwoPhases)
    {
        NUM_STATE_VARS = 6;
    }
    else
    {
        printf("Unknown type of element!\n");
        assert(0);
    }
    get_outline()->elementType=elementType;
    if(ElemTable!=nullptr)
        ElemTable->set_element_type(elementType);
    if(NodeTable!=nullptr)
        NodeTable->set_element_type(elementType);

}
void cxxTitanSimulation::set_short_speed(bool short_speed)
{
    SHORTSPEED=short_speed;
}
void cxxTitanSimulation::set_geoflow_tiny(double _geoflow_tiny)
{
    GEOFLOW_TINY=_geoflow_tiny;
}
void cxxTitanSimulation::process_input(bool start_from_restart)
{

    int i;
    int isrc;
    double doubleswap;

    MatProps *matprops_ptr=get_matprops();

    PileProps* pileprops_ptr=get_pileprops();
    FluxProps* fluxprops_ptr=get_fluxprops();

    StatProps* statprops_ptr=get_statprops();
    TimeProps* timeprops_ptr=get_timeprops();
    MapNames* mapnames_ptr=get_mapnames();
    OutLine* outline_ptr=get_outline();

    /*************************************************************************/

    int no_of_sources = fluxprops.no_of_sources;
    if(fluxprops.no_of_sources+pileprops_ptr->numpiles==0)
    {
        printf("ERROR: No material source was defined");
        assert(0);
    }

#ifdef STATISTICS_IS_NOT_TRANSFERED
    /*************************************************************************/
    //over ride regular input with statistical sample run data
    FILE* fp;
    int lhsref, lhsid, ifstatbed = 0;
    double statbed, statint;  //friction angles
    if((fp = fopen("statin.bed", "r")) != NULL)
    {
        ifstatbed = 1;
        fscanf(fp, "%d%d%lf%lf", &(statprops_ptr->lhs.refnum), &(statprops_ptr->lhs.runid), &statbed, &statint);
        fclose(fp);
        statbed *= PI / 180.0;
        statint *= PI / 180.0;
    }
    else if((fp = fopen("statin.vol", "r")) != NULL)
    {
        double volumescale;
        fscanf(fp, "%d%d%lf%lf", &(statprops_ptr->lhs.refnum), &(statprops_ptr->lhs.runid), &volumescale);
        fclose(fp);
        volumescale = pow(volumescale, 1 / 3.0);

        for(i = 0; i < numpiles; i++)
        {
            pileprops_ptr->pileheight[i] *= volumescale;
            pileprops_ptr->majorrad[i] *= volumescale;
            pileprops_ptr->minorrad[i] *= volumescale;
        }
    }

    //PCQ
    int isample = -1, ifpcqvolbed = 0;
    double pcqbedfrict;
    fp = fopen("sample.number", "r");
    if(fp != NULL)
    {
        fscanf(fp, "%d", &isample);
        fclose(fp);
        statprops_ptr->lhs.runid = isample;
        double volume;
        int intswap;
        char samplefilename[256];
        sprintf(samplefilename, "dirfluxsample.%06d", isample);
        fp = fopen(samplefilename, "r");
        if(fp != NULL)
        {
            assert(fluxprops.no_of_sources);
            fluxprops.no_of_sources = 1;

            char stringswap[4096];
            fgets(stringswap, 4096, fp);  //Nsample
            fscanf(fp, "isample=%d\n", &intswap);
            assert(intswap == isample);
            fgets(stringswap, 4096, fp);  //weight
            fgets(stringswap, 4096, fp);  //Nrand
            fgets(stringswap, 4096, fp);  //Nvol
            fgets(stringswap, 4096, fp);  //Ndir
            fgets(stringswap, 4096, fp);  //randvar
            fscanf(fp, "volume=%lf\n", &volume);
            double vel = sqrt(fluxprops.xVel[0] * fluxprops.xVel[0] + fluxprops.yVel[0] * fluxprops.yVel[0]);
            double vel_angle;
            fscanf(fp, "direction=%lf\n", &vel_angle);
            vel_angle *= PI / 180.;
            fluxprops.xVel[0] = vel * cos(vel_angle);
            fluxprops.yVel[0] = vel * sin(vel_angle);
            fluxprops.start_time[0] = 0.0;
            fscanf(fp, "flux duration=%lf\n", &(fluxprops.end_time[0]));
            fscanf(fp, "init center flux=%lf\n", &(fluxprops.influx[0]));
            printf("Vol=%g [m^3]\n", volume);
        }
    }
#endif
    if(statprops_ptr->lhs.runid>=0)
        statprops_ptr->runid = statprops_ptr->lhs.runid;

    /*************************************************************************/
    double TIME_SCALE;
    double VELOCITY_SCALE;
    if(start_from_restart==false)
    {
        matprops_ptr->set_scale(pileprops_ptr,fluxprops_ptr);
        matprops_ptr->process_input();

        TIME_SCALE = matprops_ptr->get_TIME_SCALE();
        //non-dimensionalize the inputs
        VELOCITY_SCALE = matprops_ptr->get_VELOCITY_SCALE();

        integrator->scale();

        pileprops_ptr->scale(scale_.length,scale_.height,scale_.gravity);
        fluxprops.scale(scale_.length,scale_.height,scale_.gravity);

        double smallestpileradius=min(pileprops_ptr->get_smallest_pile_radius(),fluxprops.get_smallest_source_radius());

        matprops_ptr->smallest_axis = 2.0 * smallestpileradius;
    }
    else
    {
        TIME_SCALE = matprops_ptr->get_TIME_SCALE();
        VELOCITY_SCALE = matprops_ptr->get_VELOCITY_SCALE();
    }

    //read in material map
    if(matprops_ptr->material_count > 1)
    {
        char *gis_matmap = (char *) malloc((mapnames_ptr->gis_map.size() + 5) * sizeof(char));
        strcpy(gis_matmap, mapnames_ptr->gis_map.c_str());
        strcat(gis_matmap + strlen(mapnames_ptr->gis_map.c_str()), "_Mat");

        if(Initialize_Raster_data(mapnames_ptr->gis_main.c_str(), mapnames_ptr->gis_sub.c_str(), mapnames_ptr->gis_mapset.c_str(), gis_matmap))
        {
            printf("Problem with GIS Material on processor %d\n", myid);
            assert(0);
        }
        free(gis_matmap);

        int nummat;
        Get_raster_categories(&nummat);
        if(nummat != matprops_ptr->material_count)
        {
            printf("ERROR: Input script has %d materials but material map has %d, aborting, use titan_materialnames utility to check materials\n", matprops_ptr->material_count,
                   nummat);
            assert(0);
        }
        //shuffle materials to raster order, since it can be entered out of order
        std::vector<std::string> matnames;
        std::vector<double> bedfrict;
        std::vector<double> tanbedfrict;
        //something somewhere counting from 1, so add dummy values
        matnames.push_back("nothing");
        bedfrict.push_back(12.0);
        tanbedfrict.push_back(12.0);

        char materialname[256];
        for(int imat = 1; imat <= nummat; ++imat)
        {
            Get_raster_category_name(imat, materialname);
            int imat0;
            for(imat0=1;imat0<=nummat;++imat0)
            {
                if(matprops_ptr->matnames[imat0]==materialname)
                    break;
            }
            if(imat0>nummat){
                printf("ERROR: Material %s specified  in material map is not present in input script, use titan_materialnames utility to check materials\n", materialname);
                assert(0);
            }
            matnames.push_back(matprops_ptr->matnames[imat0]);
            bedfrict.push_back(matprops_ptr->bedfrict[imat0]);
            tanbedfrict.push_back(matprops_ptr->tanbedfrict[imat0]);
        }
        matprops_ptr->matnames=matnames;
        matprops_ptr->bedfrict=bedfrict;
        matprops_ptr->tanbedfrict=tanbedfrict;
    }
    /*************************************************************************/
    //time related info
    if(start_from_restart==false)
    {
        timeprops_ptr->scale(TIME_SCALE);
    }
    /*************************************************************************/
    if(start_from_restart==false)
    {
        if(use_gis_matmap)mapnames_ptr->extramaps=1;
    }
    /*************************************************************************/
    // read in GIS information
    i = Initialize_GIS_data(mapnames_ptr->gis_main.c_str(), mapnames_ptr->gis_sub.c_str(), mapnames_ptr->gis_mapset.c_str(), mapnames_ptr->gis_map.c_str(), mapnames_ptr->gis_format);
    if(i != 0)
    {
        printf("Problem with GIS on processor %d\n", myid);
        assert(0);
    }

    if(start_from_restart==false)
    {
        matprops_ptr->calc_Vslump(pileprops_ptr,&fluxprops);
        /*************************************************************************/
        //test point information
        statprops_ptr->scale(matprops_ptr);

        /*************************************************************************/
        discharge_planes.scale(scale_.length);

        //the discharge plane section ends here
        /*************************************************************************/
    }


    /* physically we don't know how to specify a changing internal friction
     angle, we don't know how the contents of an avalanche changes as it
     picks up new material.  We don't know how to keep track of the material
     it picks up.  So we say that the internal friction angle is a constant
     */
#ifdef STATISTICS_IS_NOT_TRANSFERED
    //replace frict.data values of bed friction with ones from a
    //statistics sample run file read in above
    if(ifstatbed)
    {
        matprops_ptr->intfrict = statint;
        for(int imat = 1; imat <= matprops_ptr->material_count; imat++)
            matprops_ptr->bedfrict[imat] = statbed;
    }

    if(ifpcqvolbed)
    { //printf("PCQ bedfrict yada\n");
        matprops_ptr->bedfrict[1] = pcqbedfrict * PI / 180.0;

        double doubleswap = sqrt(tan(matprops_ptr->bedfrict[1]));
        for(int imat = 2; imat <= matprops_ptr->material_count; imat++)
        {
            matprops_ptr->bedfrict[imat] = matprops_ptr->bedfrict[1];
        }
    }
#endif
    /*************************************************************************/
    //to read in outline parameters here when it has been added
    return;
}
void cxxTitanSimulation::set_matprops(MatProps* m_matprops)
{
    FREE_VAR_IF_NOT_NULLPTR(matprops);
    matprops=m_matprops;

    //reset matprops pointer in all dependencies
    if(integrator!=nullptr)
    {
        integrator->matprops_ptr=matprops;
    }
}
void cxxTitanSimulation::set_pileprops(PileProps* m_pileprops)
{
    FREE_VAR_IF_NOT_NULLPTR(pileprops);
    pileprops=m_pileprops;

    //reset matprops pointer in all dependencies
    if(integrator!=nullptr)
    {
        integrator->pileprops_ptr=pileprops;
    }
}
void cxxTitanSimulation::set_integrator(Integrator* m_integrator)
{
    FREE_VAR_IF_NOT_NULLPTR(integrator);
    integrator=m_integrator;
}
void cxxTitanSimulation::h5write(H5::CommonFG *parent)
{
    ///////////////////////////////////////////////////////////////////////////
    //write cxxTitanSimulation
    H5::Group groupTitanSimulation(parent->createGroup("TitanSimulation"));

    TiH5_writeIntAttribute(groupTitanSimulation, REFINE_LEVEL);
    TiH5_writeIntAttribute(groupTitanSimulation, NUM_STATE_VARS);
    TiH5_writeIntAttribute(groupTitanSimulation, SHORTSPEED);
    TiH5_writeDoubleAttribute(groupTitanSimulation, GEOFLOW_TINY);
    TiH5_writeIntAttribute(groupTitanSimulation, myid);
    TiH5_writeIntAttribute(groupTitanSimulation, numprocs);
    TiH5_writeBoolAttribute(groupTitanSimulation, use_gis_matmap);
    TiH5_writeIntAttribute(groupTitanSimulation, vizoutput);
    TiH5_writeStringAttribute(groupTitanSimulation, vizoutput_prefix);
    TiH5_writeStringAttribute(groupTitanSimulation, restart_prefix);
    TiH5_writeBoolAttribute(groupTitanSimulation, restart_keep_all);
    TiH5_writeBoolAttribute(groupTitanSimulation, restart_keep_redundant_data);
    TiH5_writeIntAttribute(groupTitanSimulation, adapt);

    int titan2d_version[3]={TITAN2D_VERSION_MAJOR,TITAN2D_VERSION_MINOR,TITAN2D_VERSION_REVISION};
    TiH5_writeIntArrayAttribute(groupTitanSimulation, titan2d_version, 3);
    TiH5_writeScalarDataTypeAttribute(groupTitanSimulation, elementType, datatypeElementType);

    ///////////////////////////////////////////////////////////////////////////
    NodeTable->h5write(parent);
    ElemTable->h5write(parent);

    timeprops.h5write(parent);
    scale_.h5write(parent);
    integrator->h5write(parent);
    pileprops->h5write(parent);
    fluxprops.h5write(parent);
    discharge_planes.h5write(parent);
    matprops->h5write(parent);
    statprops->h5write(parent);
    mapnames.h5write(parent);
    outline.h5write(parent);


}
void cxxTitanSimulation::h5read(const H5::CommonFG *parent)
{
    ///////////////////////////////////////////////////////////////////////////
    //read cxxTitanSimulation
    H5::Group groupTitanSimulation(parent->openGroup("TitanSimulation"));

    TiH5_readIntAttribute(groupTitanSimulation, REFINE_LEVEL);
    TiH5_readIntAttribute(groupTitanSimulation, NUM_STATE_VARS);
    TiH5_readIntAttribute(groupTitanSimulation, SHORTSPEED);
    TiH5_readDoubleAttribute(groupTitanSimulation, GEOFLOW_TINY);
    TiH5_readIntAttribute(groupTitanSimulation, myid);
    TiH5_readIntAttribute(groupTitanSimulation, numprocs);
    TiH5_readBoolAttribute(groupTitanSimulation, use_gis_matmap);
    TiH5_readIntAttribute(groupTitanSimulation, vizoutput);
    TiH5_readStringAttribute(groupTitanSimulation, vizoutput_prefix);
    TiH5_readStringAttribute(groupTitanSimulation, restart_prefix);
    TiH5_readBoolAttribute(groupTitanSimulation, restart_keep_all);
    TiH5_readBoolAttribute(groupTitanSimulation, restart_keep_redundant_data);
    TiH5_readIntAttribute(groupTitanSimulation, adapt);
    TiH5_readScalarDataTypeAttribute(groupTitanSimulation, elementType, datatypeElementType);

    ///////////////////////////////////////////////////////////////////////////
    NodeTable=new NodeHashTable(parent);
    //@TODO do something about bc
    ElemTable=new ElementsHashTable(NodeTable,parent);

    timeprops.h5read(parent);
    scale_.h5read(parent);

    pileprops=PileProps::createPileProps(parent);
    fluxprops.h5read(parent);
    discharge_planes.h5read(parent);
    matprops=MatProps::createMatProps(parent,scale_);
    statprops=new StatProps(ElemTable,NodeTable,parent);
    mapnames.h5read(parent);
    outline.setElemNodeTable(ElemTable,NodeTable);
    outline.h5read(parent);

    integrator=Integrator::createIntegrator(parent,this);

}

void xmdfScalarAttribute(fstream &xmdf,const char *name, const char *hf5_filename, const char *h5ref, const int dim,
        const char *DataType,const int Precision, const char *Center)
{
    xmdf<< "\t\t<Attribute Type=\"Scalar\" Center=\""<<Center<<"\" Name=\""<<name<<"\">\n";
    xmdf<< "\t\t\t<DataItem DataType=\""<<DataType<<"\" Precision=\""<<Precision<<"\" Dimensions=\""<<dim<<" 1\" Format=\"HDF\">\n";
    xmdf<< "\t\t\t\t\t"<<hf5_filename<<":"<<h5ref<<"\n";
    xmdf<< "\t\t\t</DataItem>\n";
    xmdf<< "\t\t</Attribute>\n";
}
void xmdfVectorAttributeXY(fstream &xmdf,const char *name, const char *hf5_filename, const char *h5refX, const char *h5refY, const int dim,
        const char *DataType,const int Precision, const char *Center)
{
    xmdf<< "\t\t<Attribute Type=\"Vector\" Center=\""<<Center<<"\" Name=\""<<name<<"\">\n";
    xmdf<< "\t\t\t<DataItem ItemType=\"Function\" Function=\"JOIN($0,$1,0*$1)\" Dimensions=\""<<dim<<" 3\">\n";
    xmdf<< "\t\t\t\t<DataItem DataType=\""<<DataType<<"\" Precision=\""<<Precision<<"\" Dimensions=\""<<dim<<"\" Format=\"HDF\">\n";
    xmdf<< "\t\t\t\t\t"<<hf5_filename<<":"<<h5refX<<"\n";
    xmdf<< "\t\t\t\t</DataItem>\n";
    xmdf<< "\t\t\t\t<DataItem DataType=\""<<DataType<<"\" Precision=\""<<Precision<<"\" Dimensions=\""<<dim<<"\" Format=\"HDF\">\n";
    xmdf<< "\t\t\t\t\t"<<hf5_filename<<":"<<h5refY<<"\n";
    xmdf<< "\t\t\t\t</DataItem>\n";
    xmdf<< "\t\t\t</DataItem>\n";
    xmdf<< "\t\t</Attribute>\n";
}
void xmdfVectorAttributeXYZ(fstream &xmdf,const char *name, const char *hf5_filename, const char *h5refX, const char *h5refY, const char *h5refZ, const int dim,
        const char *DataType,const int Precision, const char *Center)
{
    xmdf<< "\t\t<Attribute Type=\"Vector\" Center=\""<<Center<<"\" Name=\""<<name<<"\">\n";
    xmdf<< "\t\t\t<DataItem ItemType=\"Function\" Function=\"JOIN($0,$1,$2)\" Dimensions=\""<<dim<<" 3\">\n";
    xmdf<< "\t\t\t\t<DataItem DataType=\""<<DataType<<"\" Precision=\""<<Precision<<"\" Dimensions=\""<<dim<<"\" Format=\"HDF\">\n";
    xmdf<< "\t\t\t\t\t"<<hf5_filename<<":"<<h5refX<<"\n";
    xmdf<< "\t\t\t\t</DataItem>\n";
    xmdf<< "\t\t\t\t<DataItem DataType=\""<<DataType<<"\" Precision=\""<<Precision<<"\" Dimensions=\""<<dim<<"\" Format=\"HDF\">\n";
    xmdf<< "\t\t\t\t\t"<<hf5_filename<<":"<<h5refY<<"\n";
    xmdf<< "\t\t\t\t</DataItem>\n";
    xmdf<< "\t\t\t\t<DataItem DataType=\""<<DataType<<"\" Precision=\""<<Precision<<"\" Dimensions=\""<<dim<<"\" Format=\"HDF\">\n";
    xmdf<< "\t\t\t\t\t"<<hf5_filename<<":"<<h5refZ<<"\n";
    xmdf<< "\t\t\t\t</DataItem>\n";
    xmdf<< "\t\t\t</DataItem>\n";
    xmdf<< "\t\t</Attribute>\n";
}
void xmdfNScalarAttributes(fstream &xmdf,const char *name_prefix, const char *hf5_filename, const char *h5ref_prefix, const int dim0, const int dim,
        const char *DataType,const int Precision, const char *Center, vector<string> *comp_names=NULL)
{
    for( int i=0;i<dim0;i++)
    {
        string name=name_prefix;
        string h5ref=h5ref_prefix;
        if(comp_names==NULL)
        {
            name+=to_string(i);
            h5ref+=to_string(i);
        }
        else
        {
            name+=(*comp_names)[i];
            h5ref+=(*comp_names)[i];
        }
        xmdfScalarAttribute(xmdf,name.c_str(), hf5_filename, h5ref.c_str(), dim,DataType,Precision, Center);
    }
}
void cxxTitanSimulation::xmdfWrite(const char *xmdf_filename,const char *hf5_filename,bool Quad9)
{
    //write xmdf part
    vector<string> *ElementTypesVarNames=nullptr;
    if(elementType==ElementType::SinglePhase)
        ElementTypesVarNames=&SinglePhaseVarNames;
    else if(elementType==ElementType::TwoPhases)
        ElementTypesVarNames=&TwoPhasesVarNames;

    //stringstream xmdf;
    fstream xmdf(xmdf_filename, std::fstream::in | std::fstream::out | std::fstream::ate);
    if(xmdf.good()==0||xmdf.tellg()<60)
    {
        //file do not exists write header
        xmdf.open(xmdf_filename, std::fstream::out);
        xmdf<<"<?xml version=\"1.0\" ?>\n";
        xmdf<<"<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\"[]> \n";
        xmdf<<"<Xdmf Version=\"3.0\">\n";
        xmdf<<"<Domain>\n";
        xmdf<<"<Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">\n";
    }
    else
    {
        //file exists rewind to beginning of footer
        xmdf.seekg (-60, xmdf.end);
        string line;
        std::streampos pos=xmdf.tellg();
        while(std::getline(xmdf,line))
        {
            line.erase(line.find_last_not_of(" \n\r\t")+1);
            line.erase(0, line.find_first_not_of(" \t"));
            if(line=="<!-- Footer Starts Here -->")
            {
                xmdf.seekg (pos);
                xmdf.seekp (pos);
                break;
            }
            pos=xmdf.tellg();
        }
    }

    int NumberOfElements=ElemTable->size();
    int NumberOfNodes=NodeTable->size();

    //write body
    xmdf<<"\t<Grid Name=\"TimeSeriesA\">\n";
    xmdf<<"\t\t<Time Value=\""<<timeprops.timesec()<<"\"/>\n";


    //Topology
    if(Quad9)
    {
        xmdf<<"\t\t<Topology TopologyType=\"Quad_9\" NodesPerElement=\"9\" NumberOfElements=\""<<NumberOfElements<<"\">\n";
        xmdf<<"\t\t\t<DataItem Name=\"Connections\" ItemType=\"Function\" Function=\"JOIN($0,$1,$2,$3,$4,$5,$6,$7,$8)\" Dimensions=\""<<NumberOfElements<<" 9\">\n";
        for( int i=0;i<8;i++)
        {
            xmdf<<"\t\t\t\t<DataItem DataType=\"Int\" Precision=\"4\" Dimensions=\""<<NumberOfElements<<"\" Format=\"HDF\">\n";
            xmdf<<"\t\t\t\t\t"<<hf5_filename<<":/ElemTable/node_key_ndx_"<<i<<"\n";
            xmdf<<"\t\t\t\t</DataItem>\n";
        }
        xmdf<<"\t\t\t\t<DataItem DataType=\"Int\" Precision=\"4\" Dimensions=\""<<NumberOfElements<<"\" Format=\"HDF\">\n";
        xmdf<<"\t\t\t\t\t"<<hf5_filename<<":/ElemTable/node_bubble_ndx_\n";
        xmdf<<"\t\t\t\t</DataItem>\n";
        xmdf<<"\t\t\t</DataItem>\n";
        xmdf<<"\t\t</Topology>\n";
    }
    else
    {
        xmdf<<"\t\t<Topology TopologyType=\"Quadrilateral\" NodesPerElement=\"4\" NumberOfElements=\""<<NumberOfElements<<"\" Order=\"0 3 2 1\">\n";

        xmdf<<"\t\t\t<DataItem Name=\"Connections\" ItemType=\"Function\" Function=\"JOIN($0,$1,$2,$3)\" Dimensions=\""<<NumberOfElements<<" 4\">\n";
        for( int i=0;i<4;i++)
        {
            xmdf<<"\t\t\t\t<DataItem DataType=\"Int\" Precision=\"4\" Dimensions=\""<<NumberOfElements<<"\" Format=\"HDF\">\n";
            xmdf<<"\t\t\t\t\t"<<hf5_filename<<":/ElemTable/node_key_ndx_"<<i<<"\n";
            xmdf<<"\t\t\t\t</DataItem>\n";
        }
        xmdf<<"\t\t\t</DataItem>\n";
        xmdf<<"\t\t</Topology>\n";
    }

    //Geometry
    xmdf<<"\t\t<Geometry Type=\"XYZ\">\n";
    xmdf<< "\t\t\t<DataItem Name=\"Coordinates\" ItemType=\"Function\" Function=\"JOIN($0,$1,$2)\" Dimensions=\""<<NumberOfNodes<<" 3\">\n";
    xmdf<< "\t\t\t\t<DataItem Name=\"X\" DataType=\"Float\" Precision=\"8\" Dimensions=\""<<NumberOfNodes<<"\" Format=\"HDF\">\n";
    xmdf<< "\t\t\t\t\t"<<hf5_filename<<":/NodeTable/coord_X\n";
    xmdf<< "\t\t\t\t</DataItem>\n";
    xmdf<< "\t\t\t\t<DataItem Name=\"Y\" DataType=\"Float\" Precision=\"8\" Dimensions=\""<<NumberOfNodes<<"\" Format=\"HDF\">\n";
    xmdf<< "\t\t\t\t\t"<<hf5_filename<<":/NodeTable/coord_Y\n";
    xmdf<< "\t\t\t\t</DataItem>\n";
    xmdf<< "\t\t\t\t<DataItem Name=\"Z\" DataType=\"Float\" Precision=\"8\" Dimensions=\""<<NumberOfNodes<<"\" Format=\"HDF\">\n";
    xmdf<< "\t\t\t\t\t"<<hf5_filename<<":/NodeTable/elevation_\n";
    xmdf<< "\t\t\t\t</DataItem>\n";
    xmdf<< "\t\t\t</DataItem>\n";
    xmdf<< "\t\t</Geometry>\n";
    //Attributes
    xmdfScalarAttribute(xmdf,"myprocess", hf5_filename, "/ElemTable/myprocess_", NumberOfElements, "Int",4, "Cell");
    xmdfScalarAttribute(xmdf,"generation", hf5_filename, "/ElemTable/generation_", NumberOfElements, "Int",4, "Cell");
    xmdfScalarAttribute(xmdf,"material", hf5_filename, "/ElemTable/material_", NumberOfElements, "Int",4, "Cell");

    xmdfScalarAttribute(xmdf,"shortspeed", hf5_filename, "/ElemTable/shortspeed_", NumberOfElements, "Float",8, "Cell");
    xmdfScalarAttribute(xmdf,"positive_x_side", hf5_filename, "/ElemTable/positive_x_side_", NumberOfElements, "Int",4, "Cell");

    xmdfScalarAttribute(xmdf,"stoppedflags", hf5_filename, "/ElemTable/stoppedflags_", NumberOfElements, "Int",4, "Cell");

    xmdfScalarAttribute(xmdf,"effect_bedfrict_", hf5_filename, "/ElemTable/shortspeed_", NumberOfElements, "Float",8, "Cell");
    xmdfScalarAttribute(xmdf,"effect_tanbedfrict_", hf5_filename, "/ElemTable/shortspeed_", NumberOfElements, "Float",8, "Cell");

    xmdfNScalarAttributes(xmdf,"d_state_vars_", hf5_filename, "/ElemTable/d_state_vars_",NUM_STATE_VARS * DIMENSION, NumberOfElements, "Float",8, "Cell");

    xmdfVectorAttributeXY(xmdf,"eigenvxymax_", hf5_filename, "/ElemTable/eigenvxymax_X", "/ElemTable/eigenvxymax_Y", NumberOfElements, "Float",8, "Cell");
    xmdfVectorAttributeXY(xmdf,"kactxy_", hf5_filename, "/ElemTable/kactxy_X", "/ElemTable/kactxy_Y", NumberOfElements, "Float",8, "Cell");
    xmdfVectorAttributeXY(xmdf,"zeta_", hf5_filename, "/ElemTable/zeta_X", "/ElemTable/zeta_Y", NumberOfElements, "Float",8, "Cell");
    xmdfVectorAttributeXY(xmdf,"curvature_", hf5_filename, "/ElemTable/curvature_X", "/ElemTable/curvature_Y", NumberOfElements, "Float",8, "Cell");
    xmdfVectorAttributeXYZ(xmdf,"gravity_", hf5_filename, "/ElemTable/gravity_X", "/ElemTable/gravity_Y", "/ElemTable/gravity_Z", NumberOfElements, "Float",8, "Cell");
    xmdfVectorAttributeXY(xmdf,"d_gravity_", hf5_filename, "/ElemTable/d_gravity_X", "/ElemTable/d_gravity_Y", NumberOfElements, "Float",8, "Cell");

    xmdfVectorAttributeXY(xmdf,"effect_kactxy_", hf5_filename, "/ElemTable/effect_kactxy_0", "/ElemTable/effect_kactxy_1", NumberOfElements, "Float",8, "Cell");

    xmdfScalarAttribute(xmdf,"ithelem_", hf5_filename, "/ElemTable/ithelem_", NumberOfElements, "Int",4, "Cell");
    xmdfScalarAttribute(xmdf,"iwetnode_", hf5_filename, "/ElemTable/iwetnode_", NumberOfElements, "Int",4, "Cell");
    xmdfScalarAttribute(xmdf,"Awet_", hf5_filename, "/ElemTable/Awet_", NumberOfElements, "Float",8, "Cell");
    xmdfVectorAttributeXY(xmdf,"drypoint_", hf5_filename, "/ElemTable/drypoint_0", "/ElemTable/drypoint_1", NumberOfElements, "Float",8, "Cell");
    xmdfScalarAttribute(xmdf,"Swet_", hf5_filename, "/ElemTable/Swet_", NumberOfElements, "Float",8, "Cell");

    xmdfNScalarAttributes(xmdf,"el_error_", hf5_filename, "/ElemTable/el_error_",EQUATIONS, NumberOfElements, "Float",8, "Cell");
    xmdfNScalarAttributes(xmdf,"el_solution_", hf5_filename, "/ElemTable/el_solution_",EQUATIONS, NumberOfElements, "Float",8, "Cell");

    xmdfNScalarAttributes(xmdf,"Influx_", hf5_filename, "/ElemTable/Influx_",NUM_STATE_VARS, NumberOfElements, "Float",8, "Cell",ElementTypesVarNames);
    xmdfNScalarAttributes(xmdf,"which_son_", hf5_filename, "/ElemTable/which_son_",NUM_STATE_VARS, NumberOfElements, "Int",8, "Cell",ElementTypesVarNames);

    if(Quad9)
    {
        xmdfScalarAttribute(xmdf,"node_info_", hf5_filename, "/NodeTable/info_", NumberOfNodes, "Int",4, "Node");
        xmdfNScalarAttributes(xmdf,"flux_", hf5_filename, "/NodeTable/flux_",NUM_STATE_VARS, NumberOfNodes, "Float",8, "Node",ElementTypesVarNames);
    }

    xmdfScalarAttribute(xmdf,"key_", hf5_filename, "/ElemTable/key_", NumberOfElements, "UInt",8, "Cell");

    if(elementType==ElementType::SinglePhase)
    {
        xmdfScalarAttribute(xmdf,"h", hf5_filename, "/ElemTable/state_vars_h", NumberOfElements, "Float",8, "Cell");
        xmdfVectorAttributeXY(xmdf,"hV", hf5_filename, "/ElemTable/state_vars_hVx", "/ElemTable/state_vars_hVy", NumberOfElements, "Float",8, "Cell");
    }
    else if(elementType==ElementType::TwoPhases)
    {
        xmdfNScalarAttributes(xmdf,"state_vars_", hf5_filename, "/ElemTable/state_vars_",EQUATIONS, NumberOfElements, "Float",8, "Cell",ElementTypesVarNames);
    }
    /*xmdfScalarAttribute(xmdf,"h", hf5_filename, "/ElemTable/state_vars_h", NumberOfElements, "Float",8, "Cell");
    xmdfScalarAttribute(xmdf,"hVx", hf5_filename, "/ElemTable/state_vars_hVx", NumberOfElements, "Float",8, "Cell");
    xmdfScalarAttribute(xmdf,"hVy", hf5_filename, "/ElemTable/state_vars_hVy", NumberOfElements, "Float",8, "Cell");*/


    xmdf<<"\t</Grid>\n";

    //footer
    xmdf<<"<!-- Footer Starts Here -->\n";
    xmdf<<"</Grid>\n";
    xmdf<<"</Domain>\n";
    xmdf<<"</Xdmf>\n";

    xmdf.close();

}



void cxxTitanSimulation::save_restart_file()
{
    char hf5_filename[256];
    char hf5_snapshot[256];
    char xmdf_Quad9_filename[256];
    char xmdf_Quad4_filename[256];
    char xmdf_snapshot[256];

    sprintf(hf5_filename,  "%s/snapshot_p%04d_i%08d.h5", restart_prefix.c_str(), myid, timeprops.iter);
    sprintf(hf5_snapshot,  "snapshot_p%04d_i%08d.h5", myid, timeprops.iter);
    sprintf(xmdf_snapshot, "%s/snapshot_Quad9_p%04d_i%08d.xmf", restart_prefix.c_str(), myid, timeprops.iter);
    sprintf(xmdf_Quad9_filename, "%s_Quad9_p%04d.xmf", restart_prefix.c_str(), myid);
    sprintf(xmdf_Quad4_filename, "%s_Quad4_p%04d.xmf", restart_prefix.c_str(), myid);

    //write H5
    H5::H5File file(hf5_filename, H5F_ACC_TRUNC);
    h5write(&file);
    file.close();

    if(restart_keep_all)
    {
        //@todo writes xdmf file in a such way that it is complete after each snaphsot
        xmdfWrite(xmdf_Quad9_filename,hf5_filename,true);
        xmdfWrite(xmdf_Quad4_filename,hf5_filename,false);
    }

    xmdfWrite(xmdf_snapshot,hf5_snapshot,true);

    if(restart_keep_all==false)
    {
        remove(prev_hf5_filename.c_str());
        remove(prev_xmdf_snapshot.c_str());
    }

    prev_hf5_filename=hf5_filename;
    prev_xmdf_snapshot=xmdf_snapshot;
}
void cxxTitanSimulation::load_restart(const std::string restartFilename, const std::string new_gis_main)
{
    // Create a new file using the default property lists.
    H5::H5File file(restartFilename.c_str(), H5F_ACC_RDONLY);

    h5read(&file);

    file.close();

    //change location of DEM
    if(new_gis_main!="")
    {
        if(mapnames.gis_format==MapNames::GIS_GRASS)
        {
            mapnames.gis_main=new_gis_main;
        }
        else
        {
            mapnames.gis_map=new_gis_main;
        }
    }

    //check output directory and create if necessary
    if(restart_prefix!="" && ti_dir_exists(restart_prefix.c_str())==false)
    {
        ti_mkdir(restart_prefix.c_str());
    }
    if(vizoutput_prefix!="" && ti_dir_exists(vizoutput_prefix.c_str())==false)
    {
        ti_mkdir(vizoutput_prefix.c_str());
    }
}
void cxxTitanSimulation::save_vizoutput_file(const int mode)
{
    int savefileflag=1;
    if(mode!=XDMF_ONLYCLOSE)
    {
        move_data(numprocs, myid, ElemTable, NodeTable, &timeprops);

        output_discharge(matprops, &timeprops, &discharge_planes, myid);

        if(myid == 0)
            output_summary(&timeprops, statprops, savefileflag);

        if(vizoutput & 1)
            tecplotter(elementType, ElemTable, NodeTable, matprops, &timeprops, &mapnames, statprops->vstar,vizoutput_prefix.c_str());

        if(vizoutput & 2)
            meshplotter(ElemTable, NodeTable, matprops, &timeprops, &mapnames, statprops->vstar,vizoutput_prefix.c_str());

        if(vizoutput & 8)
        {
            if(myid == 0)
                grass_sites_header_output(&timeprops,vizoutput_prefix.c_str());
            grass_sites_proc_output(ElemTable, NodeTable, myid, matprops, &timeprops,vizoutput_prefix.c_str());
        }
    }
#if HAVE_LIBHDF5
    if(vizoutput & 4)
    {
        if(elementType == ElementType::TwoPhases)
        {
            write_xdmf_two_phases(ElemTable, NodeTable, &timeprops,
                    matprops, &mapnames, mode,vizoutput_prefix.c_str());
        }
        if(elementType == ElementType::SinglePhase)
        {
            write_xdmf_single_phase(ElemTable, NodeTable, &timeprops,
                    matprops, &mapnames, mode,vizoutput_prefix.c_str());
        }
    }
#endif

    if (timeprops.iVizOutSaved > 0 && outline.enabled)
    {
    	outline.combine_results_from_threads();
//    	outline.outputST(matprops, &timeprops);
    }

    if (timeprops.iVizOutSaved > 0){

    	double temp_A = 0.0, Temp_Vol = 0.0, Temp_Vel = 0.0;
        double temp_Fx [5] = {}, temp_Fy [5] = {}, temp_FV [5] = {}, temp_P [5] = {};
        for (int i=0; i < 5; i++)
        {
        	temp_Fx[i] = 0.0;
        	temp_Fy[i] = 0.0;
        	temp_FV[i] = 0.0;
        	temp_P[i] = 0.0;
        }
        double Vx, Vy, Vm, dV, dA;
//        double teta = 38.5;

    	int no_of_buckets = ElemTable->get_no_of_buckets();
    	vector<HashEntryLine> &bucket = ElemTable->bucket;
    	tivector<Element> &elenode_ = ElemTable->elenode_;

    	//@ElementsBucketDoubleLoop
    	for (int ibuck = 0; ibuck < no_of_buckets; ibuck++) {
    		for (int ielm = 0; ielm < bucket[ibuck].ndx.size(); ielm++) {
    			Element* Curr_El = &(elenode_[bucket[ibuck].ndx[ielm]]);

			if (Curr_El->adapted_flag() > 0 && Curr_El->state_vars(0)*scale_.height >= 0.05) {

				dA = Curr_El->dx(0) * Curr_El->dx(1);
				dV = Curr_El->state_vars(0) * dA;
				Temp_Vol += dV;
				Temp_Vel += dA * sqrt(Curr_El->state_vars(1)*Curr_El->state_vars(1)+Curr_El->state_vars(2)*Curr_El->state_vars(2));


//			}
//    			if (Curr_El->adapted_flag() > 0 && Curr_El->state_vars(0)*scale_.height >= 0.05) {

    				Vx = Curr_El->state_vars(1)/Curr_El->state_vars(0);
    				Vy = Curr_El->state_vars(2)/Curr_El->state_vars(0);

    				Vm = sqrt(Vx*Vx+Vy*Vy);

    				temp_Fx[0] += Curr_El->Drivs(0);
    				temp_Fy[0] += Curr_El->Drivs(2);
    				temp_Fx[1] += Curr_El->Drivs(1);
    				temp_Fy[1] += Curr_El->Drivs(3);
        			temp_Fx[2] += Curr_El->STs(0);
        			temp_Fy[2] += Curr_El->STs(NUM_STATE_VARS);
        			temp_Fx[3] += Curr_El->STs(1);
        			temp_Fy[3] += Curr_El->STs(1 + NUM_STATE_VARS);
        			temp_Fx[4] += Curr_El->STs(2);
        			temp_Fy[4] += Curr_El->STs(2 + NUM_STATE_VARS);

        			temp_P[0] += Curr_El->Drivs(0)*Vx + Curr_El->Drivs(2)*Vy;
        			temp_P[1] += Curr_El->Drivs(1)*Vx + Curr_El->Drivs(3)*Vy;
        			temp_P[2] += Curr_El->STs(0)*Vx + Curr_El->STs(NUM_STATE_VARS)*Vy;
        			temp_P[3] += Curr_El->STs(1)*Vx + Curr_El->STs(1 + NUM_STATE_VARS)*Vy;
        			temp_P[4] += Curr_El->STs(2)*Vx + Curr_El->STs(2 + NUM_STATE_VARS)*Vy;

        			if (Vm > 0.0) {
            				temp_FV[0] += temp_P[0] / Vm;
            				temp_FV[1] += temp_P[1] / Vm;
            				temp_FV[2] += temp_P[2] / Vm;
            				temp_FV[3] += temp_P[3] / Vm;
            				temp_FV[4] += temp_P[4] / Vm;
        			}

       				temp_A += Curr_El->dx(0)*Curr_El->dx(1);
    			}
    		}
    	}

    	Area.push_back(temp_A);
    	Finertx.push_back(temp_Fx[0]);
    	Finerty.push_back(temp_Fy[0]);
    	Fconvx.push_back(temp_Fx[1]);
    	Fconvy.push_back(temp_Fy[1]);
        Fgravx.push_back(temp_Fx[2]);
        Fgravy.push_back(temp_Fy[2]);
        Fbx.push_back(temp_Fx[3]);
        Fby.push_back(temp_Fy[3]);
        Fix.push_back(temp_Fx[4]);
        Fiy.push_back(temp_Fy[4]);

        Pinert.push_back(temp_P[0]);
        Pconv.push_back(temp_P[1]);
        Pgrav.push_back(temp_P[2]);
        Pbed.push_back(temp_P[3]);
        Pint.push_back(temp_P[4]);

        FVinert.push_back(temp_FV[0]);
        FVconv.push_back(temp_FV[1]);
        FVgrav.push_back(temp_FV[2]);
        FVbed.push_back(temp_FV[3]);
        FVint.push_back(temp_FV[4]);

        Volume.push_back(Temp_Vol);
        MVelocity.push_back(Temp_Vel/Temp_Vol);
    }

}

void cxxTitanSimulation::run(bool start_from_restart)
{
    IF_MPI(MPI_Barrier (MPI_COMM_WORLD));

    int i; //-- counters
    int current_state_is_good=true;

    //-- MPI
    IF_MPI(MPI_Status status);

    TIMING1_DEFINE(start);
    TIMING1_DEFINE(end);
    TIMING1_DEFINE(t_start);
    TIMING1_DEFINE(t_start2);
    TIMING1_DEFINE(t_end);

    TIMING1_START(start);

    /* create new MPI datastructures for class objects */
    IF_MPI(MPI_New_Datatype());

    /* read original data from serial preprocessing
     code and then initialize element
     stiffness routines info */

    int xdmerr;



    double end_time = 10000.0;
    /*
     * viz_flag is used to determine which viz output to use
     * nonzero 1st bit of viz_flag means output tecplotxxxx.tec
     * nonzero 2nd bit of viz_flag means output mshplotxxxx.tec (debug purposes)
     * nonzero 3rd bit of viz_flag means output Paraview/XDMF format
     * nonzero 4th bit of viz_flag means output grass_sites files

     order_flag == 1 means use first order method
     order_flag == 2 means use second order method
     */

    //savefileflag will be flipped so first savefile will end in 0
    int savefileflag = 1;
    int Init_Node_Num, Init_Elem_Num;
    double v_star; // v/v_slump
    double nz_star; /* temporary... used for negligible velocity as stopping
     criteria paper... plan to include in v_star implicitly
     later */

    MatProps *matprops_ptr=get_matprops();
    PileProps* pileprops_ptr=get_pileprops();
    StatProps* statprops_ptr=get_statprops();
    TimeProps* timeprops_ptr=get_timeprops();
    MapNames* mapnames_ptr=get_mapnames();
    OutLine* outline_ptr=get_outline();

    process_input(start_from_restart);

    if(start_from_restart==false)
        input_summary();


    if(start_from_restart==false)
        Read_grid(myid, numprocs, &NodeTable,&ElemTable, matprops_ptr, &outline);
    
    

    if(start_from_restart==false)
    {
        setup_geoflow(ElemTable, NodeTable, myid, numprocs, matprops_ptr, &timeprops);

        move_data(numprocs, myid, ElemTable, NodeTable, &timeprops);

        AssertMeshErrorFree(ElemTable, NodeTable, numprocs, myid, -1.0);

        //initialize pile height and if appropriate perform initial adaptation
        init_piles();
    }
    else
    {
    	NodeTable->flushNodeTable();
        ElemTable->flushElemTable();
        //update temporary arrays of elements/nodes pointers
        ElemTable->updateLocalElements();
        ElemTable->updateNeighboursIndexes();
    }
    
    ElementsProperties *ElemProp=ElemTable->ElemProp;

    HAdapt hadapt(ElemTable, NodeTable, ElemProp,&timeprops,matprops_ptr,5);
    HAdaptUnrefine Unrefine(ElemTable, NodeTable, &timeprops,matprops_ptr);

    outline.setElemNodeTable(ElemTable,NodeTable);

    assert(integrator!=nullptr);

    /* for debug only, to check if exactly what's loaded will be saved again
     by doing a diff on the files.
     saverun(&BT_Node_Ptr, myid, numprocs, &BT_Elem_Ptr,
     matprops_ptr, &timeprops, &mapnames, adaptflag, order_flag,
     &statprops, &discharge_planes, &outline, &savefileflag);
     */
    if(myid == 0)
    {
        for(int imat = 1; imat <= matprops_ptr->material_count; imat++)
            printf("bed friction angle for \"%s\" is %g\n", matprops_ptr->matnames[imat].c_str(),
                   matprops_ptr->bedfrict[imat] * 180.0 / PI);

        printf("internal friction angle is %g, epsilon is %g \n method order = %i\n", integrator->int_frict * 180.0 / PI,
               scale_.epsilon, integrator->order);
        printf("REFINE_LEVEL=%d\n", REFINE_LEVEL);
    }

    IF_MPI(MPI_Barrier (MPI_COMM_WORLD));
    statprops->calc_stats(myid, matprops_ptr, &timeprops, &discharge_planes, 0.0);

    output_discharge(matprops_ptr, &timeprops, &discharge_planes, myid);

    move_data(numprocs, myid, ElemTable, NodeTable, &timeprops);
    save_vizoutput_file(XDMF_NEW);

    /*
     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

     Time Stepping Loop

     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     */

    long element_counter = 0; // for performance count elements/timestep/proc
    int ifstop = 0;
    double max_momentum = 100;  //nondimensional

    if(restart_enabled)
    {
        if(start_from_restart==false)
            save_restart_file();
    }

    /* ifend(0.5*statprops.vmean) is a hack, the original intent (when we were
     intending to use vstar as a stopping criteria) whas to have the
     calculation when vstar dropped back down below 1, instead we're
     using the ifend() function to stop the simulation when the volume
     averaged velocity falls back down below 2 meters... this hack is only
     for the colima hazard map runs, otherwise pass ifend() a constant
     valued */

    ASSERT2(ElemTable->checkPointersToNeighbours("Prestep index check",false)==0);

    PROFILING1_DEFINE(pt_start0);
    PROFILING1_DEFINE(pt_start1);
    PROFILING1_START(pt_start0);
    PROFILING1_DEFINE(pt_start);
    PROFILING1_START(pt_start);

    T_Finertx=0.0;
    T_Finerty=0.0;
    T_Fconvx=0.0;
    T_Fconvy=0.0;
    T_Fgravx=0.0;
    T_Fgravy=0.0;
    T_Fbx=0.0;
    T_Fby=0.0;
    T_Fix=0.0;
    T_Fiy=0.0;
    T_Pinert=0.0;
    T_Pconv=0.0;
    T_Pgrav=0.0;
    T_Pbed=0.0;
    T_Pint=0.0;

    titanTimingsAlongSimulation.totalTime = MPI_Wtime();
    while (!(timeprops.ifend(0)) && !ifstop)
    {
        /*
         *  mesh adaption routines
         */
        TIMING1_START(t_start);
        PROFILING1_START(pt_start);
        double TARGET = .05;
        double UNREFINE_TARGET = .01;
        int h_count = 0;
        if(timeprops.iter < 50)
        {
            scale_.frict_tiny = 0.1;
            integrator->frict_tiny=scale_.frict_tiny;
        }
        else
        {
            scale_.frict_tiny = 0.000000001;
            integrator->frict_tiny=scale_.frict_tiny;
        }

        //check for changes in topography and update if necessary
        if(timeprops.iter == 200)
        {
            update_topo(ElemTable, NodeTable, myid, numprocs, matprops_ptr, &timeprops, &mapnames);
        }
        PROFILING1_STOPADD_RESTART(tsim_iter_update_topo,pt_start);

        if((adapt != 0) && (timeprops.iter % 5 == 4))
        {
            AssertMeshErrorFree(ElemTable, NodeTable, numprocs, myid, -2.0);

            TIMING1_START(t_start2);
            ElemTable->conformation++;
            hadapt.adapt(h_count, TARGET);
            TIMING1_STOPADD(refinementTime, t_start2);
            
            
            TIMING1_START(t_start2);
            Unrefine.unrefine(UNREFINE_TARGET);

            //this move_data() here for debug... to make AssertMeshErrorFree() Work
            move_data(numprocs, myid, ElemTable, NodeTable, &timeprops);
            TIMING1_STOPADD(unrefinementTime, t_start2);

            if((numprocs > 1) && (timeprops.iter % 10 == 9))
            {
                repartition2(ElemTable, NodeTable, &timeprops);

                //this move_data() here for debug... to make AssertMeshErrorFree() Work
                move_data(numprocs, myid, ElemTable, NodeTable, &timeprops);

                ElemTable->updateLocalElements();
                ElemTable->updateNeighboursIndexes();
            }


            //update temporary arrays of elements/nodes pointers
            NodeTable->flushNodeTable();
            ElemTable->flushElemTable();
            
            ASSERT2(ElemTable->checkPointersToNeighbours("After all adaptions",false)==0);
        }
        TIMING1_STOPADD(meshAdaptionTime, t_start);
        PROFILING1_STOPADD_RESTART(tsim_iter_adapt,pt_start);

        TIMING1_START(t_start);
        integrator->step();
        //step(matprops_ptr, &timeprops, pileprops_ptr, &fluxprops, &statprops,
        //     &order, &outline, &discharge_planes, adapt);
        TIMING1_STOPADD(stepTime, t_start);
        PROFILING1_STOPADD_RESTART(tsim_iter_step,pt_start);

        //check if current state is good
        if(myid == 0)
        {
			if(statprops_ptr->statvolume==0.0){
				current_state_is_good=false;

				printf("\n###############################################################################\n");
				printf("Volume is zero, nothing more to model!\n");
				printf("Be careful with stats, some nans can meke it there!\n");
				printf("###############################################################################\n\n");

				break;
			}
			if(isnan(statprops_ptr->vmean) || isnan(statprops_ptr->vstar))
			{
				current_state_is_good=false;

				printf("\n###############################################################################\n");
				printf("Velocity is nan(not a number), something is wrong!\n");
				printf("Do not forget about restarts, the last one should be without nan.\n");
				printf("Be careful with stats, some nans can meke it there!\n");
				printf("###############################################################################\n\n");

				break;
			}
			if(isnan(statprops_ptr->statvolume))
			{
				current_state_is_good=false;

				printf("\n###############################################################################\n");
				printf("Volume is nan(not a number), something is wrong!\n");
				printf("Do not forget about restarts, the last one should be without nan.\n");
				printf("Be careful with stats, some nans can meke it there!\n");
				printf("###############################################################################\n\n");

				break;
			}
			IF_MPI(MPI_Bcast(&current_state_is_good,1,MPI_INT,0,MPI_COMM_WORLD));
        }
        else
        {
        	IF_MPI(MPI_Bcast(&current_state_is_good,1,MPI_INT,0,MPI_COMM_WORLD));
        }
        if(current_state_is_good==false)
        {
        	break;
        }
        TIMING1_START(t_start);
        /*
         * save a restart file
         */
        if(restart_enabled && timeprops.ifTimeForRestartOutput())
        {
            save_restart_file();
            //saverun(&NodeTable, myid, numprocs, &ElemTable, matprops_ptr, &timeprops, &mapnames, adapt, integrator->order,
            //        statprops, &discharge_planes, &outline, &savefileflag);
        }

        PROFILING1_STOPADD_RESTART(tsim_iter_saverestart,pt_start);
        /*
         * output results to file
         */
        if(timeprops.ifTimeForTimeSeriesOutput())
            save_vizoutput_file(XDMF_OLD);
        PROFILING1_STOPADD_RESTART(tsim_iter_output,pt_start);
#ifdef PERFTEST
        int countedvalue=timeprops.iter%2+1;
        int e_buckets=ElemTable->get_no_of_buckets();
        HashEntry* entryp;
        for(i=0; i<e_buckets; i++)
        {
            entryp = *(ElemTable->getbucketptr() + i);
            while(entryp)
            {
                Element * EmTemp = (Element*)entryp->value;
                assert(EmTemp);
                assert(EmTemp->get_counted()!=countedvalue);

                if((EmTemp->get_adapted_flag()>=NOTRECADAPTED)&&
                        (EmTemp->get_adapted_flag()<=BUFFER)
                )
                {
                    //if this element doesn't belong on this processor don't involve
                    element_counter++;
                    EmTemp->put_counted(countedvalue);
                }
                entryp = entryp->next;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        TIMING1_STOPADD(resultsOutputTime, t_start);

        if(timeprops.iter % 200 == 0)
        {
            titanTimingsAlongSimulation.totalTime = MPI_Wtime() - titanTimingsAlongSimulation.totalTime;
            if(myid == 0)
                titanTimingsAlongSimulation.print("Timings over last 200 steps (seconds):");
            titanTimingsAlongSimulation.reset();
            titanTimingsAlongSimulation.totalTime = MPI_Wtime();
        }
        PROFILING1_STOPADD_RESTART(tsim_iter_post,pt_start);

    	if((adapt != 0) && (timeprops.iter % 5 == 4))
    	{
    		N_elements.push_back(ElemTable->get_no_of_entries());
    		Time_el.push_back(timeprops.cur_time * timeprops.TIME_SCALE);
    	}


        double Vx, Vy, Vm;

    	int no_of_buckets = ElemTable->get_no_of_buckets();
    	vector<HashEntryLine> &bucket = ElemTable->bucket;
    	tivector<Element> &elenode_ = ElemTable->elenode_;

    	//@ElementsBucketDoubleLoop
    	for (int ibuck = 0; ibuck < no_of_buckets; ibuck++) {
    		for (int ielm = 0; ielm < bucket[ibuck].ndx.size(); ielm++) {
    			Element* Curr_El = &(elenode_[bucket[ibuck].ndx[ielm]]);

    			if (Curr_El->adapted_flag() > 0 && Curr_El->state_vars(0)*scale_.height >= 0.05) {

    				Vx = Curr_El->state_vars(1)/Curr_El->state_vars(0);
    				Vy = Curr_El->state_vars(2)/Curr_El->state_vars(0);

    				Vm = sqrt(Vx*Vx+Vy*Vy);

    				T_Finertx += Curr_El->Drivs(0);
    				T_Finerty += Curr_El->Drivs(2);
    				T_Fconvx += Curr_El->Drivs(1);
    				T_Fconvy += Curr_El->Drivs(3);
    				T_Fgravx += Curr_El->STs(0);
    				T_Fgravy += Curr_El->STs(NUM_STATE_VARS);
    				T_Fbx += Curr_El->STs(1);
    				T_Fby += Curr_El->STs(1 + NUM_STATE_VARS);
    				T_Fix += Curr_El->STs(2);
    				T_Fiy += Curr_El->STs(2 + NUM_STATE_VARS);

    				T_Pinert += Curr_El->Drivs(0)*Vx + Curr_El->Drivs(2)*Vy;
    				T_Pconv += Curr_El->Drivs(1)*Vx + Curr_El->Drivs(3)*Vy;
    				T_Pgrav += Curr_El->STs(0)*Vx + Curr_El->STs(NUM_STATE_VARS)*Vy;
    				T_Pbed += Curr_El->STs(1)*Vx + Curr_El->STs(1 + NUM_STATE_VARS)*Vy;
    				T_Pint += Curr_El->STs(2)*Vx + Curr_El->STs(2 + NUM_STATE_VARS)*Vy;
    			}
    		}
    	}

    	T_Finertx *= timeprops.dtime;
    	T_Finerty *= timeprops.dtime;
    	T_Fconvx *= timeprops.dtime;
    	T_Fconvy *= timeprops.dtime;
    	T_Fgravx *= timeprops.dtime;
    	T_Fgravy *= timeprops.dtime;
    	T_Fbx *= timeprops.dtime;
    	T_Fby *= timeprops.dtime;
    	T_Fix *= timeprops.dtime;
    	T_Fiy *= timeprops.dtime;

    	T_Pinert *= timeprops.dtime;
    	T_Pconv *= timeprops.dtime;
    	T_Pgrav *= timeprops.dtime;
    	T_Pbed *= timeprops.dtime;
    	T_Pint *= timeprops.dtime;
    }

    double density = 1800.0;
    double F_SCALE = density * scale_.gravity * scale_.height * scale_.length * scale_.length;
    double P_SCALE = F_SCALE * sqrt(scale_.gravity * scale_.length);
    double A_SCALE = scale_.length * scale_.length;
    double T_SCALE= timeprops.TIME_SCALE;

    FILE *fp;
    fp = fopen("elements.dat","w");
    for(i = 0; i < N_elements.size(); i++)
    	fprintf(fp,"%.0f,%e\n", N_elements[i], Time_el[i]);
    fclose(fp);

    FILE *fp1;
    fp1 = fopen("Forces.dat","w");
    for(i = 0; i < Finertx.size(); i++)
    	fprintf(fp,"%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", Finertx[i]*F_SCALE, Finerty[i]*F_SCALE, Fconvx[i]*F_SCALE, Fconvy[i]*F_SCALE, Fgravx[i]*F_SCALE, Fgravy[i]*F_SCALE, Fbx[i]*F_SCALE, Fby[i]*F_SCALE, Fix[i]*F_SCALE, Fiy[i]*F_SCALE);
    fclose(fp1);

    FILE *fp2;
    fp2 = fopen("Powers.dat","w");
    for(i = 0; i < Pinert.size(); i++)
    	fprintf(fp,"%e,%e,%e,%e,%e\n", Pinert[i]*P_SCALE, Pconv[i]*P_SCALE, Pgrav[i]*P_SCALE, Pbed[i]*P_SCALE, Pint[i]*P_SCALE);
    fclose(fp2);

    FILE *fp3;
    fp3 = fopen("FdotV.dat","w");
    for(i = 0; i < FVinert.size(); i++)
    	fprintf(fp,"%e,%e,%e,%e,%e\n", FVinert[i]*F_SCALE, FVconv[i]*F_SCALE, FVgrav[i]*F_SCALE, FVbed[i]*F_SCALE, FVint[i]*F_SCALE);
    fclose(fp3);

//    FILE *fp4;
//    fp4 = fopen("Extents.dat","w");
//    for(i = 0; i < X_el.size(); i++)
//    	fprintf(fp,"%e,%e\n", X_el[i], Y_el[i]);
//    fclose(fp4);

    FILE *fp5;
    fp5 = fopen("ForceTimeInteg.dat","w");
    fprintf(fp,"%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", T_Finertx*F_SCALE*T_SCALE, T_Finerty*F_SCALE*T_SCALE, T_Fconvx*F_SCALE*T_SCALE, T_Fconvy*F_SCALE*T_SCALE, T_Fgravx*F_SCALE*T_SCALE, T_Fgravy*F_SCALE*T_SCALE, T_Fbx*F_SCALE*T_SCALE, T_Fby*F_SCALE*T_SCALE, T_Fix*F_SCALE*T_SCALE, T_Fiy*F_SCALE*T_SCALE);
    fclose(fp5);

    FILE *fp6;
    fp6 = fopen("PowerTimeInteg.dat","w");
    fprintf(fp,"%e,%e,%e,%e,%e\n", T_Pinert*P_SCALE*T_SCALE, T_Pconv*P_SCALE*T_SCALE, T_Pgrav*P_SCALE*T_SCALE, T_Pbed*P_SCALE*T_SCALE, T_Pint*P_SCALE*T_SCALE);
    fclose(fp6);

    FILE *fp7;
    fp7 = fopen("InundatedArea.dat","w");
    for(i = 0; i < Area.size(); i++)
    	fprintf(fp,"%e\n", Area[i]*A_SCALE);
    fclose(fp7);

    FILE *fp8;
    fp8 = fopen("MeanVel_TotalVolume.dat","w");
    for(i = 0; i < Volume.size(); i++)
        fprintf(fp,"%e,%e\n", Volume[i] * A_SCALE * scale_.height, MVelocity[i] * sqrt(scale_.gravity * scale_.length));
    fclose(fp8);


    move_data(numprocs, myid, ElemTable, NodeTable, &timeprops);
    PROFILING1_STOP(tsim_iter,pt_start0);

    /*
     * save a restart file
     */
    if(current_state_is_good&&restart_enabled)
    {
        save_restart_file();
    }
    IF_MPI(MPI_Barrier(MPI_COMM_WORLD));

    output_discharge(matprops_ptr, &timeprops, &discharge_planes, myid);
    IF_MPI(MPI_Barrier(MPI_COMM_WORLD));

    if(current_state_is_good)
        save_vizoutput_file(XDMF_CLOSE);
    else
        save_vizoutput_file(XDMF_ONLYCLOSE);

    IF_MPI(MPI_Barrier(MPI_COMM_WORLD));

    // write out ending warning, maybe flow hasn't finished moving
    sim_end_warning(elementType, ElemTable, matprops_ptr, &timeprops, statprops->vstar);
    IF_MPI(MPI_Barrier(MPI_COMM_WORLD));

    //write out the final pile statistics (and run time)
    if(myid == 0)
        out_final_stats(&timeprops, statprops);

    IF_MPI(MPI_Barrier(MPI_COMM_WORLD));

    //write out stochastic simulation statistics
    if(myid == 0)
        output_stoch_stats(matprops_ptr, statprops);
    IF_MPI(MPI_Barrier(MPI_COMM_WORLD));

    //output maximum flow depth a.k.a. flow outline
    PROFILING3_START(pt_start1);
    if(outline.enabled)outline.combine_results_from_threads();
    PROFILING3_STOPADD(step_outline,pt_start1);

    if(outline.enabled)
    {
#ifdef USE_MPI
        if(numprocs > 1)
        {
            OutLine outline2;
            double dxy[2];
            dxy[0] = outline.dx;
            dxy[1] = outline.dy;
            outline2.init2(dxy, outline.xminmax, outline.yminmax);
            int NxNyout = outline.Ny * outline.stride;

    //TWO PHASES is:MPI_Reduce(*(outline.pileheight), *(outline2.pileheight), NxNyout, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(outline.pileheight, outline2.pileheight, NxNyout, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Reduce(outline.max_kinergy, outline2.max_kinergy, NxNyout, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Reduce(outline.cum_kinergy, outline2.cum_kinergy, NxNyout, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            if(myid == 0)
                outline2.output(matprops_ptr, statprops);
        }
        else
        {
            outline.output(matprops_ptr, statprops);
        }
#else //USE_MPI
        outline.output(matprops_ptr, statprops);
#endif //USE_MPI
    }

#ifdef PERFTEST
    long m = element_counter, ii;
    MPI_Allreduce ( &element_counter, &ii, 1,
            MPI_LONG, MPI_SUM, MPI_COMM_WORLD );
    end=MPI_Wtime();
    char perffilename[256];
    sprintf(perffilename,"perform%04d.%04d",numprocs,myid);
    FILE *fpperf=fopen(perffilename,"w");
    fprintf(fpperf,"%d Finished -- used %ld elements of %ld total in %e seconds, %e\n",myid,m,ii,end-start, ii/(end-start));
    fclose(fpperf);
#endif
    TIMING1_STOP(totalTime, start);
    if(myid == 0)
        titanTimings.print();
    if(myid == 0)
        titanProfiling.print();

    IF_MPI(MPI_Barrier(MPI_COMM_WORLD));
    
    return;
}

void cxxTitanSimulation::input_summary()
{
    if(myid == 0)
    {
        get_mapnames()->print0();
        get_pileprops()->print0();
        get_fluxprops()->print0();
        get_discharge_planes()->print0();
        get_integrator()->print0();
        get_matprops()->print0();
    }
    IF_MPI(MPI_Barrier(MPI_COMM_WORLD));
    return;
}
