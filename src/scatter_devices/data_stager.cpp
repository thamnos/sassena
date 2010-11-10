/*
 *  scatterdevices.cpp
 *
 *  Created on: Dec 30, 2008
 *  Authors:
 *  Benjamin Lindner, ben@benlabs.net
 *
 *  Copyright 2008,2009 Benjamin Lindner
 *
 */

// direct header
#include "scatter_devices/data_stager.hpp"

// standard header
#include <complex>
#include <fstream>

// special library headers
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/thread.hpp>

// other headers
#include "math/coor3d.hpp"
#include "decomposition/decompose.hpp"
#include <fftw3.h>
#include "control.hpp"
#include "log.hpp"
#include "sample.hpp"

using namespace std;

struct binary_max : public binary_function<size_t,size_t,size_t> {
    size_t operator() (size_t a,size_t b) {if (a>b) return a; else return b;}
};

DataStagerByFrame::DataStagerByFrame(Sample& sample,boost::mpi::communicator& allcomm,boost::mpi::communicator& partitioncomm,DivAssignment assignment) 
    : m_sample(sample),
    allcomm_(allcomm),
    partitioncomm_(partitioncomm),    
    FC_assignment(assignment)
{
    NN = allcomm_.size();    
    NF = m_sample.coordinate_sets.size();
    std::string target = Params::Inst()->scattering.target;
    NA = m_sample.atoms.selections[target]->size();

    size_t rank = allcomm_.rank();

    // check data size requirements & allocate memory
    size_t data_bytesize = FC_assignment.size()*NA*3*sizeof(coor_t);
    size_t data_bytesize_indicator_max = 0;
    boost::mpi::all_reduce(allcomm,data_bytesize,data_bytesize_indicator_max,binary_max());

    if (Params::Inst()->limits.memory.data<data_bytesize_indicator_max) {
        if (rank==0) {
            Err::Inst()->write("Insufficient Buffer size for coordinates (limits.memory.data)");
            Err::Inst()->write(string("Requested (bytes): ")+boost::lexical_cast<string>(data_bytesize_indicator_max));
        }
        throw;
    }

    p_coordinates = (coor_t*) malloc(data_bytesize);
        
    // determine number of nodes which act as file servers:
    NFN = Params::Inst()->limits.data_stager.servers;
    
    if (NFN>NN) NFN=NN; 
    if (NFN>NF) {
        NFN=NF;  
        if (allcomm_.rank()==0) {
            Info::Inst()->write("Number of data servers limited by number of frames");
        }
    }    
}

coor_t* DataStagerByFrame::stage() {
    stage_firstpartition();
    allcomm_.barrier();
    stage_fillpartitions();
    return p_coordinates;
}

void DataStagerByFrame::stage_firstpartition() {
    size_t rank = allcomm_.rank();
    
    // used by server
    //size_t frame_bytesize = NA*3*sizeof(double);
    coor_t* p_coordinates_buffer = (coor_t*) malloc(NA*3*sizeof(coor_t));

    //std::set<size_t>& assigned_frames = FS_assignment_table[rank];

    vector<DivAssignment> assignments;
    for (size_t i=0;i<NN;i++) {
        assignments.push_back(DivAssignment(NN,i,NF));
    }
    
    bool firstpartition=true;
    if (allcomm_.rank()>=partitioncomm_.size()) firstpartition=false;

    for(size_t f=0;f<NF;f++) {
        size_t s = f%NFN; // this is the responsible data server
        if (rank==s) {
            CoordinateSet* p_cset = m_sample.coordinate_sets.load(f);
            
            for(size_t n=0;n<NA;n++) {
                p_coordinates_buffer[3*n]=p_cset->c1[n];
                p_coordinates_buffer[3*n+1]=p_cset->c2[n];
                p_coordinates_buffer[3*n+2]=p_cset->c3[n];            
            }
            delete p_cset;

            size_t target_node = 0;
            for(size_t i = 0; i < assignments.size(); ++i)
            {
                if (assignments[i].contains(f)) {
                    target_node=i;
                    break;
                }
            }
            if (target_node==s) {
                coor_t* p_localdata = &(p_coordinates[FC_assignment.index(f)*NA*3]);
                memcpy(p_localdata,p_coordinates_buffer,3*NA*sizeof(coor_t));                
            } else {
                allcomm_.send(target_node,0,p_coordinates_buffer,3*NA);
            }
            
        } else if (firstpartition && FC_assignment.contains(f)) {
            coor_t* p_localdata = &(p_coordinates[FC_assignment.index(f)*NA*3]);
            allcomm_.recv(s,0,p_localdata,3*NA);            
        }
    }
    free(p_coordinates_buffer);
}


void DataStagerByFrame::stage_fillpartitions() {
    
    bool firstpartition=true;
    if (allcomm_.rank()>=partitioncomm_.size()) firstpartition=false;
    size_t partitions = allcomm_.size()/partitioncomm_.size();

    if (firstpartition) {
        for(size_t i = 1; i < partitions; ++i)
        {
            size_t target_node = i*partitioncomm_.size() + partitioncomm_.rank();
            allcomm_.send(target_node,0,p_coordinates,FC_assignment.size()*NA*3);            
        }
    } else {
        size_t source_node = partitioncomm_.rank();
        allcomm_.recv(source_node,0,p_coordinates,FC_assignment.size()*NA*3);            
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// By Atom
////////////////////////////////////////////////////////////////////////////////////////////////////


DataStagerByAtom::DataStagerByAtom(Sample& sample,boost::mpi::communicator& comm, DivAssignment assignment) 
    : m_sample(sample),
    m_comm(comm),
    FC_assignment(assignment)
{
    NN = m_comm.size();    
    NF = m_sample.coordinate_sets.size();
    std::string target = Params::Inst()->scattering.target;
    NA = m_sample.atoms.selections[target]->size();

    size_t rank = m_comm.rank();

    // check data size requirements & allocate memory
    size_t data_bytesize = FC_assignment.size()*NF*3*sizeof(coor_t);
    size_t data_bytesize_indicator_max = 0;
    boost::mpi::all_reduce(m_comm,data_bytesize,data_bytesize_indicator_max,binary_max());

    if (Params::Inst()->limits.memory.data<data_bytesize_indicator_max) {
        if (rank==0) {
            Err::Inst()->write("Insufficient Buffer size for coordinates (limits.memory.data)");
            Err::Inst()->write(string("Requested (bytes): ")+boost::lexical_cast<string>(data_bytesize_indicator_max));
        }
        throw;
    }

//    Info::Inst()->write(string("Maximum memory allocated for coordinates (bytes): ")+boost::lexical_cast<string>(data_bytesize_indicator_max));
    p_coordinates = (coor_t*) malloc(data_bytesize);
        
    // determine number of nodes which act as file servers:
    NFN = Params::Inst()->limits.data_stager.servers;
    
    if (NFN>NN) NFN=NN; 
    if (NFN>NF) {
        NFN=NF;  
        if (m_comm.rank()==0) {
            Info::Inst()->write("Number of data servers limited by number of frames");
        }
    }    

    // file FS_assignment_table, all client own this list
    for(size_t i=0;i<NF;i++) {
        size_t assigned_server = i%NFN;
        FS_assignment_table[assigned_server].insert(i);
    }

}

coor_t* DataStagerByAtom::stage() {
    stage_registration();    
    stage_data();
    return p_coordinates;
}

void DataStagerByAtom::stage_registration() {
    size_t off = FC_assignment[0];
    size_t len = FC_assignment.size();
    
    size_t NN = m_comm.size();
    size_t rank = m_comm.rank();
    
    for(size_t s = 0; s < NFN; ++s)
    {
        for(size_t c = 0; c < NN; ++c)
        {
            if (rank==s) {
                if (s==c) {            
                    add_client(off,len,c);                    
                } else {
                    size_t coff,clen;
                    m_comm.recv(c,0,&coff,1);
                    m_comm.recv(c,0,&clen,1);                    
                    add_client(coff,clen,c);
                }
            }
            if (rank==c) {
                m_comm.send(s,0,&off,1);
                m_comm.send(s,0,&len,1);                                    
            }
        }
    }
}

void DataStagerByAtom::add_client(size_t off,size_t len, size_t client) {
    bool found = false;
    for(size_t i = 0; i < AtomOffLen_to_FClist.size(); ++i)
    {
        if (AtomOffLen_to_FClist[i].first.first == off) {
            if (AtomOffLen_to_FClist[i].first.second == len) {
                AtomOffLen_to_FClist[i].second.push_back(client);
                found = true;
                break;
            }
        }
    }
    if (!found) {
        std::pair<size_t,size_t> key = make_pair(off,len);
        vector<size_t> val; val.push_back(client);
        AtomOffLen_to_FClist.push_back(make_pair(key,val));
    }
}

void DataStagerByAtom::stage_data() {

	size_t rank = m_comm.rank();

	//////////////////////////////////////////////////
	// Assignment of atoms to compute
	//////////////////////////////////////////////////    
    
    size_t buffer_bytesize = Params::Inst()->limits.memory.data_stager;
    size_t frame_bytesize = NA*3*sizeof(coor_t);
    size_t framesbuffer_maxsize = buffer_bytesize/frame_bytesize;
    
    if (framesbuffer_maxsize==0) {
        if (m_comm.rank()==0) {
            Err::Inst()->write("Cannot load trajectory into buffer.");
            Err::Inst()->write(string("limits.memory.data_stager=")+boost::lexical_cast<string>(Params::Inst()->limits.memory.data_stager));
            Err::Inst()->write(string("requested=")+boost::lexical_cast<string>(frame_bytesize));            
        }
        throw;
    }
    
    coor_t* p_coordinates_buffer = (coor_t*) malloc(framesbuffer_maxsize*NA*3*sizeof(coor_t));
    
    std::vector<std::vector<size_t> > framesbuffer(NFN);
    for(size_t f=0;f<NF;f++) {
        size_t s = f%NFN; // this is the responsible data server
        if (rank==s) {
            CoordinateSet* p_cset = m_sample.coordinate_sets.load(f);
            coor_t* p_data = &(p_coordinates_buffer[framesbuffer[s].size()*NA*3]);
            
            for(size_t n=0;n<NA;n++) {
                p_data[3*n]=p_cset->c1[n];
                p_data[3*n+1]=p_cset->c2[n];
                p_data[3*n+2]=p_cset->c3[n];            
            }
            delete p_cset;
        }
        framesbuffer[s].push_back(f);
        if (framesbuffer[s].size()==framesbuffer_maxsize) {
            distribute_coordinates(p_coordinates_buffer,framesbuffer,s);            
            framesbuffer[s].clear();
        }
    }
    
    for(size_t i = 0; i < framesbuffer.size(); ++i)
    {
        if (framesbuffer[i].size()!=0) {
            distribute_coordinates(p_coordinates_buffer,framesbuffer,i);            
            framesbuffer[i].clear();
        }
    }
    
    free(p_coordinates_buffer);
}

void DataStagerByAtom::distribute_coordinates(coor_t* p_coordinates_buffer,std::vector<std::vector<size_t> >& framesbuffer,size_t s) {
    
    size_t LNF = framesbuffer[s].size();
    size_t rank = m_comm.rank();  
    
    if (rank==s) {
        for(size_t i = 0; i < AtomOffLen_to_FClist.size(); ++i)
        {
            std::pair<size_t,size_t> offlen = AtomOffLen_to_FClist[i].first;
            std::vector<size_t> clients = AtomOffLen_to_FClist[i].second;

            size_t off = offlen.first;
            size_t len = offlen.second;
            // marshal data            
            coor_t* p_fsdata = (coor_t*) malloc(LNF*len*3*sizeof(coor_t));
            for(size_t f = 0; f < LNF; ++f)
            {
                coor_t* p_from = &(p_coordinates_buffer[f*NA*3+off*3]);
                coor_t* p_to = &(p_fsdata[f*len*3]);
                memcpy(p_to,p_from,len*3*sizeof(coor_t));
            }

            for(size_t j = 0; j < clients.size(); ++j)
            {
                size_t target_node = clients[j];
                if (target_node==s) {                    
                    fill_coordinates(p_fsdata,len,framesbuffer[s]);
                } else {
                    m_comm.send(target_node,0,p_fsdata,LNF*len*3);
                }
            }
            free(p_fsdata);
        }
    } else {
        size_t len = FC_assignment.size();
        
        coor_t* p_localdata = (coor_t*) malloc(LNF*len*3*sizeof(coor_t));
        m_comm.recv(s,0,p_localdata,LNF*len*3);
        fill_coordinates(p_localdata,len,framesbuffer[s]);
        free(p_localdata);
    }
}

void DataStagerByAtom::fill_coordinates(coor_t* p_localdata,size_t len,vector<size_t> frames) {
    for(size_t f = 0; f < frames.size(); ++f)
    {
        for(size_t n = 0; n < len; ++n)
        {
            p_coordinates[ NF*n*3 + frames[f]*3 ] = p_localdata[ len*f*3 + n*3 ];
            p_coordinates[ NF*n*3 + frames[f]*3 + 1 ] = p_localdata[ len*f*3 + n*3 + 1 ];
            p_coordinates[ NF*n*3 + frames[f]*3 + 2 ] = p_localdata[ len*f*3 + n*3 + 2 ];
        }
    }
}
