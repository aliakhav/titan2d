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
 * $Id: hashtab.h 206 2009-01-26 17:32:10Z dkumar $ 
 */

/* Hash table */
/* Every table can process NBUCKETS entries */
/* All buckets with the same entry are linked together*/
/* The pointer to the first bucket is stored in the table */
/*---------------------------------------------------------*/

#ifndef HASHTABLE_H
#define HASHTABLE_H

#include <fstream>
#include <iostream>
#include <vector>
using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>


#include "constant.h"


#define HASHTABLE_LOOKUP_LINSEARCH 8

class Element;
class Node;
class BC;
class MatProps;


class HashEntryLine
{
public:
    vector<SFC_Key> key; //!<key: object key word
    vector<ti_ndx_t> ndx; //index of record
    
    //! find position of m_key in key
    ti_ndx_t lookup_local_ndx(const SFC_Key& m_key)
    {
        int size = key.size();

        SFC_Key *keyArr = &(key[0]);
        int i;

        if(size == 0)
            return ti_ndx_doesnt_exist;
        if(m_key < keyArr[0])
            return ti_ndx_doesnt_exist;
        if(m_key > keyArr[size - 1])
            return ti_ndx_doesnt_exist;

        if(size < HASHTABLE_LOOKUP_LINSEARCH)
        {
            for(i = 0; i < size; i++)
            {
                if(m_key == keyArr[i])
                {
                    return i;
                }
            }
        }
        else
        {
            int i0, i1, i2;
            i0 = 0;
            i1 = size / 2;
            i2 = size - 1;
            while ((i2 - i0) > HASHTABLE_LOOKUP_LINSEARCH)
            {
                if(m_key > keyArr[i1])
                {
                    i0 = i1 + 1;
                    i1 = (i0 + i2) / 2;
                }
                else
                {
                    i2 = i1;
                    i1 = (i0 + i2) / 2;
                }
            }
            for(i = i0; i <= i2; i++)
            {
                if(m_key == keyArr[i])
                {
                    return i;
                }
            }
        }
        return ti_ndx_doesnt_exist;
    }
    
};

enum ContentStatus {CS_Permanent=0,CS_Removed=-1,CS_Added=1};

template<typename T>
class HashTable
{
public:
    double doublekeyrange[2];
    double hashconstant;
    double Xrange[2];
    double Yrange[2];
    double invdxrange;
    double invdyrange;
    
    int NBUCKETS;
    int ENTRIES;
    
    vector<HashEntryLine> bucket;
    
    tivector<SFC_Key> key_;
    tivector<ContentStatus> status_;
    tivector<T> elenode_;
    
public:
    HashTable(double *doublekeyrangein, int size, double XR[], double YR[],tisize_t reserved_size);
    virtual ~HashTable(){}
    
    int hash(const SFC_Key& keyi) const
    {
        //Keith made this change 20061109; and made hash an inline function
        /* NBUCKETS*2 is NBUCKETS*integer integer is empirical could be 1
         return (((int) ((key[0]*doublekeyrange[1]+key[1])/
         (doublekeyrange[0]*doublekeyrange[1]+doublekeyrange[1])*
         NBUCKETS*2+0.5) )%NBUCKETS);
         */
    #if USE_ARRAY_SFC_KEY
        return (((int) ((key.key[0] * doublekeyrange[1] + key.key[1]) * hashconstant + 0.5)) % NBUCKETS);
    #else
        unsigned oldkey[KEYLENGTH];
        SET_OLDKEY(oldkey,keyi);
        return (((int) ((oldkey[0] * doublekeyrange[1] + oldkey[1]) * hashconstant + 0.5)) % NBUCKETS);
    #endif
    }
    
    ti_ndx_t lookup_ndx(const SFC_Key& keyi);
    T* lookup(const SFC_Key& keyi);
    
    ti_ndx_t add_ndx(const SFC_Key& keyi);
    T* add(const SFC_Key& keyi);
    
    void remove(const SFC_Key& keyi);
    
    void check();
    
    
    //plain getters and setters
    int get_no_of_buckets() const {return NBUCKETS;}
    double* get_Xrange() {return Xrange;}
    double* get_Yrange() {return Yrange;}
    double* get_doublekeyrange() {return doublekeyrange;}
    double get_invdxrange() const {return invdxrange;}
    double get_invdyrange() const {return invdyrange;}
    int get_no_of_entries() const {return ENTRIES;}

    void print0() const {}
    
    void print_out(int) const {}
    
    void flush();//actually delete, removed nodes and rearrange added (sort accorting to keys)

    tisize_t size() const{return elenode_.size();}
};
template class HashTable<Node>;
template class HashTable<Element>;

////////////////////////////////////////////////////////////////////////////////
//
constexpr tisize_t elem_reserved_size=1000000;
constexpr tisize_t node_reserved_size=elem_reserved_size*9;

class NodeHashTable: public HashTable<Node>
{
public:
    NodeHashTable(double *doublekeyrangein, int size, double XR[], double YR[]);
    virtual ~NodeHashTable();
    
    
    Node* createAddNode(const SFC_Key& keyi, double *coordi, MatProps *matprops_ptr);
    Node* createAddNode(const SFC_Key& keyi, double *coordi, int inf, int ord, MatProps *matprops_ptr);
    Node* createAddNode(const SFC_Key& keyi, double* coordi, int inf, int ord, double elev, int yada);
    Node* createAddNode(FILE* fp, MatProps* matprops_ptr);
    
    void removeNode(Node* node);

    Node& node(const ti_ndx_t ndx){return elenode_[ndx];}
};
extern NodeHashTable *nodeHashTable;

////////////////////////////////////////////////////////////////////////////////
//! Hashtables for Elements
class ElementsHashTable: public HashTable<Element>
{
    
    //friend int hash(unsigned* keyi);
    friend class Element;
protected:
    NodeHashTable* NodeTable;

    vector<SFC_Key> ukeyElements;
    vector<Element*> elements;

    int NlocalElements;
    vector<SFC_Key> ukeyLocalElements;
    vector<Element*> localElements;
public:
    ElementsHashTable(double *doublekeyrangein, int size, double XR[], double YR[], NodeHashTable* nodeTable);
    virtual ~ElementsHashTable();

    Element** getElements()
    {
        return &(elements[0]);
    }
    
    void updateElements();
    //!debug function check that all allEntries are up to date, return number of mismatch
    int ckeckElementsPointers(const char *prefix);

    int getNumberOfLocalElements()
    {
        return NlocalElements;
    }
    Element** getLocalElementsValues()
    {
        return &(localElements[0]);
    }
    //only works for elements
    void updateLocalElements();
    //!debug function check that all allEntries are up to date, return number of mismatch
    int ckeckLocalElementsPointers(const char *prefix);

    //! set neighboring elements and nodes pointers in elements
    void updatePointersToNeighbours();

    //! check that neighboring elements and nodes pointers are correct
    int checkPointersToNeighbours(const char *prefix);
        
    //! default Element generator constructor, does nothing except set stoppedflags=2, this should never be used
    virtual Element* generateAddElement(const SFC_Key& keyi);

    //! constructor that creates an original element when funky is read in
    virtual Element* generateAddElement(const SFC_Key* nodekeys, const SFC_Key* neigh, int n_pro[], BC *b, int mat, int *elm_loc_in,
                             double pile_height, int myid, const SFC_Key& opposite_brother);

    //! constructor that creates a son element from its father during refinement
    virtual Element* generateAddElement(const SFC_Key* nodekeys, const SFC_Key* neigh, int n_pro[], BC *b, int gen, int elm_loc_in[],
                int *ord, int gen_neigh[], int mat, Element *fthTemp, double *coord_in, ElementsHashTable *El_Table,
                NodeHashTable *NodeTable, int myid, MatProps *matprops_ptr, int iwetnodefather, double Awetfather,
                double *drypoint_in);
    //! constructor that creates a father element from its four sons during unrefinement
    virtual Element* generateAddElement(Element *sons[], NodeHashTable *NodeTable, ElementsHashTable *El_Table, MatProps *matprops_ptr);

    //! constructor that creates/restores a saved element during restart
    virtual Element* generateAddElement(FILE* fp, NodeHashTable* NodeTable, MatProps* matprops_ptr, int myid);
    
    //!should not be used
    void removeElement(Element* elm);
    

    Element& elem(const ti_ndx_t ndx){return elenode_[ndx];}
    //here goes element content storage probably should be separate class at the end
    
};

extern ElementsHashTable *elementsHashTable;

#endif
