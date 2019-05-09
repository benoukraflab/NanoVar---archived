/*  $Id: seq_masker_ostat.hpp 462550 2015-03-19 14:07:19Z morgulis $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * Author:  Aleksandr Morgulis
 *
 * File Description:
 *   Definition of CSeqMaskerOstat class.
 *
 */

#ifndef C_WIN_MASK_USTAT_H
#define C_WIN_MASK_USTAT_H

#include <string>

//#include <corelib/ncbistre.hpp>
//#include <corelib/ncbistr.hpp>
//#include <corelib/ncbiargs.hpp>

#include "seq_masker_version.hpp"

BEGIN_NCBI_SCOPE

/**
 **\brief Base class for computing and saving unit counts data.
 **/
class NCBI_XALGOWINMASK_EXPORT CSeqMaskerOstat : public CObject
{
public:

    static char const * STAT_ALGO_COMPONENT_NAME;

    /// Version of the statistics generation algorithm.
    static CSeqMaskerVersion StatAlgoVersion;

    /**
        **\brief Exceptions that CSeqMaskerOstat can throw.
        **/
    //class CSeqMaskerOstatException : public CException
    //{
    //    public:

    //        enum EErrCode
    //        {
    //            eBadState   /**< Operation can not be performed in the current state. */
    //        };

    //        /**
    //            **\brief Get a description string for this exception.
    //            **\return C-style description string
    //            **/
    //        virtual const char * GetErrCodeString() const;

    //        NCBI_EXCEPTION_DEFAULT( CSeqMaskerOstatException, CException );
    //};

    /**
        **\brief Object constructor.
        **\param os C++ stream that should be used to save the unit counts data
        **/
    explicit CSeqMaskerOstat( 
            CNcbiOstream & os, bool alloc, string const & metadata );

    /**
        **\brief Trivial object destructor.
        **/
    virtual ~CSeqMaskerOstat() { if( alloc ) delete &out_stream; }

    /**
        **\brief Set the unit size value.
        **
        ** This method must be called before any call to setUnitCount().
        **
        **\param us new value of unit size 
        **/
    void setUnitSize( Uint1 us );

    /**
        **\brief Add count value for a particular unit.
        **
        ** This method can not be called before setUnitSize() and after
        ** the first call to setParam().
        **
        **\param unit unit value
        **\param count number of times the unit and its reverse complement
        **             occur in the genome
        **/
    void setUnitCount( Uint4 unit, Uint4 count );

    /**
        **\brief Add a comment to the unit counts file.
        **
        ** It is possible that this method is NOP for some unit counts
        ** formats.
        **
        **\param msg comment message
        **/
    void setComment( const string & msg ) { doSetComment( msg ); }

    /**
        **\brief Set a value of a WindowMasker parameter.
        **
        ** This method only can be called after all setUnitCount() calls.
        **
        **\param name the name of the parameter
        **\param value the value of the parameter
        **/
    void setParam( const string & name, Uint4 value );

    /**
        **\brief Perform any final tasks required to generate unit
        **       counts in the particular format.
        **/
    void finalize();

    /** Set the counts generation algorithm version explicitly
        (needed for convertions).
    */
    void SetStatAlgoVersion( CSeqMaskerVersion const & v ) {
        fmt_gen_algo_ver = v;
    }

    /** Get actual counts format version. */
    virtual CSeqMaskerVersion const & GetStatFmtVersion() const = 0;

protected:

    /** Algorithm parameter names. */
    static const char * PARAMS[];

    /**\name Methods used to delegate functionality to derived classes */
    /**@{*/
    virtual void doSetUnitSize( Uint4 us ) { unit_size = (Uint1)us; }
    virtual void doSetUnitCount( Uint4 unit, Uint4 count ) = 0;
    virtual void doSetComment( const string & msg ) {}
    virtual void doSetParam( const string & name, Uint4 value );
    // virtual void doSetBlank() {}
    virtual void doFinalize() {}
    /**@}*/

    /// Format algorithm parameters into a string.
    string FormatParameters() const;

    /// Combine version data and metadata into a single string.
    string FormatMetaData() const;

    /// Write metadata in binary format.
    void WriteBinMetaData( std::ostream & os ) const;

    /**
        **\brief Refers to the C++ stream that should be used to write
        **       out the unit counts data.
        **/
    CNcbiOstream& out_stream;

    /** flag indicating that the stream was allocated */
    bool alloc;

    /** metadata string */
    string metadata;

    /** unit size */
    Uint1 unit_size;

    /**\internal
        **\brief Parameter values.
        **/
    vector< Uint4 > pvalues;

    /** Unit counts */
    vector< pair< Uint4, Uint4 > > counts;

    /** version of the algorithm used to generate counts */
    CSeqMaskerVersion fmt_gen_algo_ver;

private:

    /**\name Provide reference semantics for the class */
    /**@{*/
    CSeqMaskerOstat( const CSeqMaskerOstat & );
    CSeqMaskerOstat& operator=( const CSeqMaskerOstat & );
    /**@}*/

    /**\internal
        **\brief Possible object states.
        **/
    enum 
    {
        start, /**<\internal The object has just been created. */
        ulen,  /**<\internal The unit size has been set. */
        udata, /**<\internal The unit counts data is being added. */
        thres, /**<\internal The parameters values are being set. */
        final  /**<\internal The finalization of the data structure. */
    } state;
};

END_NCBI_SCOPE

#endif
