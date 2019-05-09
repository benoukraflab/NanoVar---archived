/*  $Id: seq_masker_score_min.hpp 103491 2007-05-04 17:18:18Z kazimird $
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
 *   Header file for CSeqMaskerScoreMin class.
 *
 */

#ifndef C_SEQ_MASKER_SCORE_MIN_H
#define C_SEQ_MASKER_SCORE_MIN_H

//#include <corelib/ncbitype.h>

#include "seq_masker_score.hpp"
#include "seq_masker_istat.hpp"

BEGIN_NCBI_SCOPE

/**
 **\brief The score function object that computes maxmin of 
 **       k consecutive units in a window.
 **
 ** The score is computed as a maximum of k minimum units in
 ** the current window.
 **
 **/
class NCBI_XALGOWINMASK_EXPORT CSeqMaskerScoreMin : public CSeqMaskerScore
{
public:

    /**
     **\brief Object constructor.
     **
     **\param ustat unit length statitsics; used to initialize the
     **             base class instance
     **\param cnt the value of k
     **
     **/
    CSeqMaskerScoreMin( const CSeqMaskerIstat* ustat, Uint1 cnt = 0 )
        : CSeqMaskerScore( ustat ), count( cnt ) {}

    /**
     **\brief Object destructor.
     **/
    virtual ~CSeqMaskerScoreMin() {}

    /**
     **\brief Access the current score.
     **
     **\return the current score
     **
     **/
    virtual Uint4 operator()(); 

    /**\name Pre and postprocessing functions (trivial in this case).
     **/
    //@{
    virtual void PreAdvance( Uint4 step ) {}
    virtual void PostAdvance( Uint4 step ) {}
    //@}

protected:

    /**
     **\brief Object initialization.
     **
     ** This method is called automatically by SetWindow() method of
     ** the base class.
     **
     **/
    virtual void Init()
    {
        if( !count || count > window->NumUnits() )
            count = window->NumUnits();
    }

private:

    /**\internal
     **\brief The value of k (see the class description).
     **
     **/
    Uint1 count;
};

END_NCBI_SCOPE

#endif
