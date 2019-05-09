/*  $Id: seq_masker_util.hpp 103491 2007-05-04 17:18:18Z kazimird $
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
 *   Header file for CSeqMaskerUtil class.
 *
 */

#ifndef C_SEQ_MASKER_UTIL_H
#define C_SEQ_MASKER_UTIL_H

#include "mask_macros.h"

BEGIN_NCBI_SCOPE


/**
 **\brief Collection of various support utilities.
 **
 ** This class is used as a namespace for different utility 
 ** functions that can be used by other winmasker classes.
 **
 **/
class NCBI_XALGOWINMASK_EXPORT CSeqMaskerUtil
{
public:

    /**
     **\brief Count the bits with given value in a given bit pattern.
     **
     **\param mask the bit pattern
     **\param bit_value if 0 then 0-bits will be counted;
     **                 otherwise 1-bits will be counted
     **\return the bit count
     **
     **/
    static Uint1 BitCount( Uint4 mask, Uint1 bit_value = 1 );

    /**
     **\brief Reverse complement of a unit.
     **
     **\param seq the unit
     **\param size the unit length
     **\return the reverse complement of seq
     **
     **/
    static Uint4 reverse_complement( Uint4 seq, Uint1 size );

    /**
     **\brief Compute a hash code of a unit.
     **\param unit the target unit
     **\param k length (int bits) of the hash key
     **\param roff offset in bits from the right end of the unit
     **\return the hash code and concatenation of high and
     **         low remaining bits
     **/
    static pair< Uint4, Uint1 > hash_code( Uint4 unit, 
                                           Uint1 k, Uint1 roff )
    {
        return make_pair( (unit>>roff)&((((Uint4)1)<<k) - 1),
                          (Uint1)(((unit>>(roff + k))<<roff) 
                            + (unit&((1<<roff) - 1))) );
    }
};

END_NCBI_SCOPE

#endif
