/*  $Id: seq_masker_window_pattern_ambig.cpp 103491 2007-05-04 17:18:18Z kazimird $
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
 *   CSeqMaskerWindowPatternAmbig class member and method definitions.
 *
 */

#include "ncbi_pch.hpp"
#include <string>

#include "seq_masker_util.hpp"
#include "seq_masker_window_pattern_ambig.hpp"
//#include <objmgr/seq_vector.hpp>

BEGIN_NCBI_SCOPE
USING_SCOPE(objects);

//-------------------------------------------------------------------------
CSeqMaskerWindowPatternAmbig::
CSeqMaskerWindowPatternAmbig(const CSeqVector & arg_data,
                             Uint1 arg_unit_size, Uint1 arg_window_size,
                             Uint4 window_step, Uint4 arg_pattern,
                             TUnit arg_ambig_unit, Uint4 window_start,
                             Uint1 arg_unit_step )
    : CSeqMaskerWindowPattern( arg_data, arg_unit_size, arg_window_size,
                               window_step, arg_pattern, arg_unit_step ),
      ambig_unit( arg_ambig_unit ), ambig( false )
{
    FillWindow( window_start );
}


//-------------------------------------------------------------------------
void CSeqMaskerWindowPatternAmbig::Advance( Uint4 step )
{
    FillWindow( start + step );
}

//-------------------------------------------------------------------------
void CSeqMaskerWindowPatternAmbig::FillWindow( Uint4 winstart )
{
    first_unit = 0;
    TUnit unit = 0;
    Int4 iter = 0;
    end = winstart + unit_size - 1;

    for( ; iter < NumUnits() && end < data.size(); 
         ++iter, end += unit_step, winstart += unit_step ) {
        if( MakeUnit( winstart, unit ) ) {
            units[iter] = unit;
        } else {
            units[iter] = ambig_unit;
        }
    }

    end -= unit_step;
    end += (window_size - unit_size)%unit_step;
    start = end - window_size + 1;
    state = (iter == NumUnits());
}


END_NCBI_SCOPE
