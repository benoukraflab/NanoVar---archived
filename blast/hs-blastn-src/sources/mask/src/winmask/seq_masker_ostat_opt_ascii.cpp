/*  $Id: seq_masker_ostat_opt_ascii.cpp 462550 2015-03-19 14:07:19Z morgulis $
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
 *   Implementation of CSeqMaskerOStatOptAscii class.
 *
 */

#include "ncbi_pch.hpp"

#include "seq_masker_ostat_opt_ascii.hpp"

BEGIN_NCBI_SCOPE

#define STAT_FMT_COMPONENT_NAME "windowmasker-statistics-format-version"
#define STAT_FMT_VER_MAJOR 1
#define STAT_FMT_VER_MINOR 0
#define STAT_FMT_VER_PATCH 0
#define STAT_FMT_VER_PFX "oascii "

//------------------------------------------------------------------------------
CSeqMaskerVersion CSeqMaskerOstatOptAscii::FormatVersion(
        STAT_FMT_COMPONENT_NAME, 
        STAT_FMT_VER_MAJOR,
        STAT_FMT_VER_MINOR,
        STAT_FMT_VER_PATCH,
        STAT_FMT_VER_PFX
);

//------------------------------------------------------------------------------
CSeqMaskerOstatOptAscii::CSeqMaskerOstatOptAscii( 
        const string & name, Uint2 sz, string const & metadata )
    : CSeqMaskerOstatOpt( static_cast< CNcbiOstream& >(
        *new CNcbiOfstream( name.c_str() ) ), sz, true, metadata ) 
{ 
}

//------------------------------------------------------------------------------
CSeqMaskerOstatOptAscii::CSeqMaskerOstatOptAscii( 
        CNcbiOstream & os, Uint2 sz, string const & metadata )
    : CSeqMaskerOstatOpt( os, sz, false, metadata ) 
{ 
}

//------------------------------------------------------------------------------
void CSeqMaskerOstatOptAscii::write_out( const params & p ) const
{
    out_stream << FormatMetaData();

    // File format identifier
    out_stream << (char)65; 
    out_stream << (char)65; 
    out_stream << (char)65; 
    out_stream << (char)65 << endl; 

    out_stream << (Uint4)UnitSize() << "\n";
    out_stream << p.M << " "
               << (Uint4)p.k << " " << (Uint4)p.roff << " " << (Uint4)p.bc
               << "\n";

    for( Uint4 i = 0; i < GetParams().size(); ++i )
        out_stream << GetParams()[i] << "\n";

    for( Uint4 i = 0, sz = (1<<p.k); i < sz; ++i )
        out_stream << p.ht[i] << "\n";

    for( Uint4 i = 0; i < p.M; ++i )
        out_stream << (Uint4)(p.vt[i]) << "\n";

    out_stream << flush;
}

END_NCBI_SCOPE
