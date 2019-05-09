/*  $Id: seq_masker_score_mean.cpp 103491 2007-05-04 17:18:18Z kazimird $
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
 *   CSeqMaskerScoreMean class member and method definitions.
 *
 */

#include "ncbi_pch.hpp"
#include "seq_masker_window.hpp"
#include "seq_masker_score_mean.hpp"

BEGIN_NCBI_SCOPE


//-------------------------------------------------------------------------
CSeqMaskerScoreMean::CSeqMaskerScoreMean( 
    const CSeqMaskerIstat* ustat )
: CSeqMaskerScore( ustat ), sum( 0 ), start( 0 ), num( 0 )
{
}

//-------------------------------------------------------------------------
Uint4 CSeqMaskerScoreMean::operator()()
{ return sum/num; }

//-------------------------------------------------------------------------
void CSeqMaskerScoreMean::PreAdvance( Uint4 step )
{
    if( step == 1 && window->UnitStep() == 1 )
    {
        start = window->Start();
        sum -= *scores_start;
    }
}

//-------------------------------------------------------------------------
void CSeqMaskerScoreMean::PostAdvance( Uint4 step )
{
    if(    step == 1 
           && window->UnitStep() == 1 
           && window->Start() - start == 1 )
    {
        /*!!!!NEW CODE*/ sum -= *scores_start;
        *scores_start = (*ustat)[(*window)[num - 1]];
        sum += *scores_start;
        scores_start = (scores_start - &scores[0] == (int)(num - 1) ) 
	             ? &scores[0]
                     : scores_start + 1;
        /*!!!!NEW CODE*/ start = window->Start();
    }
    else{ FillScores(); }
}

//-------------------------------------------------------------------------
void CSeqMaskerScoreMean::Init()
{
    start = window->Start();
    num = window->NumUnits();
    scores.resize( num, 0 );
  
    FillScores();
}

//-------------------------------------------------------------------------
void CSeqMaskerScoreMean::FillScores()
{
  sum = 0;
  scores_start = &scores[0];

  for( Uint1 i = 0; i < num; ++i )
  {
    scores[i] = (*ustat)[(*window)[i]];
    sum += scores[i];
  }

  /*!!!!NEW CODE*/ start = window->Start();
}

END_NCBI_SCOPE
