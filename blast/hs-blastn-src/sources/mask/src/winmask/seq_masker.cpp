/*  $Id: seq_masker.cpp 464810 2015-04-14 16:29:42Z vakatov $
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
 *   CSeqMasker class member and method definitions.
 *
 */

#include "ncbi_pch.hpp"

#include "ncbi_limits.h"

#include "seq_masker_window.hpp"
#include "seq_masker_window_ambig.hpp"
#include "seq_masker_window_pattern.hpp"
#include "seq_masker_window_pattern_ambig.hpp"
#include "seq_masker_score_mean.hpp"
#include "seq_masker_score_min.hpp"
#include "seq_masker_score_mean_glob.hpp"
#include "seq_masker.hpp"
#include "seq_masker_util.hpp"
#include "seq_masker_istat_factory.hpp"
#include "seq_masker_cache_boost.hpp"

#include <algorithm>
#include <memory>
#include <sstream>


BEGIN_NCBI_SCOPE
USING_SCOPE(objects);

#define WIN_MASK_ALGO_NAME "window-masker-algorithm"
#define WIN_MASK_ALGO_VER_MAJOR 1
#define WIN_MASK_ALGO_VER_MINOR 0
#define WIN_MASK_ALGO_VER_PATCH 0

//-------------------------------------------------------------------------
CSeqMaskerVersion CSeqMasker::AlgoVersion(
        WIN_MASK_ALGO_NAME,
        WIN_MASK_ALGO_VER_MAJOR,
        WIN_MASK_ALGO_VER_MINOR,
        WIN_MASK_ALGO_VER_PATCH
);

//-------------------------------------------------------------------------
CSeqMasker::CSeqMasker( const string & lstat_name,
                        Uint1 arg_window_size,
                        Uint4 arg_window_step,
                        Uint1 arg_unit_step,
                        Uint4 arg_textend,
                        Uint4 arg_cutoff_score,
                        Uint4 arg_max_score,
                        Uint4 arg_min_score,
                        Uint4 arg_set_max_score,
                        Uint4 arg_set_min_score,
                        bool arg_merge_pass,
                        Uint4 arg_merge_cutoff_score,
                        Uint4 arg_abs_merge_cutoff_dist,
                        Uint4 arg_mean_merge_cutoff_dist,
                        Uint1 arg_merge_unit_step,
                        const string & arg_trigger,
                        Uint1 tmin_count,
                        bool arg_discontig,
                        Uint4 arg_pattern,
                        bool arg_use_ba )
    : ustat( CSeqMaskerIstatFactory::create( lstat_name,
                                             arg_cutoff_score,
                                             arg_textend,
                                             arg_max_score,
                                             arg_set_max_score,
                                             arg_min_score,
                                             arg_set_min_score,
                                             arg_use_ba ) ),
      score( NULL ), score_p3( NULL ), trigger_score( NULL ),
      window_size( arg_window_size ), window_step( arg_window_step ),
      unit_step( arg_unit_step ),
      merge_pass( arg_merge_pass ),
      merge_cutoff_score( arg_merge_cutoff_score ),
      abs_merge_cutoff_dist( arg_abs_merge_cutoff_dist ),
      mean_merge_cutoff_dist( arg_mean_merge_cutoff_dist ),
      merge_unit_step( arg_merge_unit_step ),
      trigger( arg_trigger == "mean" ? eTrigger_Mean
               : eTrigger_Min ),
      discontig( arg_discontig ), pattern( arg_pattern )
{
    if( window_size == 0 ) window_size = ustat->UnitSize() + 4;

    if( window_size < ustat->UnitSize() ) {
        std::ostringstream os;
        os << "window size (" << (int)window_size << ") "
              "must be greater or equal to unit size (" <<
              (int)ustat->UnitSize() << ")";
        //NCBI_THROW( CSeqMaskerException, eValidation, os.str() );
		error_and_exit(os.str());
    }

    trigger_score = score = new CSeqMaskerScoreMean( ustat );

    if( trigger == eTrigger_Min )
        trigger_score = new CSeqMaskerScoreMin( ustat, tmin_count );

    if( !score )
    {
        //NCBI_THROW( CSeqMaskerException,
        //            eScoreAllocFail,
        //            "" );
		
		std::ostringstream os;
		os << "score allocation fail.";
		error_and_exit(os.str());
    }

    if( arg_merge_pass )
    {
        score_p3 = new CSeqMaskerScoreMeanGlob( ustat );

        if( !score )
        {
            //NCBI_THROW( CSeqMaskerException,
            //            eScoreP3AllocFail,
            //            "" );
			
			std::ostringstream os;
			os << "score_p3 allocation fail.";
			error_and_exit(os.str());
        }
    }
}

//-------------------------------------------------------------------------
CSeqMasker::~CSeqMasker()
{ 
    if( trigger_score != score ) delete trigger_score;

	delete ustat;
    delete score; 
    delete score_p3;
}

//-------------------------------------------------------------------------
CSeqMasker::TMaskList *
CSeqMasker::operator()( const CSeqVector& data ) const
{ return DoMask( data, 0, data.size() ); }

//-------------------------------------------------------------------------
void 
CSeqMasker::operator()(const CSeqVector& data, CSeqMasker::TMaskList& masked_locs) const
{
	DoMask(data, 0, data.size(), masked_locs);
}

//-------------------------------------------------------------------------
CSeqMasker::TMaskList *
CSeqMasker::DoMask( 
    const CSeqVector& data, TSeqPos begin, TSeqPos stop ) const
{
    ustat->total_ = 0;
    auto_ptr<TMaskList> mask(new TMaskList);
    Uint4 cutoff_score = ustat->get_threshold();
    Uint4 textend = ustat->get_textend();
    Uint1 nbits = discontig ? CSeqMaskerUtil::BitCount( pattern ) : 0;
    Uint4 unit_size = ustat->UnitSize() + nbits;
    auto_ptr<CSeqMaskerWindow> window_ptr
        (discontig ? new CSeqMaskerWindowPattern( data, unit_size, 
                                                  window_size, window_step, 
                                                  pattern, unit_step )
         : new CSeqMaskerWindow( data, unit_size, 
                                 window_size, window_step, 
                                 unit_step, begin, stop ));
    CSeqMaskerWindow & window = *window_ptr;
    score->SetWindow( window );

    if( trigger == eTrigger_Min ) trigger_score->SetWindow( window );

    Uint4 start = 0, end = 0, cend = 0;
    Uint4 limit = textend;
    const CSeqMaskerIstat::optimization_data * od 
        = ustat->get_optimization_data();

    CSeqMaskerCacheBoost booster( window, od );

    while( window )
    {
        Uint4 ts = (*trigger_score)();
        Uint4 s = (*score)();
        Uint4 adv = window_step;

        if( s < limit )
        {
            if( end > start )
            {
                if( window.Start() > cend )
                {
                    mask->push_back( TMaskedInterval( start, end ) );
                    start = end = cend = 0;
                }
            }

            if( od != 0 && od->cba_ != 0 )
            {
                adv = window.Start();

                if( !booster.Check() )
                    break;

                adv = window_step*( 1 + window.Start() - adv );
            }
        }
        else if( ts < cutoff_score )
        {
            if( end  > start )
            {
                if( window.Start() > cend + 1 )
                {
                    mask->push_back( TMaskedInterval( start, end ) );
                    start = end = cend = 0;
                }
                else cend = window.End();
            }
        }
        else
        {
            if( end > start )
            {
                if( window.Start() > cend + 1 )
                {
                    mask->push_back( TMaskedInterval( start, end ) );
                    start = window.Start();
                }
            }
            else start = window.Start();
    
            cend = end = window.End();
        }

        
        if( adv == window_step )
            ++window;

        score->PostAdvance( adv );
    }

    if( end > start ) 
        mask->push_back( TMaskedInterval( start, end ) );

    window_ptr.reset();

    if( merge_pass )
    {
        Uint1 nbits = discontig ? CSeqMaskerUtil::BitCount( pattern ) : 0;
        Uint4 unit_size = ustat->UnitSize() + nbits;

        if( mask->size() < 2 ) return mask.release();

        TMList masked, unmasked;
        TMaskList::iterator jtmp = mask->end();

        {{
             for( TMaskList::iterator i = mask->begin(), j = --jtmp; 
                  i != j; )
             {
                 masked.push_back( mitem( i->first, i->second, unit_size, 
                                          data, *this ) );
                 Uint4 nstart = (i++)->second - unit_size + 2;
                 unmasked.push_back( mitem( nstart, i->first + unit_size - 2, 
                                            unit_size, data, *this ) );
             }

             masked.push_back( mitem( (mask->rbegin())->first,
                                      (mask->rbegin())->second, 
                                      unit_size, data, *this ) );
         }}

        Int4 count = 0;
        TMList::iterator ii = masked.begin();
        TMList::iterator j = unmasked.begin();
        TMList::iterator k = ii, l = ii;
        --k; ++l;

        for( ; ii != masked.end(); k = l = ii, --k, ++l )
        {
            Uint4 ldist = (ii != masked.begin())
                ? ii->start - k->end - 1 : 0;
            TMList::iterator tmpend = masked.end();
            --tmpend;
            Uint4 rdist = (ii != tmpend)
                ? l->start - ii->end - 1 : 0;
            double lavg = 0.0, ravg = 0.0;
            bool can_go_left =  count && ldist
                && ldist <= mean_merge_cutoff_dist;
            bool can_go_right =  rdist
                && rdist <= mean_merge_cutoff_dist;

            if( can_go_left )
            {
                TMList::iterator tmp = j; --tmp;
                lavg = MergeAvg( k, tmp, unit_size );
                can_go_left = can_go_left && (lavg >= merge_cutoff_score);
            }

            if( can_go_right )
            {
                ravg = MergeAvg( ii, j, unit_size );
                can_go_right = can_go_right && (ravg >= merge_cutoff_score);
            }

            if( can_go_right )
            {
                if( can_go_left )
                {
                    if( ravg >= lavg )
                    {
                        ++count;
                        ++ii;
                        ++j;
                    }
                    else // count must be greater than 0.
                    {
                        --count;
                        k->avg = MergeAvg( k, --j, unit_size );
                        _TRACE( "Merging " 
                                << k->start << " - " << k->end
                                << " and " 
                                << ii->start << " - " << ii->end );
                        Merge( masked, k, unmasked, j );

                        if( count )
                        {
                            ii = --k;
                            --j;
                            --count;
                        }
                        else ii = k;
                    }
                }
                else
                {
                    ++count;
                    ++ii;
                    ++j;
                }
            }
            else if( can_go_left )
            {
                --count;
                k->avg = MergeAvg( k, --j, unit_size );
                _TRACE( "Merging " 
                        << k->start << " - " << k->end
                        << " and " 
                        << ii->start << " - " << ii->end );
                Merge( masked, k, unmasked, j );

                if( count )
                {
                    ii = --k;
                    --j;
                    --count;
                }
                else ii = k;
            }
            else
            {
                ++ii;
                ++j;
                count = 0;
            }
        }

        for( ii = masked.begin(), j = unmasked.begin(), k = ii++; 
             ii != masked.end(); (k = ii++), j++ )
        {
            if( k->end + abs_merge_cutoff_dist >= ii->start )
            {
                _TRACE( "Unconditionally merging " 
                        << k->start << " - " << k->end
                        << " and " 
                        << ii->start << " - " << ii->end );
                k->avg = MergeAvg( k, j, unit_size );
                Merge( masked, k, unmasked, j );
                ii = k; 

                if( ++ii == masked.end() ) break;
            }
        }

        mask->clear();

        for( TMList::const_iterator iii = masked.begin(); iii != masked.end(); ++iii )
            mask->push_back( TMaskedInterval( iii->start, iii->end ) );
    }

    return mask.release();
}

//-------------------------------------------------------------------------
void CSeqMasker::DoMask( 
    const CSeqVector& data, TSeqPos begin, TSeqPos stop, CSeqMasker::TMaskList& masked_locs ) const
{
    ustat->total_ = 0;
    //auto_ptr<TMaskList> mask(new TMaskList);
	TMaskList* mask = &masked_locs;

    Uint4 cutoff_score = ustat->get_threshold();
    Uint4 textend = ustat->get_textend();
    Uint1 nbits = discontig ? CSeqMaskerUtil::BitCount( pattern ) : 0;
    Uint4 unit_size = ustat->UnitSize() + nbits;
    auto_ptr<CSeqMaskerWindow> window_ptr
        (discontig ? new CSeqMaskerWindowPattern( data, unit_size, 
                                                  window_size, window_step, 
                                                  pattern, unit_step )
         : new CSeqMaskerWindow( data, unit_size, 
                                 window_size, window_step, 
                                 unit_step, begin, stop ));
    CSeqMaskerWindow & window = *window_ptr;
    score->SetWindow( window );

    if( trigger == eTrigger_Min ) trigger_score->SetWindow( window );

    Uint4 start = 0, end = 0, cend = 0;
    Uint4 limit = textend;
    const CSeqMaskerIstat::optimization_data * od 
        = ustat->get_optimization_data();

    CSeqMaskerCacheBoost booster( window, od );

    while( window )
    {
        Uint4 ts = (*trigger_score)();
        Uint4 s = (*score)();
        Uint4 adv = window_step;

        if( s < limit )
        {
            if( end > start )
            {
                if( window.Start() > cend )
                {
                    mask->push_back( TMaskedInterval( start, end ) );
                    start = end = cend = 0;
                }
            }

            if( od != 0 && od->cba_ != 0 )
            {
                adv = window.Start();

                if( !booster.Check() )
                    break;

                adv = window_step*( 1 + window.Start() - adv );
            }
        }
        else if( ts < cutoff_score )
        {
            if( end  > start )
            {
                if( window.Start() > cend + 1 )
                {
                    mask->push_back( TMaskedInterval( start, end ) );
                    start = end = cend = 0;
                }
                else cend = window.End();
            }
        }
        else
        {
            if( end > start )
            {
                if( window.Start() > cend + 1 )
                {
                    mask->push_back( TMaskedInterval( start, end ) );
                    start = window.Start();
                }
            }
            else start = window.Start();
    
            cend = end = window.End();
        }

        
        if( adv == window_step )
            ++window;

        score->PostAdvance( adv );
    }

    if( end > start ) 
        mask->push_back( TMaskedInterval( start, end ) );

    window_ptr.reset();

    if( merge_pass )
    {
        Uint1 nbits = discontig ? CSeqMaskerUtil::BitCount( pattern ) : 0;
        Uint4 unit_size = ustat->UnitSize() + nbits;

        if( mask->size() < 2 ) //return mask.release();
			return;

        TMList masked, unmasked;
        TMaskList::iterator jtmp = mask->end();

        {{
             for( TMaskList::iterator i = mask->begin(), j = --jtmp; 
                  i != j; )
             {
                 masked.push_back( mitem( i->first, i->second, unit_size, 
                                          data, *this ) );
                 Uint4 nstart = (i++)->second - unit_size + 2;
                 unmasked.push_back( mitem( nstart, i->first + unit_size - 2, 
                                            unit_size, data, *this ) );
             }

             masked.push_back( mitem( (mask->rbegin())->first,
                                      (mask->rbegin())->second, 
                                      unit_size, data, *this ) );
         }}

        Int4 count = 0;
        TMList::iterator ii = masked.begin();
        TMList::iterator j = unmasked.begin();
        TMList::iterator k = ii, l = ii;
        --k; ++l;

        for( ; ii != masked.end(); k = l = ii, --k, ++l )
        {
            Uint4 ldist = (ii != masked.begin())
                ? ii->start - k->end - 1 : 0;
            TMList::iterator tmpend = masked.end();
            --tmpend;
            Uint4 rdist = (ii != tmpend)
                ? l->start - ii->end - 1 : 0;
            double lavg = 0.0, ravg = 0.0;
            bool can_go_left =  count && ldist
                && ldist <= mean_merge_cutoff_dist;
            bool can_go_right =  rdist
                && rdist <= mean_merge_cutoff_dist;

            if( can_go_left )
            {
                TMList::iterator tmp = j; --tmp;
                lavg = MergeAvg( k, tmp, unit_size );
                can_go_left = can_go_left && (lavg >= merge_cutoff_score);
            }

            if( can_go_right )
            {
                ravg = MergeAvg( ii, j, unit_size );
                can_go_right = can_go_right && (ravg >= merge_cutoff_score);
            }

            if( can_go_right )
            {
                if( can_go_left )
                {
                    if( ravg >= lavg )
                    {
                        ++count;
                        ++ii;
                        ++j;
                    }
                    else // count must be greater than 0.
                    {
                        --count;
                        k->avg = MergeAvg( k, --j, unit_size );
                        _TRACE( "Merging " 
                                << k->start << " - " << k->end
                                << " and " 
                                << ii->start << " - " << ii->end );
                        Merge( masked, k, unmasked, j );

                        if( count )
                        {
                            ii = --k;
                            --j;
                            --count;
                        }
                        else ii = k;
                    }
                }
                else
                {
                    ++count;
                    ++ii;
                    ++j;
                }
            }
            else if( can_go_left )
            {
                --count;
                k->avg = MergeAvg( k, --j, unit_size );
                _TRACE( "Merging " 
                        << k->start << " - " << k->end
                        << " and " 
                        << ii->start << " - " << ii->end );
                Merge( masked, k, unmasked, j );

                if( count )
                {
                    ii = --k;
                    --j;
                    --count;
                }
                else ii = k;
            }
            else
            {
                ++ii;
                ++j;
                count = 0;
            }
        }

        for( ii = masked.begin(), j = unmasked.begin(), k = ii++; 
             ii != masked.end(); (k = ii++), j++ )
        {
            if( k->end + abs_merge_cutoff_dist >= ii->start )
            {
                _TRACE( "Unconditionally merging " 
                        << k->start << " - " << k->end
                        << " and " 
                        << ii->start << " - " << ii->end );
                k->avg = MergeAvg( k, j, unit_size );
                Merge( masked, k, unmasked, j );
                ii = k; 

                if( ++ii == masked.end() ) break;
            }
        }

        mask->clear();

        for( TMList::const_iterator iii = masked.begin(); iii != masked.end(); ++iii )
            mask->push_back( TMaskedInterval( iii->start, iii->end ) );
    }
}


//-------------------------------------------------------------------------
double CSeqMasker::MergeAvg( TMList::iterator mi, 
                             const TMList::iterator & umi,
                             Uint4 unit_size ) const
{
    TMList::iterator tmp = mi++;
    Uint4 n1 = (tmp->end - tmp->start - unit_size + 2)/merge_unit_step;
    Uint4 n2 = (umi->end - umi->start - unit_size + 2)/merge_unit_step;
    Uint4 n3 = (mi->end - mi->start - unit_size + 2)/merge_unit_step;
    Uint4 N = (mi->end - tmp->start - unit_size + 2)/merge_unit_step;
    double a1 = tmp->avg, a2 = umi->avg, a3 = mi->avg;
    return (a1*n1 + a2*n2 + a3*n3)/N;
}

//-------------------------------------------------------------------------
void CSeqMasker::Merge( TMList & m, TMList::iterator mi, 
                        TMList & um, TMList::iterator & umi ) const
{
    TMList::iterator tmp = mi++;
    tmp->end = mi->end;
    m.erase( mi );
    umi = um.erase( umi );
}

//----------------------------------------------------------------------------
//const char * CSeqMasker::CSeqMaskerException::GetErrCodeString() const
//{
//    switch( GetErrCode() )
//    {
//    case eLstatStreamIpenFail:

//        return "can not open input stream";

//    case eLstatSyntax:

//        return "syntax error";

//    case eLstatParam:

//        return  "the following parameters could not be determined"
//                " from the unit frequency database or command line: ";

//    case eScoreAllocFail:

//        return "score function object allocation failed";

//    case eScoreP3AllocFail:

//        return "merge pass score function object allocation failed";

//    case eValidation:

//        return "validation error";

//    default: 

//        return CException::GetErrCodeString();
//    }
//}

//----------------------------------------------------------------------------
CSeqMasker::mitem::mitem( Uint4 arg_start, Uint4 arg_end, Uint1 unit_size,
                          const CSeqVector & data, const CSeqMasker & owner )
    : start( arg_start ), end( arg_end ), avg( 0.0 )
{
    const Uint1 & window_size = owner.window_size;
    const CSeqMaskerWindow::TUnit & ambig_unit = owner.ustat->AmbigUnit();
    CSeqMaskerScore * const score = owner.score_p3;
    CSeqMaskerWindow * window = NULL;

    if( owner.discontig )
        window = new CSeqMaskerWindowPatternAmbig( data, unit_size, window_size, 
                                                   owner.merge_unit_step,
                                                   owner.pattern, ambig_unit, 
                                                   start, 
                                                   owner.merge_unit_step );
    else
        window = new CSeqMaskerWindowAmbig( data, unit_size, window_size, 
                                            owner.merge_unit_step, 
                                            ambig_unit, start, 
                                            owner.merge_unit_step );

    score->SetWindow( *window );
    Uint4 step = window->Step();

    while( window->End() < end )
    {
        score->PreAdvance( step );
        ++*window;
        score->PostAdvance( step );
    }

    avg = (*score)();
    delete window;
}

//----------------------------------------------------------------------------
void CSeqMasker::MergeMaskInfo( TMaskList * dest, const TMaskList * src )
{
    if( src->empty() )
        return;

    TMaskList::const_iterator si( src->begin() );
    TMaskList::const_iterator send( src->end() );
    TMaskList::iterator di( dest->begin() );
    TMaskList::iterator dend( dest->end() );
    TMaskList res;
    TMaskedInterval seg;
    TMaskedInterval next_seg;

    if( di != dend && di->first < si->first )
        seg = *(di++);
    else seg = *(si++);

    while( true )
    {
        if( si != send ) {
            if( di != dend ) {
                if( si->first < di->first ) {
                    next_seg = *(si++);
                } else {
                    next_seg = *(di++);
                }
            } else {
                next_seg = *(si++);
            }
        } else if( di != dend ) {
            next_seg = *(di++);
        } else {
            break;
        }

        if( seg.second + 1 < next_seg.first ) {
            res.push_back( seg );
            seg = next_seg;
        } 
        else if( seg.second < next_seg.second ) {
            seg.second = next_seg.second;
        }
    }

    res.push_back( seg );
    dest->swap(res);
}

CSeqMasker* s_BuildSeqMasker(const string& lstat)
{
    Uint1 arg_window_size            = 0; // [allow setting of this field?]
    Uint4 arg_window_step            = 1;
    Uint1 arg_unit_step              = 1;
    Uint4 arg_textend                = 0; // [allow setting of this field?]
    Uint4 arg_cutoff_score           = 0; // [allow setting of this field?]
    Uint4 arg_max_score              = 0; // [allow setting of this field?]
    Uint4 arg_min_score              = 0; // [allow setting of this field?]
    Uint4 arg_set_max_score          = 0; // [allow setting of this field?]
    Uint4 arg_set_min_score          = 0; // [allow setting of this field?]
    bool  arg_merge_pass             = false;
    Uint4 arg_merge_cutoff_score     = 0;
    Uint4 arg_abs_merge_cutoff_dist  = 0;
    Uint4 arg_mean_merge_cutoff_dist = 0;
    Uint1 arg_merge_unit_step        = 0;
    const string & arg_trigger       = "mean";
    Uint1 tmin_count                 = 0;
    bool  arg_discontig              = false;
    Uint4 arg_pattern                = 0;
	
	// enable/disable some kind of optimization
	bool arg_use_ba                  = true;
	
	// get a sequence masker
	CSeqMasker* masker = NULL;
	masker = new CSeqMasker( lstat,
							 arg_window_size,
							 arg_window_step,
							 arg_unit_step,
							 arg_textend,
							 arg_cutoff_score,
							 arg_max_score,
							 arg_min_score,
							 arg_set_max_score,
							 arg_set_min_score,
							 arg_merge_pass,
							 arg_merge_cutoff_score,
							 arg_abs_merge_cutoff_dist,
							 arg_mean_merge_cutoff_dist,
							 arg_merge_unit_step,
							 arg_trigger,
							 tmin_count,
							 arg_discontig,
							 arg_pattern,
							 arg_use_ba );
	
	assert(masker != NULL);
	
	return masker;
}

END_NCBI_SCOPE
