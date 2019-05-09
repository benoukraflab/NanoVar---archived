/* $Id: greedy_align.h 392014 2013-03-13 14:36:38Z maning $
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
 * Author: Ilya Dondoshansky
 *
 */

/** @file greedy_align.h
 * Prototypes and structures for greedy gapped alignment
 */

#ifndef GREEDY_ALIGN_IMPL_H
#define	GREEDY_ALIGN_IMPL_H

#include "greedy_align.h"

#include <cmath>

/** Signal that a diagonal is invalid */
static const Int4 kInvalidOffset = -2;

static inline Int4
s_FindFirstMismatch(const Uint1* seq1, const Uint1* seq2,
                    Int4 len1, Int4 len2, Int4 seq1_index,
                    Int4 seq2_index, Boolean reverse, Uint1 rem)
{
    Int4 tmp = seq1_index;
	Uint1 c;

    if (reverse)
    {
        if (rem == 4)
        {
            while (seq1_index < len1 && seq2_index < len2 &&
                   seq1[len1 - 1 - seq1_index] < 4 &&
                   seq1[len1 - 1 - seq1_index] == seq2[len2 - 1 - seq2_index])
            {
                ++seq1_index;
                ++seq2_index;
            }
        }
        else
        {
            while (seq1_index < len1 && seq2_index < len2 &&
		   (c = Blastna2Na2(seq2[len2 - 1 - seq2_index])) < 4 &&
                   seq1[len1 - 1 - seq1_index] == c)
            {
                ++seq1_index;
                ++seq2_index;
            }
        }
    }
    else
    {
        if (rem == 4)
        {
            while (seq1_index < len1 && seq2_index < len2 &&
                   seq1[seq1_index] < 4 &&
                   seq1[seq1_index] == seq2[seq2_index])
            {
                ++seq1_index;
                ++seq2_index;
            }
        }
        else
        {
            while (seq1_index < len1 && seq2_index < len2 &&
		   (c = Blastna2Na2(seq2[seq2_index])) < 4 &&
                   seq1[seq1_index] == c)
            {
                ++seq1_index;
                ++seq2_index;
            }
        }
    }

    return seq1_index - tmp;
}

/** During the traceback for a non-affine greedy alignment,
    compute the diagonal that will result from the next
    traceback operation

    @param last_seq2_off Array of offsets into the second sequence;
                        last_seq2_off[d][k] gives the largest offset into
                        the second sequence that lies on diagonal k and
                        has distance d [in]
    @param d Starting distance [in]
    @param diag Index of diagonal that produced the starting distance [in]
    @param seq2_index The offset into the second sequence after the traceback
                operation has completed [out]
    @return The diagonal resulting from the next traceback operation
                being applied
*/
static inline Int4
s_GetNextNonAffineTback(Int4 **last_seq2_off, Int4 d,
                        Int4 diag, Int4 *seq2_index)
{
    /* choose the traceback operation that results in the
       largest seq2 offset at this point, then compute the
       new diagonal that is implied by the operation */

    if (last_seq2_off[d-1][diag-1] >
                MAX(last_seq2_off[d-1][diag], last_seq2_off[d-1][diag+1])) {
        *seq2_index = last_seq2_off[d-1][diag-1];
        return diag - 1;    /* gap in seq2 */
    }
    if (last_seq2_off[d-1][diag] > last_seq2_off[d-1][diag+1]) {
        *seq2_index = last_seq2_off[d-1][diag];
        return diag;        /* match */
    }
    *seq2_index = last_seq2_off[d-1][diag+1];
    return diag + 1;        /* gap in seq1 */
}

/** During the traceback for a greedy alignment with affine
    gap penalties, determine the next state of the traceback after
    moving upwards in the traceback array from a substitution

    @param last_seq2_off Array of offsets into the second sequence;
                        last_seq2_off[d][k] gives the largest offset into
                        the second sequence that lies on diagonal k and
                        has distance d [in]
    @param diag_lower Array of lower bounds on diagonal index [in]
    @param diag_upper Array of upper bounds on diagonal index [in]
    @param d Starting distance [in][out]
    @param diag The diagonal of the current traceback position [in]
    @param op_cost The sum of the match and mismatch scores [in]
    @param seq2_index The offset into the second sequence after the traceback
                operation has completed [out]
    @return The state for the next traceback operation
*/
static inline EGapAlignOpType
s_GetNextAffineTbackFromMatch(SGreedyOffset** last_seq2_off, Int4* diag_lower,
                           Int4* diag_upper, Int4* d, Int4 diag, Int4 op_cost,
                           Int4* seq2_index)
{
    Int4 new_seq2_index;

    /* the goal here is to choose the largest seq2 offset
       that leads to the current distance. There are three
       choices possible and each is checked in turn */

    if (diag >= diag_lower[(*d) - op_cost] &&
        diag <= diag_upper[(*d) - op_cost]) {

        new_seq2_index = last_seq2_off[(*d) - op_cost][diag].match_off;
        if (new_seq2_index >= MAX(last_seq2_off[*d][diag].insert_off,
                                  last_seq2_off[*d][diag].delete_off)) {
            *d -= op_cost;
            *seq2_index = new_seq2_index;
            return eGapAlignSub;
        }
    }
    if (last_seq2_off[*d][diag].insert_off >
                        last_seq2_off[*d][diag].delete_off) {
        *seq2_index = last_seq2_off[*d][diag].insert_off;
        return eGapAlignIns;
    }
    else {
        *seq2_index = last_seq2_off[*d][diag].delete_off;
        return eGapAlignDel;
    }
}

/** During the traceback for a greedy alignment with affine
    gap penalties, determine the next state of the traceback after
    moving upwards in the traceback array from an insertion or deletion

    @param last_seq2_off Array of offsets into the second sequence;
                        last_seq2_off[d][k] gives the largest offset into
                        the second sequence that lies on diagonal k and
                        has distance d [in]
    @param diag_lower Array of lower bounds on diagonal index [in]
    @param diag_upper Array of upper bounds on diagonal index [in]
    @param d Starting distance [in][out]
    @param diag The diagonal of the current traceback position [in]
    @param gap_open open a gap [in]
    @param gap_extend (Modified) cost to extend a gap [in]
    @param IorD The state of the traceback at present [in]
    @return The state for the next traceback operation
*/
static inline EGapAlignOpType
s_GetNextAffineTbackFromIndel(SGreedyOffset** last_seq2_off, Int4* diag_lower,
                    Int4* diag_upper, Int4* d, Int4 diag, Int4 gap_open,
                    Int4 gap_extend, EGapAlignOpType IorD)
{
    Int4 new_diag;
    Int4 new_seq2_index;
    Int4 gap_open_extend = gap_open + gap_extend;
    Int4 last_d;

    /* as with the previous routine, the traceback operation
       that leads to the current one must be determined. Either
       a gap is opened from a previous run of matches, or a
       gap is extended. If several traceback choices are
       available, choose the one that generates the largest
       seq2 offset */

    /* if the previous traceback operation is a gap, then it
       starts one diagonal away from the current diagonal... */

    if (IorD == eGapAlignIns)
        new_diag = diag - 1;
    else
        new_diag = diag + 1;

    /* ...and its seq2 offset is fixed */

    last_d = (*d) - gap_extend;
    if (new_diag >= diag_lower[last_d] &&
        new_diag <= diag_upper[last_d]) {

        if (IorD == eGapAlignIns)
            new_seq2_index =
                    last_seq2_off[last_d][new_diag].insert_off;
        else
            new_seq2_index =
                    last_seq2_off[last_d][new_diag].delete_off;
    }
    else {
        /* signal that no gap can be extended
           to achieve distance d */
        new_seq2_index = kInvalidOffset;
    }

    /* make the next traceback operation a match if it is
       valid to do so and the resulting seq2 offset exceeds the
       one derived from extending a gap */

    last_d = (*d) - gap_open_extend;
    if (new_diag >= diag_lower[last_d] &&
        new_diag <= diag_upper[last_d] &&
        new_seq2_index < last_seq2_off[last_d][new_diag].match_off) {

        *d -= gap_open_extend;
        return eGapAlignSub;
    }

    ASSERT(new_seq2_index != kInvalidOffset);
    *d -= gap_extend;
    return IorD;
}

//==============================================================================
/** Perform the greedy extension algorithm with non-affine gap penalties.
 * @param seq1 First sequence (always uncompressed) [in]
 * @param len1 Maximal extension length in first sequence [in]
 * @param seq2 Second sequence (may be compressed) [in]
 * @param len2 Maximal extension length in second sequence [in]
 * @param reverse Is extension performed in backwards direction? [in]
 * @param xdrop_threshold X-dropoff value to use in extension [in]
 * @param match_cost Match score to use in extension [in]
 * @param mismatch_cost Mismatch score to use in extension [in]
 * @param seq1_align_len Length of extension on sequence 1 [out]
 * @param seq2_align_len Length of extension on sequence 2 [out]
 * @param aux_data Structure containing all preallocated memory [in]
 * @param edit_block Edit script structure for saving traceback.
 *          Traceback is not saved if NULL is passed. [in] [out]
 * @param rem Offset within a byte of the compressed second sequence.
 *          Set to 4 if sequence is uncompressed. [in]
 * @param seed Structure to remember longest run of exact matches [out]
 * @return The minimum distance between the two sequences, i.e.
 *          the number of mismatches plus gaps in the resulting alignment
 */
Int4
BLAST_GreedyAlign (const Uint1* seq1, Int4 len1,
                   const Uint1* seq2, Int4 len2,
                   Boolean reverse, Int4 xdrop_threshold,
                   Int4 match_cost, Int4 mismatch_cost,
                   Int4* seq1_align_len, Int4* seq2_align_len,
                   SGreedyAlignMem* aux_data,
                   GapPrelimEditBlock *edit_block, Uint1 rem,
                   SGreedySeed *seed);

//==============================================================================
/** Perform the greedy extension algorithm with affine gap penalties.
 * @param seq1 First sequence (always uncompressed) [in]
 * @param len1 Maximal extension length in first sequence [in]
 * @param seq2 Second sequence (may be compressed) [in]
 * @param len2 Maximal extension length in second sequence [in]
 * @param reverse Is extension performed in backwards direction? [in]
 * @param xdrop_threshold X-dropoff value to use in extension [in]
 * @param match_cost Match score to use in extension [in]
 * @param mismatch_cost Mismatch score to use in extension [in]
 * @param in_gap_open Gap opening penalty [in]
 * @param in_gap_extend Gap extension penalty [in]
 * @param seq1_align_len Length of extension on sequence 1 [out]
 * @param seq2_align_len Length of extension on sequence 2 [out]
 * @param aux_data Structure containing all preallocated memory [in]
 * @param edit_block Edit script structure for saving traceback.
 *          Traceback is not saved if NULL is passed. [in] [out]
 * @param rem Offset within a byte of the compressed second sequence.
 *          Set to 4 if sequence is uncompressed. [in]
 * @param seed Structure to remember longest run of exact matches [out]
 * @return The score of the alignment
 */
Int4
BLAST_AffineGreedyAlign (const Uint1* seq1, Int4 len1,
                         const Uint1* seq2, Int4 len2,
                         Boolean reverse, Int4 xdrop_threshold,
                         Int4 match_cost, Int4 mismatch_cost,
                         Int4 in_gap_open, Int4 in_gap_extend,
                         Int4* seq1_align_len, Int4* seq2_align_len,
                         SGreedyAlignMem* aux_data,
                         GapPrelimEditBlock *edit_block, Uint1 rem,
                         SGreedySeed *seed);

//==============================================================================
Int4 BLAST_GreedyAlign(const Uint1* seq1, Int4 len1,
                       const Uint1* seq2, Int4 len2,
                       Boolean reverse, Int4 xdrop_threshold,
                       Int4 match_cost, Int4 mismatch_cost,
                       Int4* seq1_align_len, Int4* seq2_align_len,
                       SGreedyAlignMem* aux_data,
                       GapPrelimEditBlock *edit_block, Uint1 rem,
                       SGreedySeed *seed)
{
    Int4 seq1_index;
    Int4 seq2_index;
    Int4 index;
    Int4 d;
    Int4 k;
    Int4 diag_lower, diag_upper;
    Int4 max_dist;
    Int4 diag_origin;
    Int4 best_dist;
    Int4 best_diag;
    Int4** last_seq2_off;
    Int4* max_score;
    Int4 xdrop_offset;
    Int4 longest_match_run;
    Boolean end1_reached, end2_reached;
    //SMBSpace* mem_pool;
    MBSpaceMgr* mem_mgr;
    Boolean converged = FALSE;

    /* ordinary dynamic programming alignment, for each offset
       in seq1, walks through offsets in seq2 until an X-dropoff
       test fails, saving the best score encountered along
       the way. Instead of score, this code tracks the 'distance'
       (number of mismatches plus number of gaps) between seq1
       and seq2. Instead of walking through sequence offsets, it
       walks through diagonals that can achieve a given distance.

       Note that in what follows, the numbering of diagonals implies
       a dot matrix where increasing seq1 offsets go to the right on
       the x axis, and increasing seq2 offsets go up the y axis.
       The gapped alignment thus proceeds up and to the right in
       the graph, and diagonals are numbered increasing to the right */

    best_dist = 0;
    best_diag = 0;

    /* set the number of distinct distances the algorithm will
       examine in the search for an optimal alignment. The
       heuristic worst-case running time of the algorithm is
       O(max_dist**2 + (len1+len2)); for sequences which are
       very similar, the average running time will be sig-
       nificantly better than this */

    max_dist = aux_data->max_dist;

    /* the main loop assumes that the index of all diagonals is
       biased to lie in the middle of allocated bookkeeping
       structures */

    diag_origin = max_dist + 2;

    /* last_seq2_off[d][k] is the largest offset into seq2 that
       lies on diagonal k and has distance d */

    last_seq2_off = aux_data->last_seq2_off;

    /* Instead of tracking the best alignment score and using
       xdrop_theshold directly, track the best score for each
       unique distance and use the best score for some previously
       computed distance to implement the X-dropoff test.

       xdrop_offset gives the distance backwards in the score
       array to look */

    xdrop_offset = (xdrop_threshold + match_cost / 2) /
                           (match_cost + mismatch_cost) + 1;

    /* find the offset of the first mismatch between seq1 and seq2 */

    index = s_FindFirstMismatch(seq1, seq2, len1, len2, 0, 0,
                                reverse, rem);

    /* update the extents of the alignment, and bail out
       early if no further work is needed */

    *seq1_align_len = index;
    *seq2_align_len = index;
    seq1_index = index;

    seed->start_q = 0;
    seed->start_s = 0;
    seed->match_length = longest_match_run = index;

    if (index == len1 || index == len2) {
        /* Return the number of differences, which is zero here */

        if (edit_block != NULL)
            //GapPrelimEditBlockAdd(edit_block, eGapAlignSub, index);
            edit_block->Add(eGapAlignSub, index);
        return 0;
    }

    /* set up the memory pool */

//    mem_pool = aux_data->space;
//    if (edit_block == NULL) {
//       mem_pool = NULL;
//    }
//    else if (mem_pool == NULL) {
//       aux_data->space = mem_pool = MBSpaceNew(0);
//    }
//    else {
//        s_RefreshMBSpace(mem_pool);
//    }

    mem_mgr = aux_data->spacemgr;
    if (edit_block == NULL) {
        mem_mgr = NULL;
    } else if (mem_mgr == NULL) {
        aux_data->spacemgr = mem_mgr = new MBSpaceMgr();
    } else {
        mem_mgr->RefreshMBSpace();
    }


    /* set up the array of per-distance maximum scores. There
       are max_diags + xdrop_offset distances to track, the first
       xdrop_offset of which are 0 */

    max_score = aux_data->max_score + xdrop_offset;
    //for (index = 0; index < xdrop_offset; index++)
    //    aux_data->max_score[index] = 0;
    memset(aux_data->max_score, 0, sizeof(Int4) * xdrop_offset);

    /* fill in the initial offsets of the distance matrix */

    last_seq2_off[0][diag_origin] = seq1_index;
    max_score[0] = seq1_index * match_cost;
    diag_lower = diag_origin - 1;
    diag_upper = diag_origin + 1;
    end1_reached = end2_reached = FALSE;

    /* for each distance */
    for (d = 1; d <= max_dist; d++) {
        Int4 xdrop_score;
        Int4 curr_score;
        Int4 curr_extent = 0;
        Int4 curr_seq2_index = 0;
        Int4 curr_diag = 0;
        Int4 tmp_diag_lower = diag_lower;
        Int4 tmp_diag_upper = diag_upper;

        /* assign impossible seq2 offsets to any diagonals that
           are not in the range (diag_lower,diag_upper).
           These will serve as sentinel values for the
           inner loop */

        last_seq2_off[d - 1][diag_lower-1] = kInvalidOffset;
        last_seq2_off[d - 1][diag_lower] = kInvalidOffset;
        last_seq2_off[d - 1][diag_upper] = kInvalidOffset;
        last_seq2_off[d - 1][diag_upper+1] = kInvalidOffset;

        /* compute the score for distance d that corresponds to
           the X-dropoff criterion */

        xdrop_score = max_score[d - xdrop_offset] +
                      (match_cost + mismatch_cost) * d - xdrop_threshold;
        xdrop_score = (Int4)ceil((double)xdrop_score / (match_cost / 2));

        /* for each diagonal of interest */

        for (k = tmp_diag_lower; k <= tmp_diag_upper; k++) {

            /* find the largest offset into seq2 that increases
               the distance from d-1 to d (i.e. keeps the alignment
               from getting worse for as long as possible), then
               choose the offset into seq1 that will keep the
               resulting diagonal fixed at k

               Note that this requires kInvalidOffset+1 to be smaller
               than any valid offset into seq2, i.e. to be negative */

            seq2_index = MAX(last_seq2_off[d - 1][k + 1],
                             last_seq2_off[d - 1][k    ]) + 1;
            seq2_index = MAX(seq2_index, last_seq2_off[d - 1][k - 1]);
            seq1_index = seq2_index + k - diag_origin;

            if (seq2_index < 0 || seq1_index + seq2_index < xdrop_score) {

                /* if no valid diagonal can reach distance d, or the
                   X-dropoff test fails, narrow the range of diagonals
                   to test and skip to the next diagonal */

                if (k == diag_lower)
                    diag_lower++;
                else
                    last_seq2_off[d][k] = kInvalidOffset;
                continue;
            }
            diag_upper = k;

            /* slide down diagonal k until a mismatch
               occurs. As long as only matches are encountered,
               the current distance d will not change */

            index = s_FindFirstMismatch(seq1, seq2, len1, len2,
                                        seq1_index, seq2_index,
                                        reverse, rem);

            if (index > longest_match_run) {
                seed->start_q = seq1_index;
                seed->start_s = seq2_index;
                seed->match_length = longest_match_run = index;
            }
            seq1_index += index;
            seq2_index += index;

            /* set the new largest seq2 offset that achieves
               distance d on diagonal k */

            last_seq2_off[d][k] = seq2_index;

            /* since all values of k are constrained to have the
               same distance d, the value of k which maximizes the
               alignment score is the one that covers the most
               of seq1 and seq2 */

            if (seq1_index + seq2_index > curr_extent) {
                curr_extent = seq1_index + seq2_index;
                curr_seq2_index = seq2_index;
                curr_diag = k;
            }

            /* clamp the bounds on diagonals to avoid walking off
               either sequence. Because the bounds increase by at
               most one for each distance, diag_lower and diag_upper
               can each be of size at most max_diags+2 */

            if (seq2_index == len2) {
                diag_lower = k + 1;
                end2_reached = TRUE;
            }
            if (seq1_index == len1) {
                diag_upper = k - 1;
                end1_reached = TRUE;
            }
        }   /* end loop over diagonals */

        /* compute the maximum score possible for distance d */

        curr_score = curr_extent * (match_cost / 2) -
                        d * (match_cost + mismatch_cost);

        /* if this is the best score seen so far, update the
           statistics of the best alignment */

        if (curr_score > max_score[d - 1]) {
            max_score[d] = curr_score;
            best_dist = d;
            best_diag = curr_diag;
            *seq2_align_len = curr_seq2_index;
            *seq1_align_len = curr_seq2_index + best_diag - diag_origin;
        }
        else {
            max_score[d] = max_score[d - 1];
        }

        /* alignment has finished if the lower and upper bounds
           on diagonals to check have converged to each other */

        if (diag_lower > diag_upper) {
            converged = TRUE;
            break;
        }

        /* set up for the next distance to examine. Because the
           bounds increase by at most one for each distance,
           diag_lower and diag_upper can each be of size at
           most max_diags+2 */

        if (!end2_reached)
            diag_lower--;
        if (!end1_reached)
            diag_upper++;

        /* if no traceback is specified, the next row of
           last_seq2_off can reuse previously allocated memory */

        if (edit_block == NULL) {

            /** @todo FIXME The following assumes two arrays of
                at least max_dist+4 Int4's have already been allocated */

            last_seq2_off[d + 1] = last_seq2_off[d - 1];
        }
        else {

            /* traceback requires all rows of last_seq2_off to be saved,
               so a new row must be allocated. The allocator provides
               SThreeVal structures which are 3 times larger than Int4,
               so divide requested amount by 3 */

            /** @todo FIXME Should make allocator more general */

            //last_seq2_off[d + 1] = (Int4*) s_GetMBSpace(mem_pool,
            //                         (diag_upper - diag_lower + 7) / 3);
            last_seq2_off[d + 1] = (Int4*)mem_mgr->GetMBSpace((diag_upper - diag_lower + 7) / 3);

            /* move the origin for this row backwards */

            last_seq2_off[d + 1] = last_seq2_off[d + 1] - diag_lower + 2;
        }
    }   /* end loop over distinct distances */

    if (!converged)
        return -1;

    if (edit_block == NULL)
        return best_dist;

    /* perform traceback */

    d = best_dist;
    seq1_index = *seq1_align_len;
    seq2_index = *seq2_align_len;

    /* for all positive distances */

    while (d > 0) {
        Int4 new_diag;
        Int4 new_seq2_index;

        /* retrieve the value of the diagonal after the next
           traceback operation. best_diag starts off with the
           value computed during the alignment process */

        new_diag = s_GetNextNonAffineTback(last_seq2_off, d,
                                           best_diag, &new_seq2_index);

        if (new_diag == best_diag) {

            /* same diagonal: issue a group of substitutions */

            if (seq2_index - new_seq2_index > 0) {
                //GapPrelimEditBlockAdd(edit_block, eGapAlignSub,
                //                seq2_index - new_seq2_index);
                edit_block->Add(eGapAlignSub, seq2_index - new_seq2_index);
            }
        }
        else if (new_diag < best_diag) {

            /* smaller diagonal: issue a group of substitutions
               and then a gap in seq2 */

            if (seq2_index - new_seq2_index > 0) {
                //GapPrelimEditBlockAdd(edit_block, eGapAlignSub,
                //                seq2_index - new_seq2_index);
                edit_block->Add(eGapAlignSub, seq2_index - new_seq2_index);
            }
            //GapPrelimEditBlockAdd(edit_block, eGapAlignIns, 1);
            edit_block->Add(eGapAlignIns, 1);
        }
        else {
            /* larger diagonal: issue a group of substitutions
               and then a gap in seq1 */

            if (seq2_index - new_seq2_index - 1 > 0) {
                //GapPrelimEditBlockAdd(edit_block, eGapAlignSub,
                //                seq2_index - new_seq2_index - 1);
                edit_block->Add(eGapAlignSub, seq2_index - new_seq2_index - 1);
            }
            //GapPrelimEditBlockAdd(edit_block, eGapAlignDel, 1);
            edit_block->Add(eGapAlignDel, 1);
        }
        d--;
        best_diag = new_diag;
        seq2_index = new_seq2_index;
    }

done:
    /* handle the final group of substitutions back to distance zero,
       i.e. back to offset zero of seq1 and seq2 */

    //GapPrelimEditBlockAdd(edit_block, eGapAlignSub,
    //                      last_seq2_off[0][diag_origin]);
    edit_block->Add(eGapAlignSub, last_seq2_off[0][diag_origin]);

    return best_dist;
}

//==============================================================================
Int4 BLAST_AffineGreedyAlign (const Uint1* seq1, Int4 len1,
                              const Uint1* seq2, Int4 len2,
                              Boolean reverse, Int4 xdrop_threshold,
                              Int4 match_score, Int4 mismatch_score,
                              Int4 in_gap_open, Int4 in_gap_extend,
                              Int4* seq1_align_len, Int4* seq2_align_len,
                              SGreedyAlignMem* aux_data,
                              GapPrelimEditBlock *edit_block, Uint1 rem,
                              SGreedySeed *seed)
{
    Int4 seq1_index;
    Int4 seq2_index;
    Int4 index;
    Int4 d;
    Int4 k;
    Int4 max_dist;
    Int4 scaled_max_dist;
    Int4 diag_origin;
    Int4 best_dist;
    Int4 best_diag;
    Int4 longest_match_run;
    SGreedyOffset** last_seq2_off;
    Int4* max_score;
    Int4 xdrop_offset;
    Int4 end1_diag, end2_diag;
    //SMBSpace* mem_pool;
    MBSpaceMgr* mem_mgr;

    Int4 op_cost;
    Int4 gap_open;
    Int4 gap_extend;
    Int4 gap_open_extend;
    Int4 max_penalty;
    Int4 score_common_factor;
    Int4 match_score_half;

    Int4 *diag_lower;
    Int4 *diag_upper;
    Int4 curr_diag_lower;
    Int4 curr_diag_upper;

    Int4 num_nonempty_dist;
    const Int4 kInvalidDiag = 100000000; /* larger than any valid diag. index */
    Boolean converged = FALSE;

    /* make sure bits of match_score don't disappear if it
       is divided by 2 */

    if (match_score % 2 == 1) {
        match_score *= 2;
        mismatch_score *= 2;
        xdrop_threshold *= 2;
        in_gap_open *= 2;
        in_gap_extend *= 2;
    }

    if (in_gap_open == 0 && in_gap_extend == 0) {
       return BLAST_GreedyAlign(seq1, len1, seq2, len2, reverse,
                                xdrop_threshold, match_score,
                                mismatch_score, seq1_align_len,
                                seq2_align_len, aux_data, edit_block,
                                rem, seed);
    }

    /* ordinary dynamic programming alignment, for each offset
       in seq1, walks through offsets in seq2 until an X-dropoff
       test fails, saving the best score encountered along
       the way. Instead of score, this code tracks the 'distance'
       (number of mismatches plus number of gaps) between seq1
       and seq2. Instead of walking through sequence offsets, it
       walks through diagonals that can achieve a given distance

       Note that in what follows, the numbering of diagonals implies
       a dot matrix where increasing seq1 offsets go to the right on
       the x axis, and increasing seq2 offsets go up the y axis.
       The gapped alignment thus proceeds up and to the right in
       the graph, and diagonals are numbered increasing to the right */

    best_dist = 0;
    best_diag = 0;

    /* fill in derived scores and penalties */

    match_score_half = match_score / 2;
    op_cost = match_score + mismatch_score;
    gap_open = in_gap_open;
    gap_extend = in_gap_extend + match_score_half;
    score_common_factor = Gdb3(&op_cost, &gap_open, &gap_extend);
    gap_open_extend = gap_open + gap_extend;
    max_penalty = MAX(op_cost, gap_open_extend);

    /* set the number of distinct distances the algorithm will
       examine in the search for an optimal alignment */

    max_dist = aux_data->max_dist;
    scaled_max_dist = max_dist * gap_extend;

    /* the main loop assumes that the index of all diagonals is
       biased to lie in the middle of allocated bookkeeping structures */

    diag_origin = max_dist + 2;

    /* last_seq2_off[d][k] is the largest offset into seq2 that
       lies on diagonal k and has distance d. Unlike the non-affine
       case, the largest offset for paths ending in an insertion,
       deletion, and match must all be separately saved for
       each d and k */

    last_seq2_off = aux_data->last_seq2_off_affine;

    /* Instead of tracking the best alignment score and using
       xdrop_theshold directly, track the best score for each
       unique distance and use the best score for some previously
       computed distance to implement the X-dropoff test.

       xdrop_offset gives the distance backwards in the score
       array to look */

    xdrop_offset = (xdrop_threshold + match_score_half) /
                                      score_common_factor + 1;

    /* find the offset of the first mismatch between seq1 and seq2 */

    index = s_FindFirstMismatch(seq1, seq2, len1, len2, 0, 0,
                                reverse, rem);

    /* update the extents of the alignment, and bail out
       early if no further work is needed */

    *seq1_align_len = index;
    *seq2_align_len = index;
    seq1_index = index;

    seed->start_q = 0;
    seed->start_s = 0;
    seed->match_length = longest_match_run = index;

    if (index == len1 || index == len2) {
        /* return the score of the run of matches */
        if (edit_block != NULL)
            //GapPrelimEditBlockAdd(edit_block, eGapAlignSub, index);
            edit_block->Add(eGapAlignSub, index);
       return index * match_score;
    }

    /* set up the memory pool */

//    mem_pool = aux_data->space;
//    if (edit_block == NULL) {
//        mem_pool = NULL;
//    }
//    else if (!mem_pool) {
//       aux_data->space = mem_pool = MBSpaceNew(0);
//    }
//    else {
//        s_RefreshMBSpace(mem_pool);
//    }

    mem_mgr = aux_data->spacemgr;
    if (edit_block == NULL) {
        mem_mgr = NULL;
    } else if (!mem_mgr) {
        aux_data->spacemgr = mem_mgr = new MBSpaceMgr();
    } else {
        mem_mgr->RefreshMBSpace();
    }

    /* set up the array of per-distance maximum scores. There
       are scaled_max_dist + xdrop_offset distances to track,
       the first xdrop_offset of which are 0 */

    max_score = aux_data->max_score + xdrop_offset;
    //for (index = 0; index < xdrop_offset; index++)
    //    aux_data->max_score[index] = 0;
    memset(aux_data->max_score, 0, sizeof(Int4) * xdrop_offset);

    /* For affine greedy alignment, contributions to distance d
       can come from distances further back than d-1 (which is
       sufficient for non-affine alignment). Where non-affine
       alignment only needs to track the current bounds on diagonals
       to test, the present code must also track the upper and
       lower bounds on diagonals for max_penalty previous distances.
       These share the same preallocated array */

    diag_lower = aux_data->diag_bounds;
    diag_upper = aux_data->diag_bounds +
                 scaled_max_dist + 1 + max_penalty;

    /* the first max_penalty elements correspond to negative
       distances; initialize with an empty range of diagonals */

    for (index = 0; index < max_penalty; index++) {
        diag_lower[index] = kInvalidDiag;
        diag_upper[index] = -kInvalidDiag;
    }
    diag_lower += max_penalty;
    diag_upper += max_penalty;

    /* fill in the statistics for distance zero, i.e. the initial
       run of exact matches */

    last_seq2_off[0][diag_origin].match_off = seq1_index;
    last_seq2_off[0][diag_origin].insert_off = kInvalidOffset;
    last_seq2_off[0][diag_origin].delete_off = kInvalidOffset;
    max_score[0] = seq1_index * match_score;
    diag_lower[0] = diag_origin;
    diag_upper[0] = diag_origin;

    /* set up for distance 1 */

    curr_diag_lower = diag_origin - 1;
    curr_diag_upper = diag_origin + 1;
    end1_diag = 0;
    end2_diag = 0;
    num_nonempty_dist = 1;
    d = 1;

    /* for each distance */

    while (d <= scaled_max_dist) {
        Int4 xdrop_score;
        Int4 curr_score;
        Int4 curr_extent = 0;
        Int4 curr_seq2_index = 0;
        Int4 curr_diag = 0;
        Int4 tmp_diag_lower = curr_diag_lower;
        Int4 tmp_diag_upper = curr_diag_upper;

        /* compute the score for distance d that corresponds to
           the X-dropoff criterion */

        xdrop_score = max_score[d - xdrop_offset] +
                      score_common_factor * d - xdrop_threshold;
        xdrop_score = (Int4)ceil(1.0 * xdrop_score / match_score_half);
        if (xdrop_score < 0)
            xdrop_score = 0;

        /* for each valid diagonal */

        for (k = tmp_diag_lower; k <= tmp_diag_upper; k++) {

            /* As with the non-affine algorithm, the object is
               to find the largest offset into seq2 that can
               achieve distance d from diagonal k. Here, however,
               contributions are possible from distances < d-1 */

            /* begin by assuming the best offset comes from opening
               a gap in seq1. Since opening a gap costs gap_open_extend,
               use the offset associated with a match from that
               far back in the table. Do not use diagonal k+1 if
               it was not valid back then */

            seq2_index = kInvalidOffset;
            if (k + 1 <= diag_upper[d - gap_open_extend] &&
                k + 1 >= diag_lower[d - gap_open_extend]) {
                seq2_index = last_seq2_off[d - gap_open_extend][k+1].match_off;
            }

            /* Replace with the offset derived from extending a gap
               in seq1, if that is larger */

            if (k + 1 <= diag_upper[d - gap_extend] &&
                k + 1 >= diag_lower[d - gap_extend] &&
                seq2_index < last_seq2_off[d - gap_extend][k+1].delete_off) {
                seq2_index = last_seq2_off[d - gap_extend][k+1].delete_off;
            }

            /* save the index; if it was valid, a deletion
               (= gap in seq1) means seq2 offset slips by one */

            if (seq2_index == kInvalidOffset)
                last_seq2_off[d][k].delete_off = kInvalidOffset;
            else
                last_seq2_off[d][k].delete_off = seq2_index + 1;

            /* repeat the process assuming a gap is opened or
               extended in seq2. Gaps in seq2 do not change seq2_index */

            seq2_index = kInvalidOffset;
            if (k - 1 <= diag_upper[d - gap_open_extend] &&
                k - 1 >= diag_lower[d - gap_open_extend]) {
                seq2_index = last_seq2_off[d - gap_open_extend][k-1].match_off;
            }
            if (k - 1 <= diag_upper[d - gap_extend] &&
                k - 1 >= diag_lower[d - gap_extend] &&
                seq2_index < last_seq2_off[d - gap_extend][k-1].insert_off) {
                seq2_index = last_seq2_off[d - gap_extend][k-1].insert_off;
            }
            last_seq2_off[d][k].insert_off = seq2_index;

            /* Compare the greater of the two previous answers with
               the offset associated with a match on diagonal k. */

            seq2_index = MAX(last_seq2_off[d][k].insert_off,
                             last_seq2_off[d][k].delete_off);
            if (k <= diag_upper[d - op_cost] &&
                k >= diag_lower[d - op_cost]) {
                seq2_index = MAX(seq2_index,
                                 last_seq2_off[d - op_cost][k].match_off + 1);
            }

            /* choose the offset into seq1 so as to remain on diagonal k */

            seq1_index = seq2_index + k - diag_origin;

            /* perform the X-dropoff test; if it fails, or no
               previous cell can contribute to the current one,
               give up and try the next diagonal, adjusting the
               bounds on diagonals for distance d */

            if (seq2_index < 0 || seq1_index + seq2_index < xdrop_score) {
                if (k == curr_diag_lower)
                    curr_diag_lower++;
                else
                    last_seq2_off[d][k].match_off = kInvalidOffset;
                continue;
            }
            curr_diag_upper = k;

            /* slide down diagonal k until a mismatch
               occurs. As long as only matches are encountered,
               the current distance d will not change */

            index = s_FindFirstMismatch(seq1, seq2, len1, len2,
                                        seq1_index, seq2_index,
                                        reverse, rem);

            if (index > longest_match_run) {
                seed->start_q = seq1_index;
                seed->start_s = seq2_index;
                seed->match_length = longest_match_run = index;
            }
            seq1_index += index;
            seq2_index += index;

            /* since all values of k are constrained to have the
               same distance d, the value of k which maximizes the
               alignment score is the one that covers the most
               of seq1 and seq2 */

            last_seq2_off[d][k].match_off = seq2_index;
            if (seq1_index + seq2_index > curr_extent) {
                curr_extent = seq1_index + seq2_index;
                curr_seq2_index = seq2_index;
                curr_diag = k;
            }

            /* clamp the bounds on diagonals to avoid walking off
               either sequence */

            if (seq1_index == len1) {
                curr_diag_upper = k;
                end1_diag = k - 1;
            }
            if (seq2_index == len2) {
                curr_diag_lower = k;
                end2_diag = k + 1;
            }
        }  /* end loop over diagonals */

        /* compute the maximum score possible for distance d */

        curr_score = curr_extent * match_score_half - d * score_common_factor;

        /* if this is the best score seen so far, update the
           statistics of the best alignment */

        if (curr_score > max_score[d - 1]) {
            max_score[d] = curr_score;
            best_dist = d;
            best_diag = curr_diag;
            *seq2_align_len = curr_seq2_index;
            *seq1_align_len = curr_seq2_index + best_diag - diag_origin;
        }
        else {
            max_score[d] = max_score[d - 1];
        }

        /* save the bounds on diagonals to examine for distance d.
           Note that in the non-affine case the alignment could stop
           if these bounds converged to each other. Here, however,
           it's possible for distances less than d to continue the
           alignment even if no diagonals are available at distance d.
           Hence we can only stop if max_penalty consecutive ranges
           of diagonals are empty */

        if (curr_diag_lower <= curr_diag_upper) {
            num_nonempty_dist++;
            diag_lower[d] = curr_diag_lower;
            diag_upper[d] = curr_diag_upper;
        }
        else {
            diag_lower[d] = kInvalidDiag;
            diag_upper[d] = -kInvalidDiag;
        }

        if (diag_lower[d - max_penalty] <= diag_upper[d - max_penalty])
            num_nonempty_dist--;

        if (num_nonempty_dist == 0) {
            converged = TRUE;
            break;
        }

        /* compute the range of diagonals to test for the next
           value of d. These must be conservative, in that any
           diagonal that could possibly contribute must be allowed.
           curr_diag_lower and curr_diag_upper can each be of size at
           most scaled_max_diags+2; they can also represent an
           empty range, in which case the next value of d will never
           improve the best score */

        d++;
        curr_diag_lower = MIN(diag_lower[d - gap_open_extend],
                              diag_lower[d - gap_extend]) - 1;
        curr_diag_lower = MIN(curr_diag_lower, diag_lower[d - op_cost]);

        if (end2_diag > 0)
            curr_diag_lower = MAX(curr_diag_lower, end2_diag);

        curr_diag_upper = MAX(diag_upper[d - gap_open_extend],
                              diag_upper[d - gap_extend]) + 1;
        curr_diag_upper = MAX(curr_diag_upper,
                              diag_upper[d - op_cost]);

        if (end1_diag > 0)
            curr_diag_upper = MIN(curr_diag_upper, end1_diag);

        if (d > max_penalty) {
            if (edit_block == NULL) {

                /* if no traceback is required, the next row of
                   last_seq2_off can reuse previously allocated memory */

                last_seq2_off[d] = last_seq2_off[d - max_penalty - 1];
            }
            else {

                /* traceback requires all rows of last_seq2_off to be saved,
                   so a new row must be allocated */

//                last_seq2_off[d] = s_GetMBSpace(mem_pool,
//                                   curr_diag_upper - curr_diag_lower + 1) -
//                                   curr_diag_lower;

                  last_seq2_off[d] = mem_mgr->GetMBSpace(curr_diag_upper - curr_diag_lower + 1) - curr_diag_lower;
            }
        }
    }  /* end loop over distances */

    if (! converged) return -1;

    /* compute the traceback if necessary */

    if (edit_block != NULL) {
        EGapAlignOpType state;

        d = best_dist;
        seq2_index = *seq2_align_len;
        state = eGapAlignSub;

        while (d > 0) {
            if (state == eGapAlignSub) {
                /* substitution */
                Int4 new_seq2_index;
                state = s_GetNextAffineTbackFromMatch(last_seq2_off,
                                       diag_lower, diag_upper, &d, best_diag,
                                       op_cost, &new_seq2_index);

                ASSERT(seq2_index > new_seq2_index);
                //GapPrelimEditBlockAdd(edit_block, eGapAlignSub,
                //                    seq2_index - new_seq2_index);
                edit_block->Add(eGapAlignSub, seq2_index - new_seq2_index);
                seq2_index = new_seq2_index;
            }
            else if (state == eGapAlignIns) {
                /* gap in seq2 */
                //GapPrelimEditBlockAdd(edit_block, eGapAlignIns, 1);
                edit_block->Add(eGapAlignIns, 1);
                state = s_GetNextAffineTbackFromIndel(last_seq2_off,
                                     diag_lower, diag_upper, &d, best_diag,
                                     gap_open, gap_extend, eGapAlignIns);
                best_diag--;
            }
            else {
                /* gap in seq1 */
                //GapPrelimEditBlockAdd(edit_block, eGapAlignDel, 1);
                edit_block->Add(eGapAlignDel, 1);
                state = s_GetNextAffineTbackFromIndel(last_seq2_off,
                                     diag_lower, diag_upper, &d, best_diag,
                                     gap_open, gap_extend, eGapAlignDel);
                best_diag++;
                seq2_index--;
            }
        }

        /* write the last group of matches */
        //GapPrelimEditBlockAdd(edit_block, eGapAlignSub,
        //                      last_seq2_off[0][diag_origin].match_off);
        edit_block->Add(eGapAlignSub, last_seq2_off[0][diag_origin].match_off);
    }

    return max_score[best_dist];
}


//==============================================================================
static void s_UpdateEditScript(GapEditScript* esp, int pos, int bf, int af) {
   int op, qd, sd;

   if (bf > 0) {
      op = pos;
      qd = sd = bf;
      do {
          if (--op < 0) return;
          switch(esp->op_type[op]) {
          case eGapAlignSub:
              qd -= esp->num[op];
              sd -= esp->num[op];
              break;
          case eGapAlignIns:
              qd -= esp->num[op];
              break;
          case eGapAlignDel:
              sd -= esp->num[op];
          default:
              break;
          }
      } while (qd > 0 || sd > 0);

      esp->num[op] = -MAX(qd, sd);
      esp->op_type[op++] = eGapAlignSub;
      for (; op < pos-1; op++) esp->num[op] = 0;
      esp->num[pos] += bf;
      qd -= sd;
      esp->op_type[pos-1] = (qd>0) ? eGapAlignDel: eGapAlignIns;
      esp->num[pos-1] = (qd>0) ? qd : -qd;
   }

   if (af > 0) {
      op = pos;
      qd = sd = af;
      do {
          if (++op >= esp->size) return;
          switch(esp->op_type[op]) {
          case eGapAlignSub:
              qd -= esp->num[op];
              sd -= esp->num[op];
              break;
          case eGapAlignIns:
              qd -= esp->num[op];
              break;
          case eGapAlignDel:
              sd -= esp->num[op];
          default:
              break;
          }
      } while (qd > 0 || sd > 0);

      esp->num[op] = -MAX(qd, sd);
      esp->op_type[op--] = eGapAlignSub;
      for (; op > pos+1; op--) esp->num[op] = 0;
      esp->num[pos] += af;
      qd -= sd;
      esp->op_type[pos+1] = (qd>0) ? eGapAlignDel: eGapAlignIns;
      esp->num[pos+1] = (qd>0) ? qd : -qd;
   }
}

static void s_RebuildEditScript(GapEditScript* esp) {
   int i, j;
   for (i=0, j=-1; i<esp->size; i++) {
       if (esp->num[i] == 0) continue;
       if (j>=0 && esp->op_type[i] == esp->op_type[j]) {
           esp->num[j] += esp->num[i];
       } else if (j==-1 || esp->op_type[i] == eGapAlignSub
           || esp->op_type[j] == eGapAlignSub) {
           esp->op_type[++j] = esp->op_type[i];
           esp->num[j] = esp->num[i];
       } else {
           int d = esp->num[j] - esp->num[i];
           if (d > 0) {
              esp->num[j-1] += esp->num[i];
              esp->num[j] = d;
           } else if (d < 0) {
              esp->num[j-1] += esp->num[j];
              esp->num[j] = -d;
              esp->op_type[j] = esp->op_type[i];
           } else {
              esp->num[j-1] += esp->num[j];
              --j;
           }
       }
   }
   esp->size = ++j;
}

static void s_ReduceGaps(GapEditScript* esp, const Uint1 *q, const Uint1 *s,
                                             const Uint1 *qf,const Uint1 *sf){
   int i, j, nm1, nm2, d;
   const Uint1 *q1, *s1;

   for (q1=q, s1=s, i=0; i<esp->size; i++) {
       if (esp->num[i] == 0) continue;
       if (esp->op_type[i] == eGapAlignSub) {
           if(esp->num[i] >= 12) {
               nm1 = 1;
               if (i > 0) {
                   while (q1-nm1>=q && (*(q1-nm1) == *(s1-nm1))) ++nm1;
               }
               q1 += esp->num[i];
               s1 += esp->num[i];
               nm2 = 0;
               if (i < esp->size -1) {
                   //while (q1+1<qf && *(q1++) == *(s1++)) ++nm2;
				   while ((q1+1<qf) && (s1+1<sf) && (*(q1++) == *(s1++))) ++nm2;
               }
               if (nm1>1 || nm2>0) s_UpdateEditScript(esp, i, nm1-1, nm2);
               q1--; s1--;
           } else {
               q1 += esp->num[i];
               s1 += esp->num[i];
           }
       } else if (esp->op_type[i] == eGapAlignIns) {
           q1 += esp->num[i];
       } else {
           s1 += esp->num[i];
       }
   }
   s_RebuildEditScript(esp);

   for (i=0; i<esp->size; i++) {
       if (esp->op_type[i] == eGapAlignSub) {
           q += esp->num[i];
           s += esp->num[i];
           continue;
       }
       if (i>1 && esp->op_type[i] != esp->op_type[i-2]
               && esp->num[i-2] > 0) {
           d = esp->num[i] + esp->num[i-1] + esp->num[i-2];
           if (d == 3) {
               /* special case, no need to do further testing */
               (esp->num[i-2]) = 0;
               (esp->num[i-1]) = 2;
               (esp->num[i]) = 0;
               if (esp->op_type[i] == eGapAlignIns) {
                   ++q;
               } else {
                   ++s;
               }
           } else if (d < 12) {
               /* Try reducing this sub... */
               nm1 = 0;
               nm2 = 0;
               d = MIN(esp->num[i], esp->num[i-2]);
               q -= esp->num[i-1];
               s -= esp->num[i-1];
               q1 = q;
               s1 = s;
               if (esp->op_type[i] == eGapAlignIns) {
                   s -= d;
               } else {
                   q -= d;
               }
               for (j=0; j<esp->num[i-1]; ++j, ++q1, ++s1, ++q, ++s) {
                   if (*q1 == *s1) nm1++;
                   if (*q == *s) nm2++;
               }
               for (j=0; j<d; ++j, ++q, ++s) {
                   if (*q == *s) nm2++;
               }
               if (nm2 >= nm1 - d) {
                   (esp->num[i-2]) -= d;
                   (esp->num[i-1]) += d;
                   (esp->num[i]) -= d;
               } else {
                   q = q1;
                   s = s1;
               }
           }
       }
       if (esp->op_type[i] == eGapAlignIns) {
           q += esp->num[i];
       } else {
           s += esp->num[i];
       }
   }
   s_RebuildEditScript(esp);
}

#endif	/* GREEDY_ALIGN_IMPL_H */

