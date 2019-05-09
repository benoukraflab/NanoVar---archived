/*  $Id: win_mask_dup_table.hpp 183994 2010-02-23 20:20:11Z morgulis $
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
 *   Declaration of CheckDuplicates function.
 *
 */

#ifndef C_WIN_MASK_DUP_TABLE_HPP
#define C_WIN_MASK_DUP_TABLE_HPP

#include <string>
#include <vector>
#include <set>

//#include <objmgr/object_manager.hpp>
//#include <objmgr/scope.hpp>
//#include <objmgr/seq_entry_handle.hpp>
//#include <objmgr/bioseq_ci.hpp>
//#include <objmgr/seq_vector.hpp>
//#include <objmgr/util/sequence.hpp>
//#include <algo/winmask/win_mask_util.hpp>

// #include "win_mask_config.hpp"

BEGIN_NCBI_SCOPE

/**
 **\brief Check for possibly duplicate sequences in the input.
 **
 ** input contains the list of input file names. The files should be in
 ** the fasta format. The function checks the input sequences for
 ** duplication and reports possible duplicates to the standard error.
 ** 
 **\param input list of input file names
 **\param infmt input format
 **\param ids set of ids to check
 **\param exclude_ids set of ids to ignore
 **
 **/
//void CheckDuplicates( const vector< string > & input,
//                      const string & infmt,
//                      const CWinMaskUtil::CIdSet * ids,
//                      const CWinMaskUtil::CIdSet * exclude_ids );

END_NCBI_SCOPE

#endif
