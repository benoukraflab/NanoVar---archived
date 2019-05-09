/*  $Id: seq_masker_ostat_opt_ascii.hpp 462550 2015-03-19 14:07:19Z morgulis $
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
 *   Definition of CSeqMaskerOStatOptAscii class.
 *
 */

#ifndef C_SEQ_MASKER_OSTAT_OPT_ASCII_HPP
#define C_SEQ_MASKER_OSTAT_OPT_ASCII_HPP

#include "seq_masker_ostat_opt.hpp"

BEGIN_NCBI_SCOPE

/**
 **\brief This class is responsible for saving optimized unit counts
 **       in text format.
 **/
class NCBI_XALGOWINMASK_EXPORT
CSeqMaskerOstatOptAscii : public CSeqMaskerOstatOpt
{
    public:

        /// Format version.
        static CSeqMaskerVersion FormatVersion;

        /**
         **\brief Object constructor.
         **\param name output file name
         **\param sz requested upper limit on the size of the data structure
         **          (forwarded to CSeqMaskerOstatOpt)
         **/
        explicit CSeqMaskerOstatOptAscii( 
                const string & name, Uint2 sz, string const & metadata );

        /**
         **\brief Object constructor.
         **\param os the output stream 
         **\param sz requested upper limit on the size of the data structure
         **          (forwarded to CSeqMaskerOstatOpt)
         **/
        explicit CSeqMaskerOstatOptAscii( 
                CNcbiOstream & os, Uint2 sz, string const & metadata );

        /**
         **\brief Object destructor.
         **/
        virtual ~CSeqMaskerOstatOptAscii() {}

        /** Get actual counts format version. */
        virtual CSeqMaskerVersion const & GetStatFmtVersion() const { 
            return FormatVersion; 
        }

    protected:

        /**
         **\brief Write optimized data to the output stream.
         **\param p structure containing the hash table and threshold
         **         parameters and pointers to data tables
         **/
        virtual void write_out( const params & p ) const;
};

END_NCBI_SCOPE

#endif
