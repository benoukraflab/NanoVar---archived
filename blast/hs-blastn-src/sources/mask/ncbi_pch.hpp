#if defined(NCBI_USE_PCH)  &&  !defined(NCBI_PCH__HPP)
/*  $Id: ncbi_pch.hpp 391319 2013-03-06 22:38:49Z camacho $
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
 */

/** @file ncbi_pch.hpp
 ** Header file to be pre-compiled and speed up build of NCBI C++ Toolkit
 ** (wrapper for the "real" version under common/)
 **/
#  define NCBI_PCH__HPP
#  include <common/ncbi_pch_impl.hpp>
#endif
