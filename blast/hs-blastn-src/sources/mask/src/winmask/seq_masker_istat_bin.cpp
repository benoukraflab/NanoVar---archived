/*  $Id: seq_masker_istat_bin.cpp 462550 2015-03-19 14:07:19Z morgulis $
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
 *   Definition for CSeqMaskerIstatBin class.
 *
 */

#include "ncbi_pch.hpp"

#include "seq_masker_istat_bin.hpp"

BEGIN_NCBI_SCOPE

static const streamsize HEADER_LEN  = 2*sizeof( Uint4 );
static const streamsize TRAILER_LEN = 4*sizeof( Uint4 );

//------------------------------------------------------------------------------
//const char * CSeqMaskerIstatBin::Exception::GetErrCodeString() const
//{
//    switch( GetErrCode() )
//    {
//        case eStreamOpenFail:   return "open failed";
//        case eFormat:           return "file format error";
//        default:                return CException::GetErrCodeString();
//    }
//}

//------------------------------------------------------------------------------
CSeqMaskerIstatBin::CSeqMaskerIstatBin( const string & name,
                                        Uint4 arg_threshold,
                                        Uint4 arg_textend,
                                        Uint4 arg_max_count,
                                        Uint4 arg_use_max_count,
                                        Uint4 arg_min_count,
                                        Uint4 arg_use_min_count,
                                        Uint4 skip )
    : CSeqMaskerIstat(  arg_threshold, arg_textend, 
                        arg_max_count, arg_use_max_count,
                        arg_min_count, arg_use_min_count )
{
    streamsize iflen = 0;

    //{
    //    CFile input_file( name );

    //    if( !input_file.Exists() )
    //        NCBI_THROW( Exception, eStreamOpenFail, name + " does not exist" );

    //    iflen = (streamsize)input_file.GetLength();

    //    if( iflen < HEADER_LEN + TRAILER_LEN )
    //        NCBI_THROW( Exception, eFormat, "wrong file size" );
    //}
	
	{
		FILE* input_file = fopen( name.c_str(), "r");
		if ( !input_file )
		{
			std::ostringstream err;
			err << "file " << name << " does not exist.";
			error_and_exit(err.str());
		}
		fseek(input_file, 0ull, SEEK_END);
		iflen = ftell(input_file);
		fclose(input_file);
		if (iflen < HEADER_LEN + TRAILER_LEN)
		{
			std::string err = "wrong file size.";
			error_and_exit(err);
		}
	}

    CNcbiIfstream in_stream( name.c_str(), std::ios::binary); //IOS_BASE::binary );

    {
        char * data( new char[skip] );
        in_stream.read( data, skip );
        delete[] data;
    }

    streamsize adj( in_stream.tellg() );
    Uint4 data;
    in_stream.read( (char *)&data, sizeof( Uint4 ) );

    {
        in_stream.read( (char *)&data, sizeof( Uint4 ) );
        Uint1 us = (Uint1)data;
        
        if( us == 0 || us > 16 )
        {
			//NCBI_THROW( Exception, eFormat, "illegal unit size" );
			std::string err = "illegal unit size.";
			error_and_exit(err);
		}

        uset.set_unit_size( us );
    }

    {
        streamsize datalen = iflen - HEADER_LEN - TRAILER_LEN - adj;

        if( datalen%(2*sizeof( Uint4 )) != 0 )
		{
			//NCBI_THROW( Exception, eFormat, "wrong length" );
			std::string err = "wrong length.";
			error_and_exit(err);
		}

        Uint4 * cdata = 0;

        if( datalen > 0 )
        {
            cdata = new Uint4[datalen/sizeof( Uint4 )];
            in_stream.read( (char *)cdata, datalen );
            uset.add_info( cdata, datalen/sizeof( Uint4 ) );
        }
    }

    in_stream.read( (char *)&data, sizeof( Uint4 ) );

    set_min_count( data );

    in_stream.read( (char *)&data, sizeof( Uint4 ) );

    if( get_textend() == 0 )
        set_textend( data );

    in_stream.read( (char *)&data, sizeof( Uint4 ) );

    if( get_threshold() == 0 )
        set_threshold( data );

    in_stream.read( (char *)&data, sizeof( Uint4 ) );

    if( get_max_count() == 0 )
        set_max_count( data );

    if( get_use_min_count() == 0 )
      set_use_min_count( (get_min_count() + 1)/2 );

    if( get_use_max_count() == 0 )
      set_use_max_count( get_max_count() );
}

//------------------------------------------------------------------------------
Uint4 CSeqMaskerIstatBin::trueat( Uint4 unit ) const
{ return uset.get_info( unit ); }

//------------------------------------------------------------------------------
Uint4 CSeqMaskerIstatBin::at( Uint4 unit ) const
{
    Uint4 res = uset.get_info( unit );

    if( res == 0 || res < get_min_count() )
        return get_use_min_count();

    return (res > get_max_count()) ? get_use_max_count() : res;
}

END_NCBI_SCOPE
