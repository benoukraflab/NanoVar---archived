#ifndef CORELIB___NCBIMISC__HPP
#define CORELIB___NCBIMISC__HPP

/*  $Id: ncbimisc.hpp 371900 2012-08-13 18:07:13Z vakatov $
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
 * Author:  Denis Vakatov, Eugene Vasilchenko
 *
 *
 */

/// @file ncbimisc.hpp
/// Miscellaneous common-use basic types and functionality


//#include <corelib/ncbistl.hpp>
#ifdef HAVE_SYS_TYPES_H
#  include <sys/types.h>
#endif
#ifdef NCBI_COMPILER_ICC
// Preemptively pull in <cctype>, which breaks if we've already
// repointed is* at NCBI_is*.
#  include <cctype>
#else
#  include <ctype.h>
#endif

#if !defined(HAVE_NULLPTR)
#ifndef nullptr
#  define nullptr NULL
#endif
#endif

#if defined(_DEBUG)  &&  !defined(NCBI_NO_STRICT_CTYPE_ARGS)
#  define NCBI_STRICT_CTYPE_ARGS
#endif


/** @addtogroup AppFramework
 *
 * @{
 */

#include "mask_macros.h"


#ifndef NCBI_ESWITCH_DEFINED
#define NCBI_ESWITCH_DEFINED

extern "C" {

/*
 * ATTENTION!   Do not change this enumeration!
 *
 * It must always be kept in sync with its plain C counterpart defined in
 * "connect/ncbi_types.h". If you absolutely(sic!) need to alter this
 * type, please apply equivalent changes to both definitions.
 */

/** Aux. enum to set/unset/default various features.
 */
typedef enum ENcbiSwitch {
    eOff = 0,
    eOn,
    eDefault
} ESwitch;

} // extern "C"

#endif //!NCBI_ESWITCH_DEFINED


#ifndef NCBI_EOWNERSHIP_DEFINED
#define NCBI_EOWNERSHIP_DEFINED

extern "C" {

/*
 * ATTENTION!   Do not change this enumeration!
 *
 * It must always be kept in sync with its plain C counterpart defined in
 * "connect/ncbi_types.h". If you absolutely(sic!) need to alter this
 * type, please apply equivalent changes to both definitions.
 */

/** Ownership relations between objects.
 *
 * Can be used to define or transfer ownership of objects.
 * For example, specify if a CSocket object owns its underlying SOCK object.
 */
typedef enum ENcbiOwnership {
    eNoOwnership,       /** No ownership assumed                    */
    eTakeOwnership      /** An object can take ownership of another */
} EOwnership;

} // extern "C"

#endif //!NCBI_EOWNERSHIP_DEFINED

BEGIN_NCBI_NAMESPACE;


/// Whether a value is nullable.
enum ENullable {
    eNullable,          ///< Value can be null
    eNotNullable        ///< Value cannot be null
};


/// Signedness of a value.
enum ESign {
    eNegative = -1,     ///< Value is negative
    eZero     =  0,     ///< Value is zero
    ePositive =  1      ///< Value is positive
};


/// Whether to truncate/round a value.
enum ERound {
    eTrunc,             ///< Value must be truncated
    eRound              ///< Value must be rounded
};


/// Whether to follow symbolic links (also known as shortcuts or aliases)
enum EFollowLinks {
    eIgnoreLinks,       ///< Do not follow symbolic links
    eFollowLinks        ///< Follow symbolic links
};


/// Whether to normalize a path
enum ENormalizePath {
    eNormalizePath,     ///< Normalize a path
    eNotNormalizePath   ///< Do not normalize a path
};


/// Interrupt on signal mode
///
/// On UNIX some functions can be interrupted by a signal and EINTR errno
/// value. We can restart or cancel its execution.
enum EInterruptOnSignal {
    eInterruptOnSignal, ///< Cancel operation if interrupted by a signal
    eRestartOnSignal    ///< Restart operation if interrupted by a signal
};


/////////////////////////////////////////////////////////////////////////////
/// Support for safe bool operators
/////////////////////////////////////////////////////////////////////////////


/// Low level macro for declaring safe bool operator.
#define DECLARE_SAFE_BOOL_METHOD(Expr)                                  \
    struct SSafeBoolTag {                                               \
        void SafeBoolTrue(SSafeBoolTag*) {}                             \
    };                                                                  \
    typedef void (SSafeBoolTag::*TBoolType)(SSafeBoolTag*);             \
    operator TBoolType() const {                                        \
        return (Expr)? &SSafeBoolTag::SafeBoolTrue: 0;                  \
    }                                                                   \
    private:                                                            \
    bool operator==(TBoolType) const;                                   \
    bool operator!=(TBoolType) const;                                   \
    public:                                                             \
    friend struct SSafeBoolTag


/// Declaration of safe bool operator from boolean expression.
/// Actual operator declaration will be:
///    operator TBoolType(void) const;
/// where TBoolType is a typedef convertible to bool (member pointer).
#define DECLARE_OPERATOR_BOOL(Expr)             \
    DECLARE_SAFE_BOOL_METHOD(Expr)


/// Declaration of safe bool operator from pointer expression.
/// Actual operator declaration will be:
///    operator bool(void) const;
#define DECLARE_OPERATOR_BOOL_PTR(Ptr)          \
    DECLARE_OPERATOR_BOOL((Ptr) != 0)


/// Declaration of safe bool operator from CRef<>/CConstRef<> expression.
/// Actual operator declaration will be:
///    operator bool(void) const;
#define DECLARE_OPERATOR_BOOL_REF(Ref)          \
    DECLARE_OPERATOR_BOOL((Ref).NotNull())


/// Template used for empty base class optimization.
/// See details in the August '97 "C++ Issue" of Dr. Dobb's Journal
/// Also available from http://www.cantrip.org/emptyopt.html
/// We store usually empty template argument class together with data member.
/// This template is much like STL's pair<>, but the access to members
/// is done though methods first() and second() returning references
/// to corresponding members.
/// First template argument is represented as private base class,
/// while second template argument is represented as private member.
/// In addition to constructor taking two arguments,
/// we add constructor for initialization of only data member (second).
/// This is useful since usually first type is empty and doesn't require
/// non-trivial constructor.
/// We do not define any comparison functions as this template is intented
/// to be used internally within another templates,
/// which themselves should provide any additional functionality.

template<class Base, class Member>
class pair_base_member : private Base
{
public:
    typedef Base base_type;
    typedef Base first_type;
    typedef Member member_type;
    typedef Member second_type;
    
    pair_base_member(void)
        : base_type(), m_Member()
        {
        }
    
    explicit pair_base_member(const member_type& member_value)
        : base_type(), m_Member(member_value)
        {
        }
    
    explicit pair_base_member(const first_type& first_value,
                              const second_type& second_value)
        : base_type(first_value), m_Member(second_value)
        {
        }
    
    const first_type& first() const
        {
            return *this;
        }
    first_type& first()
        {
            return *this;
        }

    const second_type& second() const
        {
            return m_Member;
        }
    second_type& second()
        {
            return m_Member;
        }

    void Swap(pair_base_member<first_type, second_type>& p)
        {
            if (static_cast<void*>(&first()) != static_cast<void*>(&second())) {
                // work around an IBM compiler bug which causes it to perform
                // a spurious 1-byte swap, yielding mixed-up values.
                swap(first(), p.first());
            }
            swap(second(), p.second());
        }

private:
    member_type m_Member;
};


#ifdef HAVE_NO_AUTO_PTR


/////////////////////////////////////////////////////////////////////////////
///
/// auto_ptr --
///
/// Define auto_ptr if needed.
///
/// Replacement of STL's std::auto_ptr for compilers with poor "auto_ptr"
/// implementation.
/// 
/// See C++ Toolkit documentation for limitations and use of auto_ptr.

template <class X>
class auto_ptr
{
    // temporary class for auto_ptr copying
    template<class Y>
    struct auto_ptr_ref
    {
        auto_ptr_ref(auto_ptr<Y>& ptr)
            : m_AutoPtr(ptr)
            {
            }
        auto_ptr<Y>& m_AutoPtr;
    };
public:
    typedef X element_type;         ///< Define element_type

    /// Explicit conversion to auto_ptr.
    explicit auto_ptr(X* p = 0) : m_Ptr(p) {}

    /// Copy constructor with implicit conversion.
    ///
    /// Note that the copy constructor parameter is not a const
    /// because it is modified -- ownership is transferred.
    auto_ptr(auto_ptr<X>& a) : m_Ptr(a.release()) {}

    /// Assignment operator.
    auto_ptr<X>& operator=(auto_ptr<X>& a) {
        if (this != &a) {
            if (m_Ptr  &&  m_Ptr != a.m_Ptr) {
                delete m_Ptr;
            }
            m_Ptr = a.release();
        }
        return *this;
    }

    auto_ptr(auto_ptr_ref<X> ref)
        : m_Ptr(ref.m_AutoPtr.release())
    {
    }
    template <typename Y>
    operator auto_ptr_ref<Y>()
    {
        return auto_ptr_ref<Y>(*this);
    }

    /// Destructor.
    ~auto_ptr(void) {
        if ( m_Ptr )
            delete m_Ptr;
    }

    /// Deference operator.
    X&  operator*(void) const { return *m_Ptr; }

    /// Reference operator.
    X*  operator->(void) const { return m_Ptr; }

    /// Equality operator.
    int operator==(const X* p) const { return (m_Ptr == p); }

    /// Get pointer value.
    X*  get(void) const { return m_Ptr; }

    /// Release pointer.
    X* release(void) {
        X* x_Ptr = m_Ptr;  m_Ptr = 0;  return x_Ptr;
    }

    /// Reset pointer.
    void reset(X* p = 0) {
        if (m_Ptr != p) {
            delete m_Ptr;
            m_Ptr = p;
        }
    }

private:
    X* m_Ptr;               ///< Internal pointer implementation.
};

#endif /* HAVE_NO_AUTO_PTR */



/// Functor template for allocating object.
template<class X>
struct Creater
{
    /// Default create function.
    static X* Create(void)
    { return new X; }
};

/// Functor template for deleting object.
template<class X>
struct Deleter
{
    /// Default delete function.
    static void Delete(X* object)
    { delete object; }
};

/// Functor template for deleting array of objects.
template<class X>
struct ArrayDeleter
{
    /// Array delete function.
    static void Delete(X* object)
    { delete[] object; }
};

/// Functor template for the C language deallocation function, free().
template<class X>
struct CDeleter
{
    /// C Language deallocation function.
    static void Delete(X* object)
    { free(object); }
};



/////////////////////////////////////////////////////////////////////////////
///
/// AutoPtr --
///
/// Define an "auto_ptr" like class that can be used inside STL containers.
///
/// The Standard auto_ptr template from STL doesn't allow the auto_ptr to be
/// put in STL containers (list, vector, map etc.). The reason for this is
/// the absence of copy constructor and assignment operator.
/// We decided that it would be useful to have an analog of STL's auto_ptr
/// without this restriction - AutoPtr.
///
/// Due to nature of AutoPtr its copy constructor and assignment operator
/// modify the state of the source AutoPtr object as it transfers ownership
/// to the target AutoPtr object. Also, we added possibility to redefine the
/// way pointer will be deleted: the second argument of template allows
/// pointers from "malloc" in AutoPtr, or you can use "ArrayDeleter" (see
/// above) to properly delete an array of objects using "delete[]" instead
/// of "delete". By default, the internal pointer will be deleted by C++
/// "delete" operator.
///
/// @sa
///   Deleter(), ArrayDeleter(), CDeleter()

template< class X, class Del = Deleter<X> >
class AutoPtr
{
public:
    typedef X element_type;         ///< Define element type.
    typedef Del deleter_type;       ///< Alias for template argument.

    /// Constructor.
    AutoPtr(element_type* p = 0)
        : m_Ptr(p), m_Data(true)
    {
    }

    /// Constructor.
    AutoPtr(element_type* p, const deleter_type& deleter)
        : m_Ptr(p), m_Data(deleter, true)
    {
    }

    /// Constructor, own the pointed object if ownership == eTakeOwnership
    AutoPtr(element_type* p, EOwnership ownership)
        : m_Ptr(p), m_Data(ownership != eNoOwnership)
    {
    }

    /// Constructor, own the pointed object if ownership == eTakeOwnership
    AutoPtr(element_type* p, const deleter_type& deleter, EOwnership ownership)
        : m_Ptr(p), m_Data(deleter, ownership != eNoOwnership)
    {
    }

    /// Copy constructor.
    AutoPtr(const AutoPtr<X, Del>& p)
        : m_Ptr(0), m_Data(p.m_Data)
    {
        m_Ptr = p.x_Release();
    }

    /// Destructor.
    ~AutoPtr(void)
    {
        reset();
    }

    /// Assignment operator.
    AutoPtr<X, Del>& operator=(const AutoPtr<X, Del>& p)
    {
        if (this != &p) {
            bool owner = p.m_Data.second();
            reset(p.x_Release());
            m_Data.second() = owner;
        }
        return *this;
    }

    /// Assignment operator.
    AutoPtr<X, Del>& operator=(element_type* p)
    {
        reset(p);
        return *this;
    }

    /// Bool operator for use in if() clause.
    DECLARE_OPERATOR_BOOL_PTR(m_Ptr);

    // Standard getters.

    /// Dereference operator.
    element_type& operator* (void) const { return *m_Ptr; }

    /// Reference operator.
    element_type* operator->(void) const { return  m_Ptr; }

    /// Get pointer.
    element_type* get       (void) const { return  m_Ptr; }

    /// Release will release ownership of pointer to caller.
    element_type* release(void)
    {
        m_Data.second() = false;
        return m_Ptr;
    }

    /// Reset will delete old pointer, set content to new value,
    /// and accept ownership upon the new pointer.
    void reset(element_type* p = 0, EOwnership ownership = eTakeOwnership)
    {
        if ( m_Ptr != p ) {
            if (m_Ptr  &&  m_Data.second()) {
                m_Data.first().Delete(release());
            }
            m_Ptr   = p;
        }
        m_Data.second() = p != 0 && ownership == eTakeOwnership;
    }

    void Swap(AutoPtr<X, Del>& a)
    {
        swap(m_Ptr, a.m_Ptr);
        swap(m_Data, a.m_Data);
    }

private:
    element_type* m_Ptr;                  ///< Internal pointer representation.
    mutable pair_base_member<deleter_type, bool> m_Data; ///< State info.

    /// Release for const object.
    element_type* x_Release(void) const
    {
        return const_cast<AutoPtr<X, Del>*>(this)->release();
    }
};


/////////////////////////////////////////////////////////////////////////////
///
/// AutoArray --
///
/// "AutoPtr" like class for using with arrays
///
/// vector<> template comes with a performance penalty, since it always
/// initializes its content. This template is not a vector replacement,
/// it's a version of AutoPtr<> tuned for array pointers. For convenience
/// it defines array style access operator [] and size based contructor.
///
/// @sa AutoPtr
///

template< class X, class Del = ArrayDeleter<X> >
class AutoArray
{
public:
    typedef X element_type;         ///< Define element type.
    typedef Del deleter_type;       ///< Alias for template argument.

public:

    /// Construct the array using C++ new[] operator
    /// @note In this case you should use ArrayDeleter<> or compatible
    explicit AutoArray(size_t size)
        : m_Ptr(new element_type[size]), m_Data(true)
    {}

    explicit AutoArray(element_type* p = 0)
        : m_Ptr(p), m_Data(true)
    {}

    AutoArray(element_type* p, const deleter_type& deleter)
        : m_Ptr(p), m_Data(deleter, true)
    {
    }

    AutoArray(const AutoArray<X, Del>& p)
        : m_Ptr(0), m_Data(p.m_Data)
    {
        m_Ptr = p.x_Release();
    }

    ~AutoArray(void)
    {
        reset();
    }

    /// Assignment operator.
    AutoArray<X, Del>& operator=(const AutoArray<X, Del>& p)
    {
        if (this != &p) {
            bool owner = p.m_Data.second();
            reset(p.x_Release());
            m_Data.second() = owner;
        }
        return *this;
    }

    /// Assignment operator.
    AutoArray<X, Del>& operator=(element_type* p)
    {
        reset(p);
        return *this;
    }
    /// Bool operator for use in if() clause.
    DECLARE_OPERATOR_BOOL_PTR(m_Ptr);

    /// Get pointer.
    element_type* get       (void) const { return  m_Ptr; }

    /// Release will release ownership of pointer to caller.
    element_type* release(void)
    {
        m_Data.second() = false;
        return m_Ptr;
    }

    /// array style dereference (returns value)
    const element_type& operator[](size_t pos) const { return m_Ptr[pos]; }

    /// array style dereference (returns reference)
    element_type& operator[](size_t pos) { return m_Ptr[pos]; }

    /// Reset will delete old pointer, set content to new value,
    /// and accept ownership upon the new pointer.
    void reset(element_type* p = 0)
    {
        if (m_Ptr  &&  m_Data.second()) {
            m_Data.first().Delete(release());
        }
        m_Ptr   = p;
        m_Data.second() = true;
    }

    void Swap(AutoPtr<X, Del>& a)
    {
        swap(m_Ptr, a.m_Ptr);
        swap(m_Data, a.m_Data);
    }

private:
    /// Release for const object.
    element_type* x_Release(void) const
    {
        return const_cast<AutoArray<X, Del>*>(this)->release();
    }

private:
    element_type*  m_Ptr;
    mutable pair_base_member<deleter_type, bool> m_Data; ///< State info.
};



// "min" and "max" templates
//

// Always get rid of the old non-conformant min/max macros
#ifdef min
#  undef min
#endif
#ifdef max
#  undef max
#endif

#if defined(HAVE_NO_MINMAX_TEMPLATE)
#  define NOMINMAX

/// Min function template.
template <class T>
inline
const T& min(const T& a, const T& b) {
    return b < a ? b : a;
}

/// Max function template.
template <class T>
inline
const T& max(const T& a, const T& b) {
    return  a < b ? b : a;
}
#endif /* HAVE_NO_MINMAX_TEMPLATE */



// strdup()
//

//#ifndef HAVE_STRDUP
/// Supply string duplicate function, if one is not defined.
//extern char* strdup(const char* str);
//#endif



// ctype hacks
//

#ifdef NCBI_STRICT_CTYPE_ARGS

END_NCBI_NAMESPACE;

#define NCBI_CTYPEFAKEBODY \
  { return See_the_standard_on_proper_argument_type_for_ctype_macros(c); }

#ifdef isalpha
inline int NCBI_isalpha(unsigned char c) { return isalpha(c); }
inline int NCBI_isalpha(int           c) { return isalpha(c); }
template<class C>
inline int NCBI_isalpha(C c) NCBI_CTYPEFAKEBODY
#undef  isalpha
#define isalpha NCBI_isalpha
#endif

#ifdef isalnum
inline int NCBI_isalnum(unsigned char c) { return isalnum(c); }
inline int NCBI_isalnum(int           c) { return isalnum(c); }
template<class C>
inline int NCBI_isalnum(C c) NCBI_CTYPEFAKEBODY
#undef  isalnum
#define isalnum NCBI_isalnum
#endif

#ifdef isascii
inline int NCBI_isascii(unsigned char c) { return isascii(c); }
inline int NCBI_isascii(int           c) { return isascii(c); }
template<class C>
inline int NCBI_isascii(C c) NCBI_CTYPEFAKEBODY
#undef  isascii
#define isascii NCBI_isascii
#endif

#ifdef isblank
inline int NCBI_isblank(unsigned char c) { return isblank(c); }
inline int NCBI_isblank(int           c) { return isblank(c); }
template<class C>
inline int NCBI_isblank(C c) NCBI_CTYPEFAKEBODY
#undef  isblank
#define isblank NCBI_isblank
#endif

#ifdef iscntrl
inline int NCBI_iscntrl(unsigned char c) { return iscntrl(c); }
inline int NCBI_iscntrl(int           c) { return iscntrl(c); }
template<class C>
inline int NCBI_iscntrl(C c) NCBI_CTYPEFAKEBODY
#undef  iscntrl
#define iscntrl NCBI_iscntrl
#endif

#ifdef isdigit
inline int NCBI_isdigit(unsigned char c) { return isdigit(c); }
inline int NCBI_isdigit(int           c) { return isdigit(c); }
template<class C>
inline int NCBI_isdigit(C c) NCBI_CTYPEFAKEBODY
#undef  isdigit
#define isdigit NCBI_isdigit
#endif

#ifdef isgraph
inline int NCBI_isgraph(unsigned char c) { return isgraph(c); }
inline int NCBI_isgraph(int           c) { return isgraph(c); }
template<class C>
inline int NCBI_isgraph(C c) NCBI_CTYPEFAKEBODY
#undef  isgraph
#define isgraph NCBI_isgraph
#endif

#ifdef islower
inline int NCBI_islower(unsigned char c) { return islower(c); }
inline int NCBI_islower(int           c) { return islower(c); }
template<class C>
inline int NCBI_islower(C c) NCBI_CTYPEFAKEBODY
#undef  islower
#define islower NCBI_islower
#endif

#ifdef isprint
inline int NCBI_isprint(unsigned char c) { return isprint(c); }
inline int NCBI_isprint(int           c) { return isprint(c); }
template<class C>
inline int NCBI_isprint(C c) NCBI_CTYPEFAKEBODY
#undef  isprint
#define isprint NCBI_isprint
#endif

#ifdef ispunct
inline int NCBI_ispunct(unsigned char c) { return ispunct(c); }
inline int NCBI_ispunct(int           c) { return ispunct(c); }
template<class C>
inline int NCBI_ispunct(C c) NCBI_CTYPEFAKEBODY
#undef  ispunct
#define ispunct NCBI_ispunct
#endif

#ifdef isspace
inline int NCBI_isspace(unsigned char c) { return isspace(c); }
inline int NCBI_isspace(int           c) { return isspace(c); }
template<class C>
inline int NCBI_isspace(C c) NCBI_CTYPEFAKEBODY
#undef  isspace
#define isspace NCBI_isspace
#endif

#ifdef isupper
inline int NCBI_isupper(unsigned char c) { return isupper(c); }
inline int NCBI_isupper(int           c) { return isupper(c); }
template<class C>
inline int NCBI_isupper(C c) NCBI_CTYPEFAKEBODY
#undef  isupper
#define isupper NCBI_isupper
#endif

#ifdef isxdigit
inline int NCBI_isxdigit(unsigned char c) { return isxdigit(c); }
inline int NCBI_isxdigit(int           c) { return isxdigit(c); }
template<class C>
inline int NCBI_isxdigit(C c) NCBI_CTYPEFAKEBODY
#undef  isxdigit
#define isxdigit NCBI_isxdigit
#endif

#ifdef toascii
inline int NCBI_toascii(unsigned char c) { return toascii(c); }
inline int NCBI_toascii(int           c) { return toascii(c); }
template<class C>
inline int NCBI_toascii(C c) NCBI_CTYPEFAKEBODY
#undef  toascii
#define toascii NCBI_toascii
#endif

#ifdef tolower
inline int NCBI_tolower(unsigned char c) { return tolower(c); }
inline int NCBI_tolower(int           c) { return tolower(c); }
template<class C>
inline int NCBI_tolower(C c) NCBI_CTYPEFAKEBODY
#undef  tolower
#define tolower NCBI_tolower
#endif

#ifdef toupper
inline int NCBI_toupper(unsigned char c) { return toupper(c); }
inline int NCBI_toupper(int           c) { return toupper(c); }
template<class C>
inline int NCBI_toupper(C c) NCBI_CTYPEFAKEBODY
#undef  toupper
#define toupper NCBI_toupper
#endif

#undef NCBI_CTYPEFAKEBODY

BEGIN_NCBI_NAMESPACE;

#endif // NCBI_STRICT_CTYPE_ARGS

//  ITERATE
//  NON_CONST_ITERATE
//  ERASE_ITERATE
//
// Useful macro to write 'for' statements with the STL container iterator as
// a variable.
//

/// ITERATE macro to sequence through container elements.
#define ITERATE(Type, Var, Cont) \
    for ( Type::const_iterator Var = (Cont).begin(), NCBI_NAME2(Var,_end) = (Cont).end();  Var != NCBI_NAME2(Var,_end);  ++Var )

/// Non constant version of ITERATE macro.
#define NON_CONST_ITERATE(Type, Var, Cont) \
    for ( Type::iterator Var = (Cont).begin();  Var != (Cont).end();  ++Var )

/// NON_CONST_ITERATE macro for STL sets.  A special macro is required because 
/// STL sets iterate in a const manner even if you use their "iterator" type.
/// Using this macro is discouraged because if you affect an object's state
/// in such a way that it's ordering in the set should change, the set
/// could be in an invalid state.
/// The "REAL92137892" is just a random suffix to make the variable unique.
/// The number has no significance.
#define NON_CONST_SET_ITERATE( Type, Var, Cont) \
    Type::key_type *Var = NULL;  \
    Type::iterator Var##REAL92137892  = (Cont).begin(); \
    Var = ( Var##REAL92137892 != (Cont).end() ? (Type::key_type *)&*Var##REAL92137892 : NULL ); \
    for ( ; Var##REAL92137892 != (Cont).end(); \
    ( ++Var##REAL92137892 != (Cont).end() ? Var = (Type::key_type *)&*(Var##REAL92137892) : Var = NULL ) )

/// Non constant version with ability to erase current element, if container permits.
/// Use only on containers, for which erase do not ruin other iterators into the container
/// e.g., map, list, but NOT vector.
/// See also VECTOR_ERASE
#define ERASE_ITERATE(Type, Var, Cont) \
    for ( Type::iterator Var, NCBI_NAME2(Var,_next) = (Cont).begin();   \
          (Var = NCBI_NAME2(Var,_next)) != (Cont).end() &&              \
              (++NCBI_NAME2(Var,_next), true); )
/// Use this macro inside body of ERASE_ITERATE cycle to erase from
/// vector-like container. Plain erase() call would invalidate Var_next
/// iterator and would make the cycle controlling code to fail.
#define VECTOR_ERASE(Var, Cont) (NCBI_NAME2(Var,_next) = (Cont).erase(Var))

/// ITERATE macro to reverse sequence through container elements.
#define REVERSE_ITERATE(Type, Var, Cont) \
    for ( Type::const_reverse_iterator Var = (Cont).rbegin(), NCBI_NAME2(Var,_end) = (Cont).rend();  Var != NCBI_NAME2(Var,_end);  ++Var )

/// Non constant version of REVERSE_ITERATE macro.
#define NON_CONST_REVERSE_ITERATE(Type, Var, Cont) \
    for ( Type::reverse_iterator Var = (Cont).rbegin();  Var != (Cont).rend();  ++Var )


/// Type for sequence locations and lengths.
///
/// Use this typedef rather than its expansion, which may change.
typedef unsigned int TSeqPos;

/// Define special value for invalid sequence position.
const TSeqPos kInvalidSeqPos = ((TSeqPos) (-1));


/// Type for signed sequence position.
///
/// Use this type when and only when negative values are a possibility
/// for reporting differences between positions, or for error reporting --
/// though exceptions are generally better for error reporting.
/// Use this typedef rather than its expansion, which may change.
typedef int TSignedSeqPos;

/// Type for sequence GI.
///
/// Use this typedef rather than its expansion, which may change.
typedef int TGi;


/// Helper address class
class CRawPointer
{
public:
    /// add offset to object reference (to get object's member)
    static void* Add(void* object, ssize_t offset);
    static const void* Add(const void* object, ssize_t offset);
    /// calculate offset inside object
    static ssize_t Sub(const void* first, const void* second);
};


inline
void* CRawPointer::Add(void* object, ssize_t offset)
{
    return static_cast<char*> (object) + offset;
}

inline
const void* CRawPointer::Add(const void* object, ssize_t offset)
{
    return static_cast<const char*> (object) + offset;
}

inline
ssize_t CRawPointer::Sub(const void* first, const void* second)
{
    return (ssize_t)/*ptrdiff_t*/
        (static_cast<const char*> (first) - static_cast<const char*> (second));
}



/// Buffer with an embedded pre-reserved space.
///
/// It is convenient to use if you want to avoid allocation in heap
/// for smaller requested sizes, and use stack for such cases.
/// @example:
///    CFastBuffer<2048> buf(some_size);

template <size_t KEmbeddedSize, class TType = char>
class CFastBuffer
{
public:
    CFastBuffer(size_t size)
        : m_Size(size),
          m_Buffer(size <= KEmbeddedSize ? m_EmbeddedBuffer : new TType[size])
    {}
    ~CFastBuffer() { if (m_Buffer != m_EmbeddedBuffer) delete[] m_Buffer; }

    TType        operator[] (size_t pos) const { return m_Buffer[pos]; }

    TType&       operator* (void)       { return *m_Buffer; }
    const TType& operator* (void) const { return *m_Buffer; }

    TType*       operator+ (size_t offset)       { return m_Buffer + offset; }
    const TType* operator+ (size_t offset) const { return m_Buffer + offset; }

    TType*       begin(void)       { return m_Buffer; }
    const TType* begin(void) const { return m_Buffer; }

    TType*       end(void)       { return m_Buffer + m_Size; }
    const TType* end(void) const { return m_Buffer + m_Size; }

    size_t       size(void) const { return m_Size; }

private:
     size_t m_Size;
     TType* m_Buffer;
     TType  m_EmbeddedBuffer[KEmbeddedSize];
};



/// Macro used to mark a constructor as deprecated.
///
/// The correct syntax for this varies from compiler to compiler:
/// older versions of GCC (prior to 3.4) require NCBI_DEPRECATED to
/// follow any relevant constructor declarations, but some other
/// compilers (Microsoft Visual Studio 2005, IBM Visual Age / XL)
/// require it to precede any relevant declarations, whether or not
/// they are for constructors.
#if defined(NCBI_COMPILER_MSVC) || defined(NCBI_COMPILER_VISUALAGE)
#  define NCBI_DEPRECATED_CTOR(decl) NCBI_DEPRECATED decl
#else
#  define NCBI_DEPRECATED_CTOR(decl) decl NCBI_DEPRECATED
#endif

/// Macro used to mark a class as deprecated.
///
/// @sa NCBI_DEPRECATED_CTOR
#define NCBI_DEPRECATED_CLASS NCBI_DEPRECATED_CTOR(class)

END_NCBI_NAMESPACE;

BEGIN_STD_NAMESPACE;

template<class T1, class T2>
inline
void swap(NCBI_NS_NCBI::pair_base_member<T1,T2>& pair1,
          NCBI_NS_NCBI::pair_base_member<T1,T2>& pair2)
{
    pair1.Swap(pair2);
}


template<class P, class D>
inline
void swap(NCBI_NS_NCBI::AutoPtr<P,D>& ptr1,
          NCBI_NS_NCBI::AutoPtr<P,D>& ptr2)
{
    ptr1.Swap(ptr2);
}


#if (defined(NCBI_COMPILER_GCC) && NCBI_COMPILER_VERSION < 340)  ||  defined(NCBI_COMPILER_WORKSHOP)  ||  defined(NCBI_COMPILER_MIPSPRO)

#define ArraySize(array)  (sizeof(array)/sizeof((array)[0]))

#else

template<class Element, size_t Size>
inline
size_t ArraySize(const Element (&)[Size])
{
    return Size;
}

#endif

END_STD_NAMESPACE;

/// Definition of packed enum type, to save some memory
#if defined(NCBI_COMPILER_MSVC)
#  define NCBI_PACKED_ENUM_TYPE(type)  : type
#  define NCBI_PACKED_ENUM_END()
#elif (defined(NCBI_COMPILER_GCC) && NCBI_COMPILER_VERSION >= 400)  ||  defined(NCBI_COMPILER_ICC)
#  define NCBI_PACKED_ENUM_TYPE(type)
#  define NCBI_PACKED_ENUM_END()        __attribute__((packed))
#else
#  define NCBI_PACKED_ENUM_TYPE(type)
#  define NCBI_PACKED_ENUM_END()
#endif


/* @} */

#endif  /* CORELIB___NCBIMISC__HPP */
