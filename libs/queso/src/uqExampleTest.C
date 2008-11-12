#define BOOST_TEST_MODULE $Id: uqExampleTest.C $
#include <boost/test/included/unit_test.hpp>

int add( int i, int j ) { return i+j; }
int mult( int i, int j ) { return i*j; }

BOOST_AUTO_TEST_CASE( simple_test_add )
{
    BOOST_CHECK(   add( 2,2 ) == 4 );  // continues on error
    BOOST_REQUIRE( add( 2,2 ) == 4 );  // throws on error
}

BOOST_AUTO_TEST_CASE( simple_test_mult )
{
    BOOST_CHECK(   mult( 2,2 ) == 4 );  // continues on error
    BOOST_REQUIRE( mult( 2,2 ) == 4 );  // throws on error
}
