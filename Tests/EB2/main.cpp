#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include "MyTest.H"

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    MyTest mytest;

    mytest.write();

    amrex::Finalize();
}
