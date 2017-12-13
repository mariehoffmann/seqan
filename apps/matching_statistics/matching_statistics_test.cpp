// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

#include <chrono>
#include <experimental/filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

//#include <seqan3/alphabet/nucleotide/dna4.hpp>
//#include <seqan3/range/view/to_char.hpp>

#include <seqan/index.h>
#include <unidirectional_matching_statistics.h>
//#include <seqan3/indexes/matching_statistics/unidirectional_matching_statistics.hpp>

#include <sdsl/config.hpp> // for cache_config
#include <sdsl/construct.hpp>
#include <sdsl/construct_bwt.hpp>
#include <sdsl/suffix_trees.hpp>


// forward declaration
//template<typename T> LinkedList<T>::LinkedList()
using namespace seqan;

using container_t = typename std::vector<dna4>;

using test_msg_start
//template <typename container_type> struct seqan3::MS;
namespace fs = std::experimental::filesystem;

struct SetUp
{
    std::string s = "AACG";
    std::string t = "AACT";
    typedef sdsl::cst_sada<> cst_t;
};

void print_status(std::string function_name, bool start=true)
{
    if (start == true)
        printf("run test %s ...", function_name);
    else
        printf(" successful!\n");
}

// unidirectional matching statistics: default constructor
void default_construction)
{
    print_status("default_construction");
    // default constructors
    MS<container_t> ms{};
    // copy constructor
    MS<container_t> ms2{ms};
    // assignment construction
    MS<container_t> ms3 = ms2;
    print_status("default_construction", false);
}

/*
// unidirectional matching statistics: non-default constructors
TEST_F(matching_statistics_test_fixture, construct_by_sequence)
{
    // construct by sequence and create tmp files where sequences are written to.
    MS<std::string> ms(s, t);
    auto paths = ms.get_source_files();
    EXPECT_TRUE(fs::exists(paths.first) && fs::exists(paths.second));
    std::string line;
    std::ifstream file1(paths.first);
    if (file1.is_open()) {getline(file1,line); file1.close();}
    EXPECT_EQ(s, line);
    std::ifstream file2(paths.second);
    if (file2.is_open()) {getline(file2,line); file2.close();}
    EXPECT_EQ(t, line);
}

// construct by file path
TEST_F(matching_statistics_test_fixture, construct_by_file)
{
    // create files
    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::system_clock::now().time_since_epoch()).count();
    fs::path tmp_dir = fs::current_path() / fs::path{"tmp"} / fs::path{std::to_string(now)};
    fs::create_directory(tmp_dir);
    fs::path file1 = fs::path{"file1.txt"};
    fs::path file2 = fs::path{"file2.txt"};
    std::ofstream(tmp_dir / file1) << s;
    std::ofstream(tmp_dir / file2) << t;
    // construct from file
    MS<container_t> ms(tmp_dir, file1, file2);
    auto paths = ms.get_source_files();
    EXPECT_EQ(tmp_dir / file1, paths.first);
    EXPECT_EQ(tmp_dir / file2, paths.second);
}

// unidirectional matching statistics: sdsl compressed suffix tree construction
TEST_F(matching_statistics_test_fixture, sdsl_cst)
{
    MS<std::string> ms(s, t);
    cst_t cst_s;
    ms.construct_cst(0, cst_s);
    auto root = cst_s.root();
    std::cout << "cst.size = " << cst_s.size() << std::endl;
    std::cout << "num children of root: " << cst_s.children(root).size() << std::endl;
    std::cout << "degree root node: " << cst_s.degree(0) << std::endl;
    EXPECT_EQ(4u, cst_s.degree(0));
}

// construct bwt
TEST_F(matching_statistics_test_fixture, sdsl_bwt)
{
    MS<std::string> ms(s, t);
    sdsl::int_vector<8> bwt_vec;
    EXPECT_TRUE(ms.construct_bwt(0, bwt_vec));
}
*/

int main(/*int argc, char** argv*/)
{
    default_construction();
    return 0;
}
