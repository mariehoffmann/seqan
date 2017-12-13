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

/*!\file
 * \brief Matching statistics (e.g. shortest unique substrings) for sequence comparison.
 * \ingroup container
 * \author Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
 */

#pragma once
//#ifndef __MS_HPP
//#define __MS_HPP

#include <algorithm>
#include <cassert>
#include <experimental/filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <sys/types.h>
#include <type_traits>
#include <utility>
#include <vector>

#include <range/v3/utility/iterator_traits.hpp>
#include <range/v3/view/reverse.hpp>

#include <sdsl/config.hpp> // for cache_config
#include <sdsl/construct.hpp>
#include <sdsl/suffix_trees.hpp>
//#include <sdsl/util.hpp>

//#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
//#include <seqan3/range/concept.hpp>
#include <seqan3/range/view/to_char.hpp>

// Note: in seqan3/sdsl-lite/external cmake CMakeLists.txt has to be executed once to create e.g. needed header file divsufsort.h
namespace fs = std::experimental::filesystem;

typedef std::map<std::string, std::string> tMSS;

using namespace seqan3;
using namespace seqan3::literal;

#define STRING_S 0
#define STRING_T 1

//template <typename alphabet_type>
//auto const seqan3::view::char_to<alphabet_type>;

// [1] Belazzougui et al. "Indexed Matching Statistics and Shortest Unique Substrings"

// TODO: requires that container's value type fits into a byte! like uint8_t in int_types.hpp
template <typename container_t=std::string>
struct MS
{
private:

protected:
    using index_t = signed int; //ranges::v3::size_type_t<container_type>;
    // Algorithm specific value type
    using value_t = unsigned int;
    // Type of container values
    using alphabet_t = ranges::v3::value_type_t<container_t>;

    //!\brief Threshold for number of substring occurrences.
    unsigned short int tau{1};

    sdsl::bit_vector ms;
    sdsl::select_support_mcl<1,1> ss = sdsl::select_support_mcl<1,1>(&ms);

    //!\brief Type of compressed suffix tree.
    typedef sdsl::cst_sada<> cst_t;
    //!\brief Configuration object for sdsl.
    sdsl::cache_config test_config;

    struct source_t
    {
        fs::path filename;
        fs::path filename_rev = "";
        typename std::add_pointer_t<container_t> sequence{nullptr};
    };

    // tmp file directory
    fs::path tmp_dir = fs::current_path() / fs::path{"./tmp"};
    // default data source file path (when creating files from input sequences) and no file paths are given.
    // Else it will be overwritten by construction by file.
    fs::path source_path = tmp_dir;

    std::array<source_t, 2> srcs = {source_t{fs::path{}}, source_t{fs::path{}}};

public:
    //cst_t cst;  // TODO: move to private section, derive test class from this for access
    // TODO: use alphabet vectors and crate strings from them

    std::pair<container_t, container_t> get_sequences() const
    {
        return std::make_pair(*srcs[0]->sequence, *srcs[1]->sequence);
    }

    fs::path get_base_path() const
    {
        return source_path;
    }

    // return absolute paths of source files
    std::pair<fs::path, fs::path> get_source_files() const
    {
        return std::make_pair(source_path / srcs[0].filename, source_path / srcs[1].filename);
    }

    void write_files()
    {
        assert(srcs[0].sequence->size() > 0 && srcs[1].sequence->size() > 0);
        // create reverse strings
//        srcs[STRING_S_REV].sequence = &std::string(srcs[STRING_S].sequence->rbegin(), srcs[STRING_S].sequence->rend());
        // create directory if empty
        if (!fs::exists(tmp_dir))
            fs::create_directory(tmp_dir);
        // create temporary files if not existent
        if (srcs[0].filename.empty() || !srcs[0].filename.has_extension()){
            srcs[0].filename = "s.txt";
            srcs[0].filename_rev = "s_rev.txt";
        }
        if (srcs[1].filename.empty() || !srcs[1].filename.has_extension())
        {
            srcs[1].filename = "t.txt";
            srcs[1].filename_rev = "t_rev.txt";
        }
        // write cached strings and their reverse into files
        for (unsigned short int i = 0; i < srcs.size(); ++i)
        {
            std::ofstream(source_path / srcs[i].filename) << *(srcs[i].sequence);
            std::cout << "get_absolute_path(srcs[i].filename) exists: " << fs::exists(source_path / srcs[i].filename) << std::endl;
            // TODO: use view like auto v = ssss | view::char_to<dna4>;
            std::string tmp(srcs[i].sequence->rbegin(), srcs[i].sequence->rend());
            std::ofstream(source_path / srcs[i].filename_rev) << tmp;
            std::cout << "get_absolute_path(srcs[i].filename_rev) exists: " << fs::exists(source_path / srcs[i].filename_rev) << std::endl;
        }
    }

    //!\brief construct compressed suffix tree of string s according to Sadakane
    void construct_cst(size_t const i, cst_t & cst, bool const reverse=false)
    {
        fs::path src_file = (reverse ? srcs[i].filename_rev : srcs[i].filename);
        assert(fs::exists(source_path / src_file));
        sdsl::construct(cst, source_path / src_file, 1);
        auto root = cst.root();
        for (auto child: cst.children(root)) {
            std::cout << "sada id = " << cst.id(child) << std::endl;
            if (cst.id(child) > 0)
                std::cout << "1st char of edge = '" << cst.edge(child, 1) << "'" << std::endl;
        }
    }

    //!\brief Construct Burrows-Wheeler Transform (BWT) of string from source file array index i.
    // TODO: bug when using "std::array<source_t, 2>::size_type" instead of size_t
    // see lcp_construct_test
    bool construct_bwt(size_t const i, sdsl::int_vector<8> & bwt, bool const reverse=false)
    {
        fs::path log_file = source_path / fs::path{"log_file.txt"};
        fs::path src_file = (reverse ? srcs[i].filename_rev : srcs[i].filename);
        assert(fs::exists(source_path / src_file));
        tMSS f_file_map	= { {sdsl::conf::KEY_TEXT, (source_path / src_file).string()},
                {sdsl::conf::KEY_SA, (source_path / src_file).string()}};
        //sdsl::cache_config config{/*f_delete_files*/ true, /*f_dir*/ tmp_dir.string(), /*f_id*/ "", f_file_map};

        std::ofstream(log_file) << "call construct_bwt" << std::endl;

        // Prepare Input
        bool success;
        test_config = sdsl::cache_config(false, tmp_dir, std::to_string(sdsl::util::pid()));
        uint8_t num_bytes = 1;
        sdsl::int_vector<8> text;
        std::string test_file = (source_path / src_file).string();

        success = sdsl::load_vector_from_file(text, test_file, num_bytes);
        std::ofstream(log_file, std::ios::app) << "load_vector_from_file success: " << success << std::endl;
        success = sdsl::contains_no_zero_symbol(text, test_file);
        std::ofstream(log_file, std::ios::app) << "contains_no_zero_symbol yes: " << success << std::endl;
        sdsl::append_zero_symbol(text);
        success = store_to_cache(text, sdsl::conf::KEY_TEXT, test_config);
        std::ofstream(log_file, std::ios::app) << "store_to_cache yes: " << success << std::endl;

        // Construct SA
        sdsl::int_vector<> sa(text.size(), 0, sdsl::bits::hi(text.size())+1);
        sdsl::algorithm::calculate_sa((const unsigned char*)text.data(), text.size(), sa);
        success = sdsl::store_to_cache(sa, sdsl::conf::KEY_SA, test_config);

        std::ofstream(log_file, std::ios::app) << "store_to_cache successfull: " << success << std::endl;
        std::ofstream(log_file, std::ios::app) << "sa before construct_bwt call: " << std::endl;
        for (auto it = sa.begin(); it != sa.end(); ++it)
            std::ofstream(log_file, std::ios::app) << *it;

        sdsl::construct_bwt<8>(test_config);
        std::ofstream(log_file, std::ios::app) << "sa after construct_bwt call: " << std::endl;
        for (auto it = sa.begin(); it != sa.end(); ++it)
            std::ofstream(log_file, std::ios::app) << *it;
        std::ofstream(log_file, std::ios::app) << std::endl;
        // content written to config.file_map[key] = file_name;
        std::string file_bwt = test_config.file_map[sdsl::conf::KEY_BWT];
        std::ofstream(log_file, std::ios::app) << "bwt array written to " << file_bwt << std::endl;

        std::ofstream(log_file, std::ios::app) << "load bwt from cached file ..." << std::endl;
        success = sdsl::load_from_cache<sdsl::int_vector<8>>(bwt, sdsl::conf::KEY_BWT, test_config);
        std::ofstream(log_file, std::ios::app) << "loading was successfull: " << success << std::endl;

        for (auto it = bwt.begin(); it != bwt.end(); ++it)
            std::ofstream(log_file, std::ios::app) << *it;
        std::ofstream(log_file, std::ios::app) << std::endl;


        return success;
    }


    //!\brief Default constructor.
    constexpr MS() = default;

    //!\brief Copy constructor.
    constexpr MS(MS const &) = default;

    //!\brief Copy construction via assignment.
    constexpr MS & operator=(MS const &) = default;

    //!\brief Move constructor.
    constexpr MS (MS &&) = default;

    //!\brief Move assignment.
    constexpr MS & operator=(MS &&) = default;

    // init with strings and write to file
    constexpr MS(container_t & _s, container_t & _t, unsigned short int const _tau=1)
    {
        srcs[STRING_S].sequence = &_s;
        srcs[STRING_T].sequence = &_t;
        write_files();
    }

    // _filex has to be the absolute path to the sequence files. If not existent
    // an assertion will be thrown.
    // base_path folder containing source files file1 and file2
    constexpr MS(fs::path _source_path, fs::path _file1, fs::path _file2)
    {
        assert(fs::exists(_source_path / _file1) && fs::exists(_source_path / _file1));
        source_path = _source_path;
        srcs[0].filename = _file1;
        srcs[1].filename = _file2;
    }

    //!\brief Use default deconstructor.
    ~MS() = default;

    constexpr value_t select(index_t const i)
    {
        assert(i >= 0);
        return ss(i);
    }

    constexpr value_t operator[](index_t const i) const
    {
        assert(i >= -1);
        if (i == -1)
            return 1;
        return select(i) - 2*i;
    }

     /*
      * Definition "unidirectional matching statistics":
      * Given two strings s and t and a threshold tau > 0,the unidirectional matching
      * statistics MS(t,s,tau) of t with respect to s is a vector of length |t| that
      * stores at index i in [0..|t| − 1] the length of the longest prefix of
      * t[i..|t| − 1] that occurs at least tau times in s.
      * t = AACT, s = AACG, MS = 3210, ms = 000111
      */
      void compute()
      {
         // build suffix tree (ST) from string s and its reverse s_rev
         // TODO: named attribute not possible like reverse=true
         cst_t cst_s, cst_s_rev;
         construct_cst(0, cst_s);
         construct_cst(0, cst_s_rev, true);

         // construct Burrows-Wheeler transform (BWT) from string s
         sdsl::int_vector<8> bwt_s, bwt_s_rev;
         construct_bwt(0, bwt_s);
         construct_bwt(0, bwt_s_rev);

        ms.resize(2*srcs[1].sequence->size());
        // set last bit explicitly to 0 for the case only the first 2|t|-1 bits
        // will be set (see page 4 of [1])
        ms[2*srcs[1].sequence->size()-1] = 0;
        index_t tmp = -1; // MS[-1]
        index_t j = 0;
        // 1st pass: compute consecutive 1s for ms

        sdsl::bit_vector runs = sdsl::bit_vector(srcs[1].sequence->size()-1, 0);
        for (index_t i = 0; i < srcs[1].sequence->size(); ++i)
        {

            if (j >= ms.size()){
                std::cout << "error, j is exceeding allocated bit_vector size with j = " << j << std::endl;
                break;
            }
        }
    }
};

//#endif

//}  // namespace seqan3
