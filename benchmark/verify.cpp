/*******************************************************************************
 * benchmark/verify.cpp
 *
 * R-MAT graph generator test
 *
 * Copyright (C) 2019 Lorenz Hübschle-Schneider <lorenz@4z2.de>
 ******************************************************************************/

#include <benchmark/util.hpp>
#include <rmat/degree_dist.hpp>
#include <rmat/generators/select.hpp>
#include <rmat/graph_generator.hpp>
#include <rmat/parallel_do.hpp>
#include <rmat/rmat.hpp>
#include <rmat/timer.hpp>

#include <tlx/cmdline_parser.hpp>
#include <tlx/logger.hpp>

#include <atomic>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <unistd.h>
#include <vector>

static constexpr bool debug = true;

struct arguments {
    size_t num_edges;
    int log_nodes, depth;
    double par_a, par_b, par_c;
    size_t block_size;
};

template <typename RNG, typename rmat_t>
void run(const arguments &args, size_t seed,
         const std::string &outfn) {
    LOG << "";
    LOG << "R-MAT generation using " << rmat_t::name << " method";

    RNG gen(seed);
    rmat_t r(gen, args.log_nodes, args.par_a, args.par_b, args.par_c);
    r.init(args.depth);

    const bool debug = true;
    rmat::timer timer;
    graph_generator<rmat_t, RNG, true> graphgen(r, args.block_size);
    graphgen.get_edges(args.num_edges, seed, debug);

    double duration = timer.get();

    auto &deg_stats = graphgen.get_stats();
    std::string outfn_template = outfn + '_' +  rmat_t::name + "_XXXXXX";
    char* outfn_gen = new char[outfn_template.length() + 1];
    strncpy(outfn_gen, outfn_template.c_str(), outfn_template.length() + 1);
    int fd = mkstemp(outfn_gen);
    close(fd); //  we'll open it again
    LOG1 << "writing degree stats to: " << outfn_gen;
    deg_stats.write_histogram(outfn_gen);
    delete[] outfn_gen;

    LOG << "RESULT type=sampling"
        << " total=" << duration
        << " edges=" << args.num_edges
        << " logn=" << args.log_nodes
        << " scramble=" << rmat_t::scramble_ids
        << " method=" << rmat_t::name
        << " paths=" << r.table_size()
        << " seed=" << seed
        << " threads=" << rmat::get_num_threads();

    LOG << "RMAT: " << std::fixed << std::setprecision(3) << duration;
}

int main(int argc, char *argv[]) {
    tlx::CmdlineParser clp;

    size_t edges = 10'000'000, seed = 0, num_threads = 0, block_size = 0;
    int log_n = 20, depth = 9;
    double par_a = 0.57, par_b = 0.19, par_c = 0.19;
    bool scramble = false;
    std::string verify_fn;
    clp.add_bool('p', "scramble", scramble, "permute node IDs (scramble)");
    clp.add_int('n', "logn", log_n, "log2 of number of nodes");
    clp.add_size_t('m', "edges", edges, "number of edges to generate");
    clp.add_double('a', "par_a", par_a, "R-MAT parameter a");
    clp.add_double('b', "par_b", par_b, "R-MAT parameter b");
    clp.add_double('c', "par_c", par_c, "R-MAT parameter c");
    clp.add_int('d', "depth", depth, "depth of generated paths (variant 2)");
    clp.add_size_t('s', "seed", seed, "random generator seed (0 for random)");
    clp.add_size_t('g', "blocksize", block_size, "block size for parallelisation");
    clp.add_size_t('t', "threads", num_threads, "number of threads");
    clp.add_string('v', "verify", verify_fn, "output degree distribution base name");

    if (!clp.process(argc, argv))
        return -1;

    if (par_a < 0 || par_b < 0 || par_c < 0 || par_a + par_b + par_c > 1) {
        LOG1 << "Condition a + b + c < 1, a > 0, b > 0, c > 0 violated";
        return -1;
    }
    if (log_n > 62) {
        LOG1 << "More than 2^62 nodes are not supported";
        return -1;
    }
    if (seed == 0) {
        seed = std::random_device{}();
    }
    if (num_threads == 0) {
        num_threads = std::thread::hardware_concurrency();
    }
    if (block_size == 0) {
        block_size = size_t{1} << 16;
    }
    if (verify_fn == "") {
        verify_fn = "/tmp/rmat";
    }
    clp.print_result();

    rmat::init_threads(num_threads);

    using RNG = rmat::generators::select_t;
    LOG << "Selected " << RNG::name << " generator";

    const arguments args { edges, log_n, depth,
                           par_a, par_b, par_c, block_size};

    if (scramble)
        run<RNG, rmat::rmat<true>>(args, seed, verify_fn);
    else
        run<RNG, rmat::rmat<false>>(args, seed, verify_fn);

    rmat::release_threads();
    return 0;
}
