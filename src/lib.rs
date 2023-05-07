/* MIT License
 *
 * Copyright (c) 2023 Andrew Smith
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use, copy,
 * modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */



use std::process;
use std::fs::File;
use std::io::{prelude::*,BufReader};

use serde::{Serialize, Deserialize};
use msite::MSite;

use statrs::distribution::Normal;
use statrs::distribution::ContinuousCDF;


// ADS: (TODO) test this!
fn wilson_ci_for_binomial(
    alpha: f64,
    n: u64,
    p_hat: f64,
    lower: &mut f64,
    upper: &mut f64,
) {
    let g = Normal::new(0.0, 1.0).unwrap();
    let z: f64 = g.inverse_cdf(1.0 - alpha/2.0);

    let n = n as f64;
    let zz = z*z;
    let denom: f64 = 1.0 + zz/n;
    let first_term: f64 = (p_hat + zz/(2.0*n))/denom;
    let discriminant: f64 = p_hat*(1.0 - p_hat)/n + zz/(4.0*n*n);
    let second_term: f64 = z*discriminant.sqrt()/denom;

    *lower = first_term - second_term;
    if *lower < 0.0 {
        *lower = 0.0;
    }
    *upper = first_term + second_term;
    if *upper > 1.0 {
        *upper = 1.0;
    }
}

#[derive(Default,Serialize,Deserialize)]
pub struct LevelsCounter {
    pub total_sites: u64,
    pub sites_covered: u64,
    pub total_c: u64,
    pub total_t: u64,
    pub max_depth: u64,
    pub mutations: u64,
    pub called_meth: u64,
    pub called_unmeth: u64,
    pub mean_agg: f64,

    // derived values
    pub coverage: u64,
    pub sites_covered_fraction: f64,
    pub mean_depth: f64,
    pub mean_depth_covered: f64,
    pub mean_meth: f64,
    pub mean_meth_weighted: f64,
    pub fractional_meth: f64,
}



impl LevelsCounter {
    pub fn update(&mut self, s: &MSite) {
        if s.is_mutated() {
            self.mutations += 1;
        }
        else if s.n_reads > 0 {
            self.sites_covered += 1;
            self.max_depth = std::cmp::max(self.max_depth, s.n_reads);
            self.total_c += s.n_meth();
            self.total_t += s.n_reads - s.n_meth();
            self.mean_agg += s.meth;
            let mut lower: f64 = 0.0;
            let mut upper: f64 = 0.0;
            // ADS: (below) replace this with some way to avoid all
            // the re-evaluation of the Gaussian stuff.
            wilson_ci_for_binomial(LevelsCounter::ALPHA, s.n_reads, s.meth,
                                   &mut lower, &mut upper);
            if lower > 0.5 {
                self.called_meth += 1;
            }
            if upper < 0.5 {
                self.called_unmeth += 1;
            }
        }
        self.total_sites += 1;
    }
    pub fn get_coverage(&self) -> u64 {
        self.total_c + self.total_t
    }
    pub fn get_total_called(&self) -> u64 {
        self.called_meth + self.called_unmeth
    }
    pub fn get_mean_meth_weighted(&self) -> f64 {
        (self.total_c as f64)/(self.get_coverage() as f64)
    }
    pub fn get_fractional_meth(&self) -> f64 {
        (self.called_meth as f64)/(self.get_total_called() as f64)
    }
    pub fn get_mean_meth(&self) -> f64 {
        self.mean_agg/(self.sites_covered as f64)
    }

    const ALPHA: f64 = 0.05;

    pub fn set_derived_values(&mut self) {
        self.coverage = self.get_coverage();
        self.sites_covered_fraction =
            (self.sites_covered as f64)/(self.total_sites as f64);
        self.mean_depth = (self.coverage as f64)/(self.total_sites as f64);
        self.mean_depth_covered =
            (self.coverage as f64)/(self.sites_covered as f64);

        if self.sites_covered != 0 {
            self.mean_meth = self.get_mean_meth();
            self.mean_meth_weighted = self.get_mean_meth_weighted();
        }
        else {
            self.mean_meth = 0.0;
            self.mean_meth_weighted = 0.0;
        }

        if self.get_total_called() > 0 {
            self.fractional_meth = self.get_fractional_meth();
        }
        else {
            self.fractional_meth = 0.0;
        }
    }
}


impl std::fmt::Display for LevelsCounter {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let yaml = serde_yaml::to_string(&self).unwrap();
        write!(f, "{yaml}")
    }
}


// ADS: I'm using this struct because it makes it more convenient to
// have the desired indentation
#[derive(Default,Serialize,Deserialize)]
struct LC {
    cytosine: LevelsCounter,
    cpg: LevelsCounter,
    cpg_symmetric: LevelsCounter,
    chh: LevelsCounter,
    ccg: LevelsCounter,
    cxg: LevelsCounter,
}

pub fn run_mlevels(
    verbose: bool,
    input: &String,
    output: &String,
) {

    // setup the input file
    let in_file = File::open(input).unwrap_or_else(|err| {
        eprintln!("{err}");
        process::exit(1);
    });
    let in_file = BufReader::new(in_file);

    // setup the output stream
    let mut out = File::create(output).unwrap_or_else(|err| {
        eprintln!("{err}");
        process::exit(1);
    });

    // includes all counters
    let mut lc: LC = Default::default();

    let mut prev_site = MSite::new();
    let mut prev_is_cpg = false;

    // iterate over lines in the counts file
    for line in in_file.lines() {
        let line = line.unwrap();

        // make the current line into a site
        let site = MSite::build(&line).unwrap_or_else(|_err| {
            eprintln!("failed parsing site: {}", line);
            process::exit(1);
        });

        if site.chrom != prev_site.chrom && verbose {
            eprintln!("PROCESSING:\t{}",
                      std::str::from_utf8(&site.chrom).unwrap());
        }

        // do this first, because the "add" for CpG can invalidate it
        lc.cytosine.update(&site);

        if site.is_cpg() {
            lc.cpg.update(&site);
            if prev_is_cpg && prev_site.is_mate_of(&site) {
                prev_site.add(&site);
                lc.cpg_symmetric.update(&prev_site);
                prev_is_cpg = false;
            }
            else {
                prev_is_cpg = true;
            }
        }
        else if site.is_chh() {
            lc.chh.update(&site);
        }
        else if site.is_ccg() {
            lc.ccg.update(&site);
        }
        else if site.is_cxg() {
            lc.cxg.update(&site);
        }
        else {
            eprintln!("bad site type: {}", line);
            process::exit(1);
        }

        prev_site = site;
    }

    lc.cytosine.set_derived_values();
    lc.cpg.set_derived_values();
    lc.cpg_symmetric.set_derived_values();
    lc.chh.set_derived_values();
    lc.ccg.set_derived_values();
    lc.cxg.set_derived_values();

    write!(out, "{}", serde_yaml::to_string(&lc).unwrap()).unwrap();
}
