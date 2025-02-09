#!/usr/bin/env python

"""
-------------------------------------------------------------------------------------------------------------------------------------------
MIT License

Copyright (c) 2024 George P. Tiley

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Contact: george.tiley@gmail.com
-------------------------------------------------------------------------------------------------------------------------------------------

estploidy uses a command-line interface to interact with the estploidy package. A successful installation should allow you to access all functions and their help with:

estploidy --help
"""

import click
import time
import datetime
import logging

#--------------------------------#
# CLI ENTRY POINT
#--------------------------------#

@click.group()
def cli():
    current_time = datetime.now()
    time_string = current_time.strftime('%Y-%m-%d-%H-%M-%S')
    logfile = time_string + 'estploidy.log'
    print(f'Start estploidy at {current_time}')
    logging.basicConfig(
       filename = logfile,
       format='%(asctime)s - %(levelname)s - %(message)s',
       filemode = 'w' 
    )
    print(f'Created logfile {logfile}')



@cli.command(context_settings={'help_option_names': ['-h','--help']})

@click.argument('sample_sheet',type=str)
@click.option('-v', '--vcf_file', type=str, default='dummy.vcf',
              help = 'name of the input vcf file. should not be compressed.'
)
@click.option('-i', '--imputation_method', type=str, default='drop',
              help = 'decide how to impute missing data if at all. options are: drop, mean, and popmean'
)
@click.option('-o', '--output_file', type=str, default='samples.vcf',
              help = 'name of the output text file. will be a matrix of allele frequencies'
)

def individual_frequencies(sample_sheet, vcf_file, imputation_method, output_file):
    from estploidy.utils import map_individuals
    from estploidy.calculate_frequencies.calculate_frequencies import get_ind_freqs

    start_time = time.process_time()
    logging.info(f'Begin at {start_time}')
    logging.info(f'Checking all individuals in {sample_sheet} are present in {vcf}')
    ind_map = map_individuals(sample_sheet, vcf_file)
    logging.info(f'Calculating individual allele frequencies from {vcf}')
    get_ind_freqs(ind_map, vcf_file, output_file)
    logging.info(f'Matrix of allele frequencies written out to {output_file}')
    end_time = time.process_time()
    logging.info(f'End at {end_time}')
    compute_time = (end_time - start_time) / 60
    logging.info(f'Total compute time was {compute_time} minutes')

##----------------
## Calculate population-level allele frequencies for fst and genotype-environment association analyses
##----------------
#@click.argument('sample_sheet',type=str)
#@click.option('-v', '--vcf_file', type=str, default='dummy.vcf',
#              help = 'name of the input vcf file. should not be compressed.'
#)
#@click.option('-o', '--output_file', type=str, default='samples.vcf',
#              help = 'name of the output text file. will be a matrix of allele frequencies'
#)
#def population_frequencies(sample_sheet, vcf_file, output_file):
#    pass
